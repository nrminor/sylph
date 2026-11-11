use crate::cmdline::*;
use std::path::Path;
use std::io::prelude::*;
use std::io;
use std::io::BufWriter;
use fxhash::FxHashMap;
use crate::constants::*;
use crate::sketch::*;
use crate::twostage_db::TwoStageDb;
use crate::types::*;
use crate::profile::ProfileError;
use crate::profile_core::{self, Parameters};
use log::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::BufReader;
use std::sync::Mutex;

impl From<&ContainArgs> for Parameters {
    fn from(args: &ContainArgs) -> Self {
        Self {
            minimum_ani: args.minimum_ani,
            pseudotax: args.pseudotax,
            min_number_kmers: args.min_number_kmers,
            min_contain: args.min_contain,
            min_count_correct: args.min_count_correct,
            ratio: args.ratio,
            mme: args.mme,
            nb: args.nb,
            mle: args.mle,
            no_adj: args.no_adj,
            no_ci: args.no_ci,
            mean_coverage: args.mean_coverage,
            seq_id: args.seq_id,
            estimate_unknown: args.estimate_unknown,
            estimate_read_counts: args.estimate_read_counts,
            redundant_ani: args.redundant_ani,
            log_reassignments: args.log_reassignments,
        }
    }
}

fn calculation_failure(error: ProfileError) -> ! {
    log::error!("Profiling calculation failed: {}", error);
    std::process::exit(1)
}

fn print_ani_result(ani_result: &AniResult, pseudotax: bool, writer: &mut Box<dyn Write + Send>) {
    let print_final_ani = format!("{:.2}", f64::min(ani_result.final_est_ani * 100., 100.));
    let lambda_print;
    if let AdjustStatus::Lambda(lambda) = ani_result.lambda {
        lambda_print = format!("{:.3}", lambda);
    } else if ani_result.lambda == AdjustStatus::High {
        lambda_print = format!("HIGH");
    } else {
        lambda_print = format!("LOW");
    }
    let low_ani = ani_result.ani_ci.0;
    let high_ani = ani_result.ani_ci.1;
    let low_lambda = ani_result.lambda_ci.0;
    let high_lambda = ani_result.lambda_ci.1;

    let ci_ani;
    if low_ani.is_none() || high_ani.is_none() {
        ci_ani = "NA-NA".to_string();
    } else {
        ci_ani = format!(
            "{:.2}-{:.2}",
            low_ani.unwrap() * 100.,
            high_ani.unwrap() * 100.
        );
    }

    let ci_lambda;
    if low_lambda.is_none() || high_lambda.is_none() {
        ci_lambda = "NA-NA".to_string();
    } else {
        ci_lambda = format!("{:.2}-{:.2}", low_lambda.unwrap(), high_lambda.unwrap());
    }


    //"Sample_file\tQuery_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\tEff_cov\tANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tNaive_ANI\tContig_name",

    if !pseudotax{
        writeln!(writer, 
            "{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{:.0}\t{:.3}\t{}/{}\t{:.2}\t{}",
            ani_result.seq_name,
            ani_result.gn_name,
            print_final_ani,
            ani_result.final_est_cov,
            ci_ani,
            lambda_print,
            ci_lambda,
            ani_result.median_cov,
            ani_result.mean_cov,
            ani_result.containment_index.0,
            ani_result.containment_index.1,
            ani_result.naive_ani * 100.,
            ani_result.contig_name,
        ).expect("Error writing to file");
    }
    else{
        writeln!(writer,
            "{}\t{}\t{:.4}\t{:.4}\t{}\t{:.3}\t{}\t{}\t{}\t{:.0}\t{:.3}\t{}/{}\t{:.2}\t{}\t{}",
            ani_result.seq_name,
            ani_result.gn_name,
            ani_result.rel_abund.unwrap(),
            ani_result.seq_abund.unwrap(),
            print_final_ani,
            ani_result.final_est_cov,
            ci_ani,
            lambda_print,
            ci_lambda,
            ani_result.median_cov,
            ani_result.mean_cov,
            ani_result.containment_index.0,
            ani_result.containment_index.1,
            ani_result.naive_ani * 100.,
            ani_result.kmers_lost.unwrap(),
            ani_result.contig_name,
        ).expect("Error writing to file");

    }
}

fn get_chunks(indices: &Vec<usize>, steps: usize) -> Vec<Vec<usize>>{
    let mut start = 0;
    let mut end = steps;
    let len = indices.len();
    let mut return_chunks = vec![];

    while start < len {
        if end > len {
            end = len;
        }

        let chunk: Vec<usize> = (start..end).collect();
        start = end;
        end += steps;
        return_chunks.push(chunk);
    }
    return_chunks
}

/// Two-stage stage 1 + 2 against one `.syl2db`: screen `sequence_sketch`
/// against the database's pooled `screen_index` (a single inverted pass over
/// the sample at the sparse, per-database `screen_c`), then decode and return
/// *dense* sketches for only the genomes that pass. The returned set is one
/// database's contribution to the combined active-genome list the (expensive)
/// dense profiling pass runs against.
fn compute_dense_survivors(
    args: &ContainArgs,
    db: &TwoStageDb,
    sequence_sketch: &SequencesSketch,
) -> Vec<GenomeSketch> {
    let screen_index = &db.screen_index;
    // Stage 1: cheap, permissive screen (query-like settings, no CIs).
    let dense_args = Parameters::from(args);
    let mut screen_args = dense_args.clone();
    screen_args.pseudotax = false;
    screen_args.minimum_ani = Some(args.screen_ani);
    screen_args.no_ci = true;
    // The screen scores against sketches sub-sampled to `db.screen_c`, so the
    // dense `-M/--min-number-kmers` floor would over-reject here: a genome has
    // only ~length/screen_c sparse k-mers, so the dense floor M corresponds to a
    // length of M*screen_c, silently dropping smaller genomes (viruses,
    // plasmids, short contigs) that single-stage profiling would report. Scale
    // the floor to the screen resolution so a genome that could clear the dense
    // floor also clears the screen; genomes truly below the floor are still
    // rejected at the dense stage. `db.screen_c` is the *stored* (possibly
    // adaptively-densified, see `write_two_stage_db`) rate, so this scaling
    // automatically becomes less aggressive too when a small genome forced a
    // denser effective rate. `--min-contain` is deliberately deferred to the
    // decoded dense sketch below: applying an absolute dense evidence floor to
    // this much sparser screen could discard a valid dense candidate.
    screen_args.min_number_kmers = args.min_number_kmers * db.c as f64 / db.screen_c as f64;

    // Stage 1 = "Path B": one inverted pass over the sample produces, per genome,
    // the same matched-coverage multiset the per-genome `get_stats` loop would
    // collect; feeding it to the same `finalize_stats` + `min_number_kmers`
    // checks reproduces the survivor set exactly, in O(sample) rather than
    // O(reference) work. (See experiments/7_mphf_screen_again.)
    let hits: Vec<(u32, Vec<u32>)> = screen_index
        .gather_hits(sequence_sketch)
        .into_iter()
        .collect();
    // Optional per-survivor dump: (genome, matched k-mers, total screen k-mers,
    // naive ANI, adjusted ANI, median coverage).
    let dump: Mutex<Vec<(String, usize, usize, f64, f64, f64)>> = Mutex::new(vec![]);
    let mut survivors: Vec<usize> = hits
        .into_par_iter()
        .filter_map(|(g, covs)| {
            let g = g as usize;
            let n_kmers = screen_index.sparse_count[g] as usize;
            // Mirror get_stats: reject genomes below the (scaled) k-mer floor.
            if (n_kmers as f64) < screen_args.min_number_kmers {
                return None;
            }
            let contain_count = covs.len();
            // finalize_stats applies the screen min-ANI gate (returns None below it).
            let fin = profile_core::finalize_stats(&screen_args, db.k, n_kmers, contain_count, covs, None)
                .unwrap_or_else(|error| calculation_failure(error))?;
            if args.screen_dump.is_some() {
                dump.lock().unwrap().push((
                    db.genome_file_name(g as u32).to_string(),
                    contain_count,
                    n_kmers,
                    fin.naive_ani * 100.,
                    fin.final_est_ani * 100.,
                    fin.median_cov,
                ));
            }
            Some(g)
        })
        .collect();
    survivors.sort_unstable();
    log::info!(
        "{}: stage-1 screen (c={}, min-ANI {}) kept {} / {} candidate genomes",
        sequence_sketch.file_name, db.screen_c, args.screen_ani, survivors.len(), screen_index.num_genomes()
    );
    if let Some(path) = &args.screen_dump {
        let mut f = BufWriter::new(File::create(path).expect("could not create --screen-dump file"));
        writeln!(f, "Genome_file\tscreen_matched_kmers\tscreen_total_kmers\tnaive_ani\tscreen_adjusted_ani\tscreen_median_cov").unwrap();
        for (g, m, t, na, ea, mc) in dump.into_inner().unwrap() {
            writeln!(f, "{}\t{}\t{}\t{:.3}\t{:.3}\t{}", g, m, t, na, ea, mc).unwrap();
        }
        log::info!("Wrote stage-1 screen dump to {}", path);
    }

    // Stage 2: decode each survivor's dense block into a short-lived sketch,
    // run the (prefetch-friendly) pass-1 profiling on it, and keep it only if
    // it passes -- so the discarded majority is freed immediately and never
    // cached, keeping peak RAM proportional to the genomes that survive rather
    // than to the (much larger) screen-survivor set.
    let dense: Mutex<Vec<GenomeSketch>> = Mutex::new(vec![]);
    survivors.par_iter().for_each(|i| {
        match db.decode_dense(*i as u32) {
            Ok(g) => {
                // Keep only genomes that pass pass-1 profiling; drop (free) the rest.
                if profile_core::get_stats(&dense_args, &g, sequence_sketch, None, false)
                    .unwrap_or_else(|error| calculation_failure(error)).is_some() {
                    dense.lock().unwrap().push(g);
                }
            }
            Err(e) => warn!("Could not decode dense block for genome index {}: {}", i, e),
        }
    });
    let dense = dense.into_inner().unwrap();
    log::info!(
        "{}: stage-2 dense profiling (c={}) against {} genomes",
        sequence_sketch.file_name, db.c, dense.len()
    );
    dense
}

pub fn contain(mut args: ContainArgs, pseudotax_in: bool) {

    if pseudotax_in{
        args.pseudotax = true;
    }

    let level;
    if args.trace {
        level = log::LevelFilter::Trace;
    } else if args.debug {
        level = log::LevelFilter::Debug;
    }
    else{
        level = log::LevelFilter::Info;
    }
    
    simple_logger::SimpleLogger::new()
        .with_level(level)
        .init()
        .unwrap();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    let out_writer = match args.out_file_name {
        Some(ref x) => {
            let path = Path::new(&x);
            Box::new(BufWriter::new(File::create(&path).unwrap())) as Box<dyn Write + Send>
        }
        None => Box::new(BufWriter::new(io::stdout())) as Box<dyn Write + Send>,
    };

    if args.estimate_read_counts{
        args.estimate_unknown = true;
        log::info!("--estimate-read-counts detected, also enabling -u. Sequence_abundance column will be set to estimated read counts, not abundance. This is still experimental.");
    }
    let parameters = Parameters::from(&args);

    log::info!("Obtaining sketches...");
    let mut genome_sketch_files = vec![];
    let mut genome_files = vec![];
    let mut two_stage_db_files = vec![];
    let mut read_sketch_files = vec![];
    let mut read_files = vec![];

    let mut all_files = args.files.clone();
    all_files.extend(args.databases.clone());

    if let Some(ref newline_file) = args.file_list{
        let file = File::open(newline_file).unwrap();
        let reader = BufReader::new(file);
        for line in reader.lines() {
            all_files.push(line.unwrap());
        }

    }


    for file in all_files.iter(){

        let mut genome_sketch_good_suffix = false;
        for suff in QUERY_FILE_SUFFIX_VALID{
            if file.ends_with(suff){
                genome_sketch_good_suffix = true;
                break
            }
        }

        let mut sample_sketch_good_suffix = false;
        for suff in SAMPLE_FILE_SUFFIX_VALID{
            if file.ends_with(suff){
                sample_sketch_good_suffix = true;
                break
            }
        }

        if file.ends_with(TWO_STAGE_DB_SUFFIX){
            two_stage_db_files.push(file);
        } else if genome_sketch_good_suffix{
            genome_sketch_files.push(file);
        } else if sample_sketch_good_suffix{
            read_sketch_files.push(file);
        } else if is_fasta(&file) {
            genome_files.push(file);
        } else if is_fastq(&file) {
            read_files.push(vec![file]);
        } else {
            warn!(
                "{} file extension is not a sketch or a fasta/fastq file.",
                &file
            );
        }
    }

    if args.first_pair.len() != args.second_pair.len() {
        error!("Different number of paired sequences (-1, -2) for sketching. Exiting.");
        std::process::exit(1);
    }

    // zip together the first and second pair files, push them to read_files
    for (first, second) in args.first_pair.iter().zip(args.second_pair.iter()) {
        read_files.push(vec![first,second]);
    }

    for read in args.reads.iter() {
        read_files.push(vec![read]);
    }

    if genome_sketch_files.is_empty() && genome_files.is_empty() && two_stage_db_files.is_empty(){
        log::error!("No genome files found; see sylph query/profile -h for help. Exiting");
        std::process::exit(1);
    }

    if read_sketch_files.is_empty() && read_files.is_empty(){
        log::error!("No read files found; see sylph query/profile -h for help. Exiting");
        std::process::exit(1);
    }

    // Any number of .syl2db files can be combined with each other and with
    // plain .syldb/fasta genomes in one call: each .syl2db is screened and
    // decoded independently per sample (see the per-sample loop below), and
    // its survivors are unioned with the (always fully-loaded, unscreened)
    // plain genomes before profiling/reassignment.
    let two_stage_dbs: Vec<TwoStageDb> = two_stage_db_files.iter().map(|f| {
        log::info!("Opening two-stage database {} (loading stage-1 sparse index)...", f);
        crate::twostage_db::open_file(f)
            .unwrap_or_else(|e| panic!("{} is not a valid two-stage database: {}", f, e))
    }).collect();
    for db in &two_stage_dbs {
        if db.is_empty() {
            log::error!("Two-stage database contains no genomes. Exiting");
            std::process::exit(1);
        }
    }

    let genome_sketches = get_genome_sketches(&args, &genome_sketch_files, &genome_files);
    log::info!("Finished obtaining genome sketches.");

    if two_stage_dbs.is_empty() && genome_sketches.is_empty() {
        log::error!("No genome sketches found; see sylph query/profile -h for help. Exiting");
        std::process::exit(1);
    }
    if !genome_sketches.is_empty()
        && genome_sketches.first().unwrap().pseudotax_tracked_nonused_kmers.is_none() && args.pseudotax {
        log::error!("Attempting profiling, but *.syldb was sketched with the --disable-profiling option. Exiting");
        std::process::exit(1);
    }
    // (No equivalent pseudotax-completeness check is needed for two_stage_dbs:
    // db-convert already requires 100% pseudotax coverage at conversion time.)

    // Duplicate genomes loaded across multiple sources (e.g. the same genome
    // present in both a .syl2db and a plain .syldb) collide as identical
    // GenomeSketch entries downstream in reassignment/dereplication, causing
    // nondeterministic output (whichever copy a racing thread inserts first
    // "wins"). Warn loudly rather than silently producing inconsistent runs.
    {
        let mut seen: FxHashMap<(&str, &str), usize> = FxHashMap::default();
        for db in &two_stage_dbs {
            for g in 0..db.len() {
                *seen.entry((db.genome_file_name(g as u32), db.genome_first_contig_name(g as u32))).or_insert(0) += 1;
            }
        }
        for gs in &genome_sketches {
            *seen.entry((gs.file_name.as_str(), gs.first_contig_name.as_str())).or_insert(0) += 1;
        }
        let dups: Vec<&(&str, &str)> = seen.iter().filter(|(_, &n)| n > 1).map(|(k, _)| k).collect();
        if !dups.is_empty() {
            log::warn!(
                "{} genome(s) appear more than once across the combined input sources \
                 (same file name + contig name loaded from two different databases): {:?}. \
                 This can cause nondeterministic reassignment/dereplication output -- consider \
                 removing the duplicate(s) from one of the input databases.",
                dups.len(), dups
            );
        }
    }

    // k must match exactly across every loaded source.
    let mut ks: Vec<usize> = two_stage_dbs.iter().map(|db| db.k).collect();
    if !genome_sketches.is_empty() { ks.push(genome_sketches[0].k); }
    ks.dedup();
    if ks.len() > 1 {
        log::error!("Inconsistent -k across combined genome sources ({:?}). Exiting", ks);
        std::process::exit(1);
    }
    let db_k = ks[0];

    // The sample sketch must not be sparser than any active genome's own rate
    // (the shared calculation rejects mismatches): use the minimum -c across every loaded
    // source (two-stage dense rate and plain -c alike).
    let effective_genome_c = two_stage_dbs.iter().map(|db| db.c)
        .chain(genome_sketches.iter().map(|g| g.c))
        .min()
        .expect("at least one genome source must be present (checked above)");

    log::info!(
        "Loaded {} two-stage database(s) ({} genomes total) and {} plain genome sketch(es).",
        two_stage_dbs.len(),
        two_stage_dbs.iter().map(|d| d.len()).sum::<usize>(),
        genome_sketches.len(),
    );

    let num_raw_read_files = read_files.len();
    let step;
    if let Some(sample_threads) = args.sample_threads{
        if sample_threads > 0{
            step = sample_threads;
        }
        else{
            step = 1;
        }
    }
    else{
        // Historically `profile` (pseudotax) reserved fewer samples-in-flight
        // than `query` via the `threads/3 + 1` floor, leaving headroom for the
        // per-sample reassignment pass. That floor was inert whenever
        // `num_raw_read_files >= threads/3 + 1` (the common case), since with
        // single-threaded-per-file sketching `step` only mattered via
        // `min(step, num_raw_read_files)` in `get_chunks` below -- but it
        // becomes an active bottleneck now that `step` also sizes
        // `threads_per_file` for the multi-threaded sketch pipeline: with few
        // files and many threads, inflating `step` starves each file's
        // pipeline of threads for no benefit (nested rayon parallelism in the
        // reassignment pass doesn't need reserved threads; it work-steals from
        // the same global pool). Match `query`'s formula.
        step = usize::max(1, usize::min(num_raw_read_files, args.threads))
    }

    // Threads available per file's own sketching pipeline, given `step` files
    // are processed concurrently out of the shared `args.threads` budget.
    let threads_per_file = (args.threads / step.max(1)).max(1);

    let read_sketch_files_as_vec = read_sketch_files.clone().into_iter().map(|x| vec![x]).collect::<Vec<Vec<&String>>>();
    read_files.extend(read_sketch_files_as_vec);
    let sequence_index_vec = (0..read_files.len()).collect::<Vec<usize>>();
    let out_writer:Mutex<Box<dyn Write + Send>> = Mutex::new(out_writer);

    let chunks = get_chunks(&sequence_index_vec, step);

    print_header(args.pseudotax,&mut *out_writer.lock().unwrap(), args.estimate_unknown);
    chunks.into_iter().for_each(|chunk| {
        chunk.into_par_iter().for_each(|j|{
            let is_sketch = j >= read_files.len() - read_sketch_files.len();
            let sequence_sketch = get_seq_sketch(&args, &read_files[j], is_sketch, effective_genome_c, db_k, threads_per_file);
            if sequence_sketch.is_some(){
                let first_read_file = read_files[j][0];
                let sequence_sketch = sequence_sketch.unwrap();

                // Screen each two-stage database, then union its survivors with
                // the (always fully-loaded, unscreened) plain genomes into one
                // combined active-genome list for this sample.
                let mut dense_local: Vec<GenomeSketch> = Vec::new();
                for db in &two_stage_dbs {
                    dense_local.extend(compute_dense_survivors(&args, db, &sequence_sketch));
                }
                let active_sketches: Vec<&GenomeSketch> =
                    dense_local.iter().chain(genome_sketches.iter()).collect();

                let stats_vec_seq = profile_core::calculate(&parameters, &active_sketches, &sequence_sketch)
                    .unwrap_or_else(|error| calculation_failure(error));
                if args.pseudotax {
                    log::info!("{} has {} genomes passing profiling threshold. ", &first_read_file, stats_vec_seq.len());
                }
                
                let mut out_writer = out_writer.lock().unwrap();
                for res in stats_vec_seq{
                    print_ani_result(&res, args.pseudotax, &mut *out_writer);
                }
            }
            if read_files[j].len() > 1{
                log::info!("Finished paired sample {}.", &read_files[j][0]);
            }
            else{
                log::info!("Finished sample {}.", &read_files[j][0]);
            }
        });
    });

    log::info!("sylph finished.");
}

fn print_header(pseudotax: bool, writer: &mut Box<dyn Write + Send>, estimate_unknown: bool) {
    if !pseudotax{
        writeln!(writer,
            //"Sample_file\tQuery_file\tAdjusted_ANI\tNaive_ANI\tANI_5-95_percentile\tEff_cov\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tContig_name",
            "Sample_file\tGenome_file\tAdjusted_ANI\tEff_cov\tANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tNaive_ANI\tContig_name",
            ).expect("Error writing to file.");
    }
    else{
        let cov_head;
        if estimate_unknown{
            cov_head = "True_cov";
        }
        else{
            cov_head = "Eff_cov";
        }
        writeln!(writer,
            "Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\t{}\tANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tNaive_ANI\tkmers_reassigned\tContig_name", cov_head
            ).expect("Error writing to file.");
    }
}

fn get_genome_sketches(
    args: &ContainArgs,
    genome_sketch_files: &Vec<&String>,
    genome_files: &Vec<&String>,
) -> Vec<GenomeSketch> {
    let mut lowest_genome_c = None;
    let mut current_k = None;

    let genome_sketches = Mutex::new(vec![]);

    for genome_sketch_file in genome_sketch_files {
        let file = File::open(genome_sketch_file).expect(&format!("The sketch `{}` could not be opened. Exiting", genome_sketch_file));
        let genome_reader = BufReader::with_capacity(10_000_000, file);
        let genome_sketches_vec: Vec<GenomeSketch> = bincode::deserialize_from(genome_reader)
            .expect(&format!(
                "The sketch `{}` is not a valid sketch. Perhaps it is an older, incompatible version ",
                &genome_sketch_file
            ));
        if genome_sketches_vec.is_empty() {
            continue;
        }
        let c = genome_sketches_vec.first().unwrap().c;
        let k = genome_sketches_vec.first().unwrap().k;
        if lowest_genome_c.is_none() {
            lowest_genome_c = Some(c);
        } else if lowest_genome_c.unwrap() < c {
            lowest_genome_c = Some(c);
        }
        if current_k.is_none() {
            current_k = Some(genome_sketches_vec.first().unwrap().k);
        } else if current_k.unwrap() != k {
            error!("Query sketches have inconsistent -k. Exiting.");
            std::process::exit(1);
        }
        genome_sketches.lock().unwrap().extend(genome_sketches_vec);
    }

    genome_files.into_par_iter().for_each(|genome_file|{
        if lowest_genome_c.is_some() && lowest_genome_c.unwrap() < args.c{
            error!("Value of -c for contain is {} -- greater than the smallest value of -c for a genome sketch {}. Continuing without sketching.", args.c, lowest_genome_c.unwrap());
        }
        else if current_k.is_some() && current_k.unwrap() != args.k{
            error!("-k {} is not equal to -k {} found in sketches. Continuing without sketching.", args.k, current_k.unwrap());
        }
        else {
            if args.individual{
            let indiv_gn_sketches = sketch_genome_individual(args.c, args.k, genome_file, args.min_spacing_kmer, args.pseudotax);
                genome_sketches.lock().unwrap().extend(indiv_gn_sketches);

            }
            else{
                let genome_sketch_opt = sketch_genome(args.c, args.k, &genome_file, args.min_spacing_kmer, args.pseudotax);
                if genome_sketch_opt.is_some() {
                    genome_sketches.lock().unwrap().push(genome_sketch_opt.unwrap());
                }
            }
        }
    });

    return genome_sketches.into_inner().unwrap();
}

fn get_seq_sketch(
    args: &ContainArgs,
    read_file: &Vec<&String>,
    is_sketch_file: bool,
    genome_c: usize,
    genome_k: usize,
    threads: usize,
) -> Option<SequencesSketch> {
    if is_sketch_file {
        let read_file = read_file[0];
        let read_sketch_file = read_file;
        let file = File::open(read_sketch_file).expect(&format!(
            "The sketch `{}` could not be opened",
            &read_sketch_file
        ));
        let read_reader = BufReader::with_capacity(10_000_000, file);
        let read_sketch: SequencesSketch = bincode::deserialize_from(read_reader).expect(
            &format!("The sketch `{}` is not a valid sketch. Perhaps it is an older incompatible version ", read_sketch_file),
        );
        if read_sketch.c > genome_c {
            error!("{} is a pre-sketched sample at -c {}; this is sparser than required by the loaded genome sketch(es) (-c {}). A pre-sketched sample can't be re-sketched denser after the fact -- re-sketch it at -c <= {} first, or use a raw fastq input instead. Continuing without sketching.", read_file, read_sketch.c, genome_c, genome_c);
            return None;
        }
        else if read_sketch.c < genome_c{
            info!("{} is a pre-sketched sample at -c {}, denser than required by the loaded genome sketch(es) (-c {}) -- using it as-is.", read_file, read_sketch.c, genome_c);
        }

        return Some(read_sketch);
    } else {
        // Raw reads are sketched fresh here, so unlike the pre-sketched-file
        // branch above there's no unrecoverable case: if the requested -c is
        // sparser than what the loaded genome(s) need, just sketch at the
        // denser rate genome_c instead -- never sparser than requested, and
        // always safe/correct (a sample can't be sparser than any genome it's
        // compared against).
        let effective_c = args.c.min(genome_c);
        if args.c > genome_c{
            info!("{} value of -c for reads is {}; this is sparser than required by the loaded genome sketch(es) (-c {}). Sketching at the denser -c {} instead.", read_file[0], args.c, genome_c, effective_c);
        }
        if genome_k != args.k {
            error!(
                "{} -k {} is not equal to -k {} found in sketches. Continuing without sketching.",
                read_file[0], args.k, genome_k
            );
            return None;
        } else {
            let pipeline_params = crate::parallel_sketch::PipelineParams {
                batch_records: args.sketch_batch_size,
                batch_max_bytes: args.sketch_batch_max_bytes,
                channel_depth: args.sketch_channel_depth,
                num_shards: args.sketch_shards,
            };
            if read_file.len() == 1{
                let read_sketch_opt = if crate::parallel_sketch::should_use_pipeline(args.no_sketch_pipeline, threads, &read_file[0]) {
                    crate::parallel_sketch::sketch_sequences_needle_parallel(&read_file[0], effective_c, args.k, None, false, threads, &pipeline_params)
                } else {
                    sketch_sequences_needle(&read_file[0], effective_c, args.k, None, false)
                };
                return read_sketch_opt;
            }
            else if read_file.len() == 2{
                let read_sketch_opt = if crate::parallel_sketch::should_use_pipeline(args.no_sketch_pipeline, threads, &read_file[0]) {
                    crate::parallel_sketch::sketch_pair_sequences_parallel(&read_file[0], &read_file[1], effective_c, args.k, None, false, DEFAULT_FPR, threads, &pipeline_params)
                } else {
                    sketch_pair_sequences(&read_file[0], &read_file[1], effective_c, args.k, None, false, DEFAULT_FPR)
                };
                return read_sketch_opt;
            }
            else{
                panic!("Internal Error: read_file has length {}. Something went wrong...", read_file.len());
            }
        }
    }
}


#[cfg(all(test, feature = "cli"))]
mod min_number_kmers_tests {
    use super::*;
    use clap::Parser;

    fn parse_contain_args(subcommand: &str, extra_args: &[&str]) -> ContainArgs {
        let mut argv = vec!["sylph", subcommand];
        argv.extend_from_slice(extra_args);
        argv.extend(["database.syldb", "sample.sylsp"]);

        let cli = Cli::try_parse_from(argv).unwrap();
        match cli.mode {
            Mode::Profile(mut args) => {
                args.pseudotax = true;
                args
            }
            Mode::Query(args) => args,
            _ => panic!("expected profile or query arguments"),
        }
    }

    fn synthetic_genome_passes(args: &ContainArgs, genome_kmers: usize, contained: usize) -> bool {
        let genome = GenomeSketch {
            genome_kmers: (0..genome_kmers as u64).collect(),
            pseudotax_tracked_nonused_kmers: Some(vec![]),
            file_name: "genome".to_string(),
            first_contig_name: "contig".to_string(),
            c: 200,
            k: 31,
            gn_size: 1_000,
            min_spacing: 0,
        };
        let mut sequence = SequencesSketch::new("sample".to_string(), 200, 31, false, None, 150.0);
        sequence
            .kmer_counts
            .extend((0..contained as u64).map(|kmer| (kmer, 1)));

        profile_core::get_stats(&Parameters::from(args), &genome, &sequence, None, false)
            .expect("valid synthetic coverage evidence").is_some()
    }

    #[test]
    fn min_number_kmers_defaults_to_10() {
        let args = parse_contain_args("profile", &[]);

        assert_eq!(args.min_number_kmers, 10.0);
    }

    #[test]
    fn explicit_min_number_kmers_overrides_default() {
        let args = parse_contain_args("profile", &["-M", "32.5"]);

        assert_eq!(args.min_number_kmers, 32.5);
    }

    #[test]
    fn min_number_kmers_discards_small_genome_sketches() {
        let mut args = parse_contain_args("profile", &[]);
        args.minimum_ani = Some(0.0);
        args.no_ci = true;

        assert!(!synthetic_genome_passes(&args, 9, 9));
    }

    #[test]
    fn min_contain_defaults_to_7_for_all_genome_sizes() {
        let mut args = parse_contain_args("profile", &[]);
        assert_eq!(args.min_contain, 7);

        args.minimum_ani = Some(0.0);
        args.no_ci = true;

        assert!(!synthetic_genome_passes(&args, 10, 6));
        assert!(synthetic_genome_passes(&args, 10, 7));
        assert!(!synthetic_genome_passes(&args, 500, 6));
        assert!(synthetic_genome_passes(&args, 500, 7));
    }

    #[test]
    fn explicit_min_contain_overrides_default() {
        let mut args = parse_contain_args("profile", &["--min-contain", "3"]);
        assert_eq!(args.min_contain, 3);

        args.minimum_ani = Some(0.0);
        args.no_ci = true;
        assert!(!synthetic_genome_passes(&args, 500, 2));
        assert!(synthetic_genome_passes(&args, 500, 3));
    }
}
