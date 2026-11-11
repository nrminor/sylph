use std::sync::Mutex;

use fxhash::FxHashMap;
use rayon::prelude::*;
use statrs::distribution::{DiscreteCDF, Poisson};

use crate::constants::*;
use crate::inference::*;
use crate::profile::ProfileError;
use crate::types::*;

#[derive(Clone)]
pub(crate) struct Parameters {
    pub(crate) minimum_ani: Option<f64>,
    pub(crate) pseudotax: bool,
    pub(crate) min_number_kmers: f64,
    pub(crate) min_contain: usize,
    pub(crate) min_count_correct: f64,
    pub(crate) ratio: bool,
    pub(crate) mme: bool,
    pub(crate) nb: bool,
    pub(crate) mle: bool,
    pub(crate) no_adj: bool,
    pub(crate) no_ci: bool,
    pub(crate) mean_coverage: bool,
    pub(crate) seq_id: Option<f64>,
    pub(crate) estimate_unknown: bool,
    pub(crate) estimate_read_counts: bool,
    pub(crate) redundant_ani: f64,
    pub(crate) log_reassignments: bool,
}

type ConfidenceInterval = (Option<f64>, Option<f64>);

pub(crate) struct Finalized {
    pub(crate) naive_ani: f64,
    pub(crate) final_est_ani: f64,
    pub(crate) final_est_cov: f64,
    pub(crate) median_cov: f64,
    pub(crate) mean_cov: f64,
    pub(crate) lambda: AdjustStatus,
    pub(crate) ani_ci: ConfidenceInterval,
    pub(crate) lambda_ci: ConfidenceInterval,
}

fn finite(value: f64, what: &'static str) -> Result<f64, ProfileError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ProfileError::NonFiniteEstimate(what))
    }
}

fn checked_sum(covs: &[u32]) -> Result<u32, ProfileError> {
    covs.iter()
        .try_fold(0u32, |sum, &x| sum.checked_add(x))
        .ok_or(ProfileError::ArithmeticOverflow("coverage sum"))
}

fn minimum_ani_fraction(args: &Parameters) -> f64 {
    if let Some(value) = args.minimum_ani {
        value / 100.0
    } else if args.pseudotax {
        MIN_ANI_P_DEF
    } else {
        MIN_ANI_DEF
    }
}

pub(crate) fn calculate<'a>(
    args: &Parameters,
    genomes: &[&'a GenomeSketch],
    sample: &SequencesSketch,
) -> Result<Vec<AniResult<'a>>, ProfileError> {
    let results = Mutex::new(Vec::new());
    genomes
        .par_iter()
        .try_for_each(|genome| -> Result<(), ProfileError> {
            if let Some(result) = get_stats(args, genome, sample, None, args.log_reassignments)? {
                results.lock().expect("result mutex poisoned").push(result);
            }
            Ok(())
        })?;
    let mut results = results.into_inner().expect("result mutex poisoned");
    if args.pseudotax {
        // First-pass coverage correction is unused by reassignment and dereplication.
        log::info!(
            "{} taxonomic profiling; reassigning k-mers for {} genomes...",
            sample.file_name,
            results.len()
        );
        let winners = winner_table(&results, args.log_reassignments);
        let remaining: Vec<_> = results.iter().map(|r| r.genome_sketch).collect();
        let reassigned = Mutex::new(Vec::new());
        remaining
            .into_par_iter()
            .try_for_each(|genome| -> Result<(), ProfileError> {
                if let Some(result) =
                    get_stats(args, genome, sample, Some(&winners), args.log_reassignments)?
                {
                    reassigned
                        .lock()
                        .expect("result mutex poisoned")
                        .push(result);
                }
                Ok(())
            })?;
        results = derep(
            &results,
            reassigned.into_inner().expect("result mutex poisoned"),
            args.redundant_ani,
            sample.k,
        )?;
        if results.is_empty() {
            return Ok(results);
        }
    }
    if !results.is_empty() && args.estimate_unknown {
        let identity = match args.seq_id {
            Some(id) => (id / 100.0).powf(sample.k as f64),
            None => get_kmer_identity(sample)?,
        };
        if !identity.is_finite() || identity <= 0.0 {
            return Err(ProfileError::InsufficientEvidence(
                "read identity correction must be finite and positive",
            ));
        }
        log::debug!(
            "{} has estimated identity {:.3}.",
            sample.file_name,
            identity.powf(1.0 / sample.k as f64) * 100.0
        );
        let denominator = sample.mean_read_length - sample.k as f64 + 1.0;
        if !sample.mean_read_length.is_finite()
            || sample.mean_read_length <= 0.0
            || !denominator.is_finite()
            || denominator <= 0.0
        {
            return Err(ProfileError::InsufficientEvidence(
                "read length correction denominator must be finite and positive",
            ));
        }
        let multiplier = finite(
            sample.mean_read_length / denominator,
            "read length correction",
        )?;
        for result in &mut results {
            result.final_est_cov = finite(
                result.final_est_cov / identity * multiplier,
                "corrected coverage",
            )?;
        }
    }
    if args.pseudotax {
        let bases_explained = if args.estimate_unknown {
            let fraction = estimate_covered_bases(&results, sample)?;
            log::info!(
                "{} has {:.2}% of reads detected in database by profile",
                sample.file_name,
                fraction * 100.0
            );
            fraction
        } else {
            1.0
        };
        let total_cov = finite(
            results.iter().map(|x| x.final_est_cov).sum(),
            "total coverage",
        )?;
        let total_seq_cov = finite(
            results
                .iter()
                .map(|x| x.final_est_cov * x.genome_sketch.gn_size as f64)
                .sum(),
            "total sequence coverage",
        )?;
        if total_cov <= 0.0 || (!args.estimate_read_counts && total_seq_cov <= 0.0) {
            return Err(ProfileError::InsufficientEvidence(
                "abundance denominator must be positive",
            ));
        }
        if args.estimate_read_counts
            && (!sample.mean_read_length.is_finite() || sample.mean_read_length <= 0.0)
        {
            return Err(ProfileError::InsufficientEvidence(
                "read length must be finite and positive",
            ));
        }
        for result in &mut results {
            result.rel_abund = Some(finite(
                result.final_est_cov / total_cov * 100.0,
                "taxonomic abundance",
            )?);
            let seq = if args.estimate_read_counts {
                (result.final_est_cov * result.genome_sketch.gn_size as f64
                    / sample.mean_read_length
                    * bases_explained)
                    .round()
            } else {
                result.final_est_cov * result.genome_sketch.gn_size as f64 / total_seq_cov
                    * 100.0
                    * bases_explained
            };
            result.seq_abund = Some(finite(seq, "sequence abundance")?);
        }
        results.sort_by(|x, y| {
            y.rel_abund
                .unwrap()
                .partial_cmp(&x.rel_abund.unwrap())
                .expect("finite abundance")
        });
    } else {
        results.sort_by(|x, y| {
            y.final_est_ani
                .partial_cmp(&x.final_est_ani)
                .expect("finite ANI")
        });
    }
    Ok(results)
}

fn derep<'a>(
    old: &[AniResult],
    new: Vec<AniResult<'a>>,
    ani: f64,
    k: usize,
) -> Result<Vec<AniResult<'a>>, ProfileError> {
    let mut lookup = FxHashMap::default();
    for result in old {
        lookup.insert(result.genome_sketch, result);
    }
    let threshold = finite((ani / 100.0).powf(k as f64), "reassignment threshold")?;
    let mut retained = Vec::new();
    for result in new {
        let previous = lookup
            .get(result.genome_sketch)
            .expect("reassigned genome came from first pass");
        let lost = previous
            .containment_index
            .0
            .checked_sub(result.containment_index.0)
            .ok_or(ProfileError::InvalidSketch(
                "reassignment increased containment",
            ))? as f64;
        let limit = threshold * result.containment_index.1 as f64;
        if lost < limit {
            retained.push(result);
        } else {
            log::debug!(
                "genome {} had num k-mers reassigned = {}, threshold was {}, removing.",
                result.gn_name,
                lost,
                limit
            );
        }
    }
    Ok(retained)
}

fn estimate_covered_bases(
    results: &[AniResult],
    sample: &SequencesSketch,
) -> Result<f64, ProfileError> {
    let multiplier = sample.mean_read_length / (sample.mean_read_length - sample.k as f64 + 1.0);
    let covered = finite(
        results
            .iter()
            .map(|r| r.genome_sketch.gn_size as f64 * r.final_est_cov)
            .sum(),
        "covered bases",
    )?;
    let counts = sample
        .kmer_counts
        .values()
        .try_fold(0usize, |sum, &count| sum.checked_add(count as usize))
        .ok_or(ProfileError::ArithmeticOverflow("sample marker count"))?;
    let bases = sample
        .c
        .checked_mul(counts)
        .ok_or(ProfileError::ArithmeticOverflow("sample covered bases"))?;
    let tentative = finite(bases as f64 * multiplier, "sample covered bases")?;
    if tentative == 0.0 {
        return Ok(0.0);
    }
    finite(
        f64::min(covered / tentative, 1.0),
        "fraction of covered bases",
    )
}

fn winner_table<'a>(
    results: &'a Vec<AniResult>,
    log_reassign: bool,
) -> FxHashMap<Kmer, (f64, &'a GenomeSketch, bool)> {
    let mut map = FxHashMap::default();
    for res in results {
        let genome = res.genome_sketch;
        for kmer in genome
            .genome_kmers
            .iter()
            .chain(genome.pseudotax_tracked_nonused_kmers.iter().flatten())
        {
            let winner = map
                .entry(*kmer)
                .or_insert((res.final_est_ani, genome, false));
            if res.final_est_ani > winner.0 {
                *winner = (res.final_est_ani, genome, true);
            }
        }
    }
    if log_reassign {
        log::info!("------------- Logging k-mer reassignments -----------------");
        let mut indices = FxHashMap::default();
        for (i, res) in results.iter().enumerate() {
            log::info!(
                "Index\t{}\t{}\t{}",
                i,
                res.genome_sketch.file_name,
                res.genome_sketch.first_contig_name
            );
            indices.insert(res.genome_sketch, i);
        }
        (0..results.len()).into_par_iter().for_each(|i| {
            let res = &results[i];
            let mut edges = FxHashMap::default();
            for kmer in &res.genome_sketch.genome_kmers {
                let winner = map[kmer].1;
                if winner != res.genome_sketch {
                    *edges.entry((indices[winner], i)).or_insert(0usize) += 1;
                }
            }
            for (key, count) in edges {
                if count > 10 {
                    log::info!("{}->{}\t{}\tkmers reassigned", key.0, key.1, count);
                }
            }
        });
    }
    map
}

pub(crate) fn get_stats<'a>(
    args: &Parameters,
    genome: &'a GenomeSketch,
    sample: &SequencesSketch,
    winner_map: Option<&FxHashMap<Kmer, (f64, &GenomeSketch, bool)>>,
    log_reassign: bool,
) -> Result<Option<AniResult<'a>>, ProfileError> {
    if genome.k != sample.k || genome.c < sample.c {
        return Err(ProfileError::IncompatibleSketch {
            sample_c: sample.c,
            sample_k: sample.k,
            reference_c: genome.c,
            reference_k: genome.k,
        });
    }
    if genome.k == 0 || genome.c == 0 || sample.c == 0 {
        return Err(ProfileError::InvalidSketch("k and c must be positive"));
    }
    let kmers = &genome.genome_kmers;
    if (kmers.len() as f64) < args.min_number_kmers {
        return Ok(None);
    }
    let mut covs = Vec::new();
    let mut lost = 0usize;
    for kmer in kmers {
        if let Some(&count) = sample.kmer_counts.get(kmer) {
            if count == 0 {
                continue;
            }
            if let Some(map) = winner_map {
                let winner = map
                    .get(kmer)
                    .ok_or(ProfileError::InvalidSketch("missing reassignment winner"))?;
                if winner.1 != genome {
                    lost += 1;
                    continue;
                }
            }
            covs.push(count);
        }
    }
    let contained = covs.len();
    if contained < args.min_contain {
        log::debug!(
            "Discarding {}/{}: {} contained k-mers is below the minimum evidence floor {}",
            genome.file_name,
            genome.first_contig_name,
            contained,
            args.min_contain
        );
        return Ok(None);
    }
    let reassign_log = if winner_map.is_some() && log_reassign {
        Some((
            genome.file_name.as_str(),
            genome.first_contig_name.as_str(),
            lost,
        ))
    } else {
        None
    };
    let Some(fin) = finalize_stats(args, genome.k, kmers.len(), contained, covs, reassign_log)?
    else {
        return Ok(None);
    };
    Ok(Some(AniResult {
        naive_ani: fin.naive_ani,
        final_est_ani: fin.final_est_ani,
        final_est_cov: fin.final_est_cov,
        seq_name: sample
            .sample_name
            .clone()
            .unwrap_or_else(|| sample.file_name.clone()),
        gn_name: genome.file_name.as_str(),
        contig_name: genome.first_contig_name.as_str(),
        mean_cov: fin.mean_cov,
        median_cov: fin.median_cov,
        containment_index: (contained, kmers.len()),
        lambda: fin.lambda,
        ani_ci: fin.ani_ci,
        lambda_ci: fin.lambda_ci,
        genome_sketch: genome,
        rel_abund: None,
        seq_abund: None,
        kmers_lost: winner_map.map(|_| lost),
    }))
}

pub(crate) fn finalize_stats(
    args: &Parameters,
    k: usize,
    n_kmers: usize,
    contain_count: usize,
    mut covs: Vec<u32>,
    reassign_log: Option<(&str, &str, usize)>,
) -> Result<Option<Finalized>, ProfileError> {
    if covs.is_empty() {
        return Ok(None);
    }
    if k == 0
        || n_kmers == 0
        || contain_count > n_kmers
        || contain_count != covs.len()
        || covs.contains(&0)
    {
        return Err(ProfileError::InvalidSketch("invalid coverage evidence"));
    }
    let naive_ani = finite(
        (contain_count as f64 / n_kmers as f64).powf(1.0 / k as f64),
        "naive ANI",
    )?;
    covs.sort();
    let median_cov = covs[covs.len() / 2] as f64;
    let pois = Poisson::new(median_cov)
        .map_err(|_| ProfileError::InvalidSketch("invalid median coverage"))?;
    let mut max_cov = f64::MAX;
    if median_cov < 30.0 {
        for &cov in &covs[covs.len() / 2..] {
            if pois.cdf(cov.into()) < CUTOFF_PVALUE {
                max_cov = cov as f64;
            } else {
                break;
            }
        }
    }
    let mut full_covs = vec![0; n_kmers - contain_count];
    for &cov in &covs {
        if cov as f64 <= max_cov {
            full_covs.push(cov);
        }
    }
    if full_covs.is_empty() {
        return Err(ProfileError::InsufficientEvidence(
            "no coverage after filtering",
        ));
    }
    let sum = checked_sum(&full_covs)? as f64;
    let geq1_mean_cov = sum / covs.len() as f64;
    let use_lambda = if median_cov > MEDIAN_ANI_THRESHOLD {
        AdjustStatus::High
    } else {
        match estimate_lambda(args, &full_covs, k as f64)? {
            None => AdjustStatus::Low,
            Some(lambda) => AdjustStatus::Lambda(finite(lambda, "lambda")?),
        }
    };
    let final_est_cov = if let AdjustStatus::Lambda(lam) = use_lambda {
        lam
    } else if median_cov < MAX_MEDIAN_FOR_MEAN_FINAL_EST || args.mean_coverage {
        geq1_mean_cov
    } else {
        median_cov
    };
    let opt_lambda = if use_lambda == AdjustStatus::Low || use_lambda == AdjustStatus::High {
        None
    } else {
        Some(final_est_cov)
    };
    let estimate = ani_from_lambda(opt_lambda, k as f64, &full_covs);
    let final_est_ani = finite(
        if opt_lambda.is_none() || args.no_adj {
            naive_ani
        } else {
            estimate.unwrap_or(naive_ani)
        },
        "adjusted ANI",
    )?;
    let min_ani = minimum_ani_fraction(args);
    if final_est_ani < min_ani {
        if let Some((gn, ctg, lost)) = reassign_log {
            log::info!("Genome/contig {}/{} has ANI = {} < {} after reassigning {} k-mers ({} contained k-mers after reassign)", gn, ctg, final_est_ani * 100., min_ani * 100., lost, contain_count);
        }
        return Ok(None);
    }
    let (ani_ci, lambda_ci) = if !args.no_ci && opt_lambda.is_some() {
        bootstrap_interval(&full_covs, k as f64, args)?
    } else {
        ((None, None), (None, None))
    };
    Ok(Some(Finalized {
        naive_ani,
        final_est_ani,
        final_est_cov: finite(final_est_cov, "coverage")?,
        median_cov,
        mean_cov: finite(geq1_mean_cov, "mean coverage")?,
        lambda: use_lambda,
        ani_ci,
        lambda_ci,
    }))
}

fn estimate_lambda(
    args: &Parameters,
    covs: &Vec<u32>,
    k: f64,
) -> Result<Option<f64>, ProfileError> {
    // Legacy inference uses signed histogram counters. Bound their frequency,
    // not individual u32 coverage values; bootstrap samples pass this gate too.
    if covs.len() > i32::MAX as usize {
        return Err(ProfileError::ArithmeticOverflow("inference marker count"));
    }
    Ok(if args.ratio {
        ratio_lambda(covs, args.min_count_correct)
    } else if args.mme {
        mme_lambda(covs)
    } else if args.nb {
        binary_search_lambda(covs)
    } else if args.mle {
        mle_zip(covs, k)
    } else {
        ratio_lambda(covs, args.min_count_correct)
    })
}

fn ani_from_lambda(lambda: Option<f64>, k: f64, covs: &[u32]) -> Option<f64> {
    let lambda = lambda?;
    let contained = covs.iter().filter(|&&x| x != 0).count();
    let ani = (contained as f64 / (1.0 - f64::exp(-lambda)) / covs.len() as f64).powf(1.0 / k);
    if ani < 0.0 || ani.is_nan() {
        None
    } else {
        Some(ani)
    }
}

fn bootstrap_interval(
    covs: &[u32],
    k: f64,
    args: &Parameters,
) -> Result<(ConfidenceInterval, ConfidenceInterval), ProfileError> {
    let mut rng = fastrand::Rng::with_seed(DEFAULT_RNG_SEED);
    let mut anis = Vec::new();
    let mut lambdas = Vec::new();
    for _ in 0..100 {
        let mut sampled = Vec::with_capacity(covs.len());
        for _ in 0..covs.len() {
            sampled.push(covs[rng.usize(..covs.len())]);
        }
        checked_sum(&sampled)?;
        let lambda = estimate_lambda(args, &sampled, k)?;
        let ani = ani_from_lambda(lambda, k, &sampled);
        if let (Some(ani), Some(lambda)) = (ani, lambda) {
            if !ani.is_nan() && !lambda.is_nan() {
                anis.push(Some(finite(ani, "bootstrap ANI")?));
                lambdas.push(Some(finite(lambda, "bootstrap lambda")?));
            }
        }
    }
    anis.sort_by(|x, y| x.partial_cmp(y).expect("finite bootstrap ANI"));
    lambdas.sort_by(|x, y| x.partial_cmp(y).expect("finite bootstrap lambda"));
    if anis.len() < 50 {
        return Ok(((None, None), (None, None)));
    }
    let n = anis.len();
    Ok((
        (anis[n * 5 / 100 - 1], anis[n * 95 / 100 - 1]),
        (lambdas[n * 5 / 100 - 1], lambdas[n * 95 / 100 - 1]),
    ))
}

fn get_kmer_identity(sample: &SequencesSketch) -> Result<f64, ProfileError> {
    let mut median = 0u32;
    let mut moving = 0.0;
    let mut n = 1.0;
    for &count in sample.kmer_counts.values() {
        if count > 1 {
            if count > median {
                median += 1;
            } else {
                median -= 1;
            }
            moving += median as f64;
            n += 1.0;
        }
    }
    moving /= n;
    log::debug!(
        "Estimated continuous median k-mer count for {} is {:.3}",
        sample.file_name,
        moving
    );
    let mut ones = 0usize;
    let mut others = 0u32;
    for &count in sample.kmer_counts.values() {
        if count == 1 {
            ones = ones
                .checked_add(1)
                .ok_or(ProfileError::ArithmeticOverflow("singleton marker count"))?;
        } else {
            others = others
                .checked_add(count)
                .ok_or(ProfileError::ArithmeticOverflow(
                    "non-singleton marker count",
                ))?;
        }
    }
    let eps = others as f64 / (others as f64 + ones as f64 + 0.1);
    if moving < MED_KMER_FOR_ID_EST && sample.mean_read_length < 400.0 {
        log::info!("{} short-read sample has high diversity compared to sequencing depth (approx. avg depth < 3). Using 99.5% as read accuracy estimate instead of automatic detection for --estimate-unknown.", sample.file_name);
        Ok(0.995f64.powf(sample.k as f64))
    } else if eps < 1.0 {
        Ok(eps)
    } else {
        Ok(1.0)
    }
}
