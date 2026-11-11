//! Independent driver: scientific work is delegated to the pinned dependency.
use std::{collections::BTreeMap, error::Error, fs, path::Path, process::Command};

use clap::Parser;
use serde::Serialize;
use sylph::cmdline::{Cli, Mode};

#[path = "../profile_options.rs"]
mod profile_options;

const REVISION: &str = "1f7fcb53aa20b592d8b3b97f75d9971a60eb06e9";

fn dna(mut state: u64, length: usize) -> String {
    (0..length)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[(state >> 62) as usize] as char
        })
        .collect()
}

fn reads(reference: &str, step: usize, mutate: bool) -> Vec<String> {
    (0..=reference.len() - 510)
        .step_by(step)
        .map(|offset| {
            let mut read = reference.as_bytes()[offset..offset + 510].to_vec();
            if mutate {
                for index in (75..510).step_by(150) {
                    read[index] = match read[index] {
                        b'A' => b'C',
                        b'C' => b'G',
                        b'G' => b'T',
                        b'T' => b'A',
                        _ => unreachable!(),
                    };
                }
            }
            String::from_utf8(read).unwrap()
        })
        .collect()
}

fn fastq(reads: &[String]) -> String {
    reads
        .iter()
        .enumerate()
        .map(|(index, read)| format!("@read{index}\n{read}\n+\n{}\n", "I".repeat(read.len())))
        .collect()
}

fn run(directory: &Path, arguments: &[&str]) -> Result<(), Box<dyn Error>> {
    let output = Command::new(std::env::current_exe()?)
        .current_dir(directory)
        .args(arguments)
        .output()?;
    if !output.status.success() {
        return Err(format!("{arguments:?}: {}", String::from_utf8_lossy(&output.stderr)).into());
    }
    fs::write(
        directory.join(format!("{}.stderr", arguments.last().unwrap())),
        &output.stderr,
    )?;
    Ok(())
}

fn generate(directory: &Path) -> Result<(), Box<dyn Error>> {
    fs::create_dir_all(directory)?;
    let directory = fs::canonicalize(directory)?;
    let shared = dna(101, 2000);
    let a = format!("{shared}{}", dna(103, 6030));
    let b = format!("{shared}{}", dna(107, 8030));
    fs::write(directory.join("a.fa"), format!(">synthetic_a\n{a}\n"))?;
    fs::write(directory.join("b.fa"), format!(">synthetic_b\n{b}\n"))?;
    let mut sample = vec![dna(109, 510)];
    let first = sample.len();
    sample.extend(reads(&a, 350, false));
    let second = sample.len();
    sample.extend(reads(&b, 80, true));
    let ends = [first, second, sample.len()];
    // Pinned AVX2 processes floor(window_count / 4) windows per lane.
    // These inputs leave no scalar-only trailing windows (not an x86 run).
    for sequence in [&a, &b].into_iter().chain(sample.iter()) {
        assert_eq!((sequence.len() - 31 + 1) % 4, 0);
    }
    run(
        &directory,
        &[
            "sketch",
            "a.fa",
            "b.fa",
            "-c",
            "20",
            "-t",
            "1",
            "--no-sketch-pipeline",
            "-o",
            "reference",
        ],
    )?;
    for (index, end) in ends.into_iter().enumerate() {
        let input = format!("prefix{index}.fq");
        fs::write(directory.join(&input), fastq(&sample[..end]))?;
        for repeat in 0..5 {
            let output = format!("prefix{index}.repeat{repeat}.tsv");
            run(
                &directory,
                &[
                    "profile",
                    "reference.syldb",
                    &input,
                    "-c",
                    "20",
                    "-t",
                    "1",
                    "--no-sketch-pipeline",
                    "--log-reassignments",
                    "-o",
                    &output,
                ],
            )?;
        }
    }
    fs::write(
        directory.join("inputs.yaml"),
        serde_yaml::to_string(&sample)?,
    )?;
    fs::write(directory.join("provenance.yaml"), format!("revision: {REVISION}\ndriver: pinned-public-entrypoints\nread_length: 510\nk: 31\nc: 20\nprefix_ends: {ends:?}\nrepetitions: 5\n"))?;
    run(
        &directory,
        &[
            "query",
            "reference.syldb",
            "prefix2.fq",
            "-c",
            "20",
            "-t",
            "1",
            "--no-sketch-pipeline",
            "-o",
            "initial-query.tsv",
        ],
    )?;
    // Reversing reference insertion order is a tie diagnostic, not a claimed
    // permutation-invariance contract for profiling arbitrary inputs.
    run(
        &directory,
        &[
            "sketch",
            "b.fa",
            "a.fa",
            "-c",
            "20",
            "-t",
            "1",
            "--no-sketch-pipeline",
            "-o",
            "reversed",
        ],
    )?;
    run(
        &directory,
        &[
            "profile",
            "reversed.syldb",
            "prefix2.fq",
            "-c",
            "20",
            "-t",
            "1",
            "--no-sketch-pipeline",
            "--log-reassignments",
            "-o",
            "reversed.tsv",
        ],
    )?;
    Ok(())
}

#[derive(Serialize)]
struct Prefix {
    end: usize,
    input_name: String,
    arguments: Vec<String>,
    rows: Vec<BTreeMap<String, String>>,
}

#[derive(Serialize)]
struct Fixture {
    revision: &'static str,
    oracle: &'static str,
    rustc: String,
    target: &'static str,
    extraction_backend: &'static str,
    precision_and_limits: &'static str,
    checksums_sha256: String,
    database_arguments: Vec<String>,
    genome_fastas: BTreeMap<String, String>,
    c: usize,
    k: usize,
    no_dedup: bool,
    reads: Vec<String>,
    repetitions: usize,
    reversed_reference_diagnostic_equal: bool,
    reassignment_diagnostics: BTreeMap<String, Vec<String>>,
    initial_query_rows: Vec<BTreeMap<String, String>>,
    prefixes: Vec<Prefix>,
}

fn rows(path: &Path) -> Result<Vec<BTreeMap<String, String>>, Box<dyn Error>> {
    let text = fs::read_to_string(path)?;
    let mut lines = text.lines();
    let headers: Vec<_> = lines
        .next()
        .ok_or("missing TSV header")?
        .split('\t')
        .collect();
    lines
        .map(|line| {
            let values: Vec<_> = line.split('\t').collect();
            assert_eq!(values.len(), headers.len(), "malformed oracle TSV");
            Ok(headers
                .iter()
                .zip(values)
                .map(|(key, value)| (key.to_string(), value.to_string()))
                .collect())
        })
        .collect()
}

fn command_stdout(command: &mut Command) -> Result<String, Box<dyn Error>> {
    let output = command.output()?;
    if !output.status.success() {
        return Err(String::from_utf8_lossy(&output.stderr).into_owned().into());
    }
    Ok(String::from_utf8(output.stdout)?)
}

fn frozen_prefixes(source: &Path) -> Result<Vec<Prefix>, Box<dyn Error>> {
    let mut prefixes = Vec::new();
    for (index, end) in [1, 23, 143].into_iter().enumerate() {
        let expected = rows(&source.join(format!("prefix{index}.repeat0.tsv")))?;
        for repeat in 1..5 {
            assert_eq!(
                expected,
                rows(&source.join(format!("prefix{index}.repeat{repeat}.tsv")))?,
                "unstable pinned oracle"
            );
        }
        prefixes.push(Prefix {
            end,
            input_name: format!("prefix{index}.fq"),
            arguments: [
                "profile",
                "reference.syldb",
                &format!("prefix{index}.fq"),
                "-c",
                "20",
                "-t",
                "1",
                "--no-sketch-pipeline",
                "--log-reassignments",
                "-o",
                &format!("prefix{index}.repeat0.tsv"),
            ]
            .into_iter()
            .map(str::to_owned)
            .collect(),
            rows: expected,
        });
    }
    assert_eq!(prefixes[0].rows.len(), 0);
    assert_eq!(prefixes[1].rows.len(), 1);
    assert_eq!(prefixes[2].rows.len(), 2);
    assert_eq!(
        prefixes[2].rows,
        rows(&source.join("reversed.tsv"))?,
        "reference-order diagnostic changed"
    );
    Ok(prefixes)
}

fn reassignment_diagnostics(
    source: &Path,
) -> Result<BTreeMap<String, Vec<String>>, Box<dyn Error>> {
    let mut diagnostics = BTreeMap::new();
    for name in ["prefix1.repeat0.tsv", "prefix2.repeat0.tsv", "reversed.tsv"] {
        let log = fs::read_to_string(source.join(format!("{name}.stderr")))?;
        let events: Vec<String> = log
            .lines()
            .filter_map(|line| line.split_once("[sylph::contain] "))
            .map(|(_, event)| event)
            .filter(|event| {
                event.starts_with("Index\t")
                    || event.contains("kmers reassigned")
                    || event.contains("reassigning k-mers")
                    || event.contains("passing profiling threshold")
            })
            .map(str::to_owned)
            .collect();
        diagnostics.insert(name.into(), events);
    }
    let reversed = &diagnostics["reversed.tsv"];
    assert!(reversed
        .iter()
        .any(|line| line == "Index\t0\tb.fa\tsynthetic_b"));
    assert!(reversed
        .iter()
        .any(|line| line == "Index\t1\ta.fa\tsynthetic_a"));
    // Pinned winner_table replaces an existing winner only on strict `>`.
    // A visited second still wins all 45 shared markers: this cohort has a
    // strict winner, despite both initial ANI TSV values being capped at 100.
    assert!(reversed
        .iter()
        .any(|line| line == "1->0\t45\tkmers reassigned"));
    Ok(diagnostics)
}

fn freeze(source: &Path, destination: &Path) -> Result<(), Box<dyn Error>> {
    let reads: Vec<String> =
        serde_yaml::from_str(&fs::read_to_string(source.join("inputs.yaml"))?)?;
    assert_eq!(
        reads.len(),
        143,
        "cohort changed: re-characterize explicitly"
    );
    let prefixes = frozen_prefixes(source)?;
    let mut checksum = Command::new("shasum");
    checksum.current_dir(source).args([
        "-a",
        "256",
        "reference.syldb",
        "a.fa",
        "b.fa",
        "inputs.yaml",
        "prefix0.fq",
        "prefix1.fq",
        "prefix2.fq",
    ]);
    for name in ["Cargo.toml", "Cargo.lock", "src/bin/profile-reference.rs"] {
        checksum.arg(Path::new(env!("CARGO_MANIFEST_DIR")).join(name));
    }
    let mut genome_fastas = BTreeMap::new();
    for name in ["a.fa", "b.fa"] {
        genome_fastas.insert(name.into(), fs::read_to_string(source.join(name))?);
    }
    let fixture = Fixture {
        revision: REVISION,
        oracle: "Independent driver calling unchanged pinned public Cli, sketch::sketch and contain::contain; not an original upstream binary. Synthetic DNA has no intended biological function.",
        rustc: command_stdout(Command::new("rustc").arg("-vV"))?,
        target: env!("REFERENCE_TARGET"),
        extraction_backend: extraction_backend(),
        precision_and_limits: "TSV semantic values only: abundance 4 decimals; coverage/lambda/mean 3 decimals; ANI and CI endpoints 2 decimals. Adjusted ANI is capped at 100 in TSV, not internally. No unrounded full-profile oracle. Five fresh-process invocations agree; not a general order-independence guarantee, particularly for -u. All sequences have (length-k+1)%4=0, avoiding pinned AVX2 trailing-window omission; no x86 execution claimed. All reads are 510 bases, bypassing <=400-base signature dedup and short-read -u fallback. Reversed-reference diagnostic agrees but is not an arbitrary permutation contract.",
        checksums_sha256: command_stdout(&mut checksum)?,
        database_arguments: ["sketch", "a.fa", "b.fa", "-c", "20", "-t", "1", "--no-sketch-pipeline", "-o", "reference"].into_iter().map(str::to_owned).collect(),
        genome_fastas,
        c: 20, k: 31, no_dedup: false, reads,
        repetitions: 5,
        reversed_reference_diagnostic_equal: true,
        reassignment_diagnostics: reassignment_diagnostics(source)?,
        initial_query_rows: rows(&source.join("initial-query.tsv"))?,
        prefixes,
    };
    fs::create_dir_all(destination)?;
    fs::copy(
        source.join("reference.syldb"),
        destination.join("reference.syldb"),
    )?;
    fs::write(
        destination.join("default.yaml"),
        serde_yaml::to_string(&fixture)?,
    )?;
    Ok(())
}

fn extraction_backend() -> &'static str {
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") {
        return "avx2";
    }
    "scalar"
}

fn main() -> Result<(), Box<dyn Error>> {
    let arguments: Vec<_> = std::env::args_os().collect();
    if arguments
        .get(1)
        .is_some_and(|value| value == "generate-options")
    {
        assert_eq!(
            arguments.len(),
            3,
            "usage: profile-reference generate-options SOURCE"
        );
        return profile_options::generate(Path::new(&arguments[2]));
    }
    if arguments
        .get(1)
        .is_some_and(|value| value == "freeze-options")
    {
        assert_eq!(
            arguments.len(),
            4,
            "usage: profile-reference freeze-options SOURCE OUTPUT.yaml"
        );
        return profile_options::freeze(Path::new(&arguments[2]), Path::new(&arguments[3]));
    }
    if arguments.get(1).is_some_and(|value| value == "generate") {
        assert_eq!(
            arguments.len(),
            3,
            "usage: profile-reference generate DIRECTORY"
        );
        return generate(Path::new(&arguments[2]));
    }
    if arguments.get(1).is_some_and(|value| value == "freeze") {
        assert_eq!(
            arguments.len(),
            4,
            "usage: profile-reference freeze SOURCE DESTINATION"
        );
        return freeze(Path::new(&arguments[2]), Path::new(&arguments[3]));
    }
    match Cli::parse().mode {
        Mode::Sketch(args) => sylph::sketch::sketch(args),
        Mode::Profile(args) => sylph::contain::contain(args, true),
        Mode::Query(args) => sylph::contain::contain(args, false),
        _ => return Err("driver supports only sketch, profile, query, generate".into()),
    }
    Ok(())
}
