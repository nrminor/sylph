use std::{collections::BTreeMap, error::Error, fs, path::Path, process::Command};

use serde::Serialize;

use super::{command_stdout, extraction_backend, rows, run, REVISION};

fn cases() -> [(&'static str, &'static [&'static str]); 7] {
    [
        ("minimum_ani", &["--minimum-ani", "99.9"]),
        ("min_contain", &["--min-contain", "163"]),
        ("min_number_kmers", &["--min-number-kmers", "190"]),
        ("confidence_intervals_disabled", &["--no-ci"]),
        ("unknown_auto_identity", &["-u"]),
        ("estimated_read_counts", &["--estimate-read-counts"]),
        ("unknown_explicit_identity", &["-u", "-I", "99"]),
    ]
}

fn arguments<'a>(controls: &'a [&'a str], output: &'a str) -> Vec<&'a str> {
    let mut arguments = vec![
        "profile",
        "reference.syldb",
        "prefix2.fq",
        "-c",
        "20",
        "-t",
        "1",
        "--no-sketch-pipeline",
        "--log-reassignments",
    ];
    arguments.extend_from_slice(controls);
    arguments.extend_from_slice(&["-o", output]);
    arguments
}

pub fn generate(source: &Path) -> Result<(), Box<dyn Error>> {
    // Existing cohort only: never rebuild its database or regenerate defaults.
    let frozen = Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures/profile/reference.syldb");
    assert_eq!(
        fs::read(source.join("reference.syldb"))?,
        fs::read(frozen)?,
        "wrong frozen cohort database"
    );
    for (name, controls) in cases() {
        for repeat in 0..5 {
            let output = format!("options-{name}.repeat{repeat}.tsv");
            run(source, &arguments(controls, &output))?;
        }
    }
    Ok(())
}

#[derive(Serialize)]
struct Case {
    name: &'static str,
    arguments: Vec<String>,
    rows: Vec<BTreeMap<String, String>>,
}

#[derive(Serialize)]
struct OptionsFixture {
    revision: &'static str,
    oracle: &'static str,
    fixed_cohort: &'static str,
    prefix_end: usize,
    repetitions: usize,
    rustc: String,
    target: &'static str,
    extraction_backend: &'static str,
    precision_and_limits: &'static str,
    fixed_cohort_sha256: String,
    raw_input_sha256: String,
    generator_sha256: String,
    cases: Vec<Case>,
}

fn frozen_cases(source: &Path) -> Result<Vec<Case>, Box<dyn Error>> {
    let default = rows(&source.join("prefix2.repeat0.tsv"))?;
    let mut frozen = Vec::new();
    for (name, controls) in cases() {
        let output = format!("options-{name}.repeat0.tsv");
        let expected = rows(&source.join(&output))?;
        for repeat in 1..5 {
            let repeated = source.join(format!("options-{name}.repeat{repeat}.tsv"));
            assert_eq!(expected, rows(&repeated)?, "unstable pinned option {name}");
        }
        assert_ne!(expected, default, "option {name} has no observable effect");
        frozen.push(Case {
            name,
            arguments: arguments(controls, &output)
                .into_iter()
                .map(str::to_owned)
                .collect(),
            rows: expected,
        });
    }
    Ok(frozen)
}

pub fn freeze(source: &Path, output: &Path) -> Result<(), Box<dyn Error>> {
    let cases = frozen_cases(source)?;
    let manifest = Path::new(env!("CARGO_MANIFEST_DIR"));
    let fixture = OptionsFixture {
        revision: REVISION,
        oracle: "Independent driver calling unchanged pinned public CLI parser and contain::contain, not an original upstream binary. No adapted scientific bodies.",
        fixed_cohort: "default.yaml",
        prefix_end: 143,
        repetitions: 5,
        rustc: command_stdout(Command::new("rustc").arg("-vV"))?,
        target: env!("REFERENCE_TARGET"),
        extraction_backend: extraction_backend(),
        precision_and_limits: "Same rounded-semantic precision and backend limitations as default.yaml. Five fresh raw-input processes per configuration agree. Every read is 510 bases, excluding the order-sensitive short-read fallback; this does not establish general -u map-order invariance. No unrounded full-profile oracle. Sequence_abundance is percent except estimated_read_counts, which is a rounded experimental read-count estimate printed to four decimals. Effective coverage includes unknown correction; Eff_lambda and its interval do not.",
        fixed_cohort_sha256: command_stdout(Command::new("shasum").current_dir(manifest.join("../fixtures/profile")).args(["-a", "256", "default.yaml", "reference.syldb"]))?,
        raw_input_sha256: command_stdout(Command::new("shasum").current_dir(source).args(["-a", "256", "prefix2.fq"]))?,
        generator_sha256: command_stdout(Command::new("shasum").current_dir(manifest).args(["-a", "256", "Cargo.toml", "Cargo.lock", "src/bin/profile-reference.rs", "src/profile_options.rs"]))?,
        cases,
    };
    fs::write(output, serde_yaml::to_string(&fixture)?)?;
    Ok(())
}
