use std::{collections::BTreeMap, error::Error, fs, path::PathBuf, process::Command};

use serde::Serialize;
use sylph::{seeding::fmh_seeds, sketch::sketch_sequences_needle};

const REVISION: &str = "1f7fcb53aa20b592d8b3b97f75d9971a60eb06e9";
const FILE_NAME: &str = "reads.fastq";

#[derive(Serialize)]
struct Reference {
    revision: String,
    rustc: String,
    target: String,
    extraction_backend: String,
    cases: Vec<Case>,
}

#[derive(Serialize)]
struct Case {
    name: String,
    c: usize,
    k: usize,
    no_dedup: bool,
    reads: Vec<String>,
    prefixes: Vec<Expected>,
}

#[derive(Serialize)]
struct Expected {
    kmer_counts: BTreeMap<u64, u32>,
    c: usize,
    k: usize,
    file_name: String,
    sample_name: Option<String>,
    paired: bool,
    mean_read_length_bits: u64,
}

// Deterministic arithmetic input construction, not a sketch/hash implementation.
fn dna(mut state: u64, length: usize) -> String {
    (0..length)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[(state >> 62) as usize] as char
        })
        .collect()
}

fn boundary_read(seed: u64, length: usize) -> String {
    let mut read = dna(seed, length - 24);
    read.push_str(&"A".repeat(24));
    read
}

fn inputs() -> Vec<(&'static str, usize, usize, Vec<String>)> {
    let duplicate = dna(17, 200);
    let duplicate_k31 = dna(17, 202);
    let short = dna(31, 20);
    let below = boundary_read(53, 65);
    let at = boundary_read(79, 66);
    let long = boundary_read(101, 401);
    let mut signatures: Vec<String> = b"ACGT"
        .iter()
        .map(|base| {
            let mut read = duplicate.as_bytes().to_vec();
            // Each signature samples one of these two positions. The rest of
            // the read is shared, exposing the count-four dedup escape.
            read[0] = *base;
            read[1] = *base;
            String::from_utf8(read).unwrap()
        })
        .collect();
    signatures.push(signatures[0].clone());
    vec![
        (
            "duplicate_across_boundary",
            1,
            21,
            vec![duplicate.clone(), duplicate],
        ),
        (
            "short_and_signature_boundary_65_66",
            8,
            21,
            vec![short, below.clone(), below, at.clone(), at],
        ),
        (
            "dedup_boundary_400_401",
            16,
            21,
            vec![long[..400].into(), long[..400].into(), long.clone(), long],
        ),
        ("count_four_escape", 8, 21, signatures),
        (
            "duplicate_across_boundary_k31",
            1,
            31,
            vec![duplicate_k31.clone(), duplicate_k31],
        ),
    ]
}

fn oracle(
    reads: &[String],
    c: usize,
    k: usize,
    no_dedup: bool,
) -> Result<Option<Expected>, Box<dyn Error>> {
    let mut fastq = String::new();
    for (index, read) in reads.iter().enumerate() {
        fastq.push_str(&format!(
            "@read{index}\n{read}\n+\n{}\n",
            "I".repeat(read.len())
        ));
    }
    fs::write(FILE_NAME, fastq)?;
    Ok(
        sketch_sequences_needle(FILE_NAME, c, k, Some("fixture".into()), no_dedup).map(|sketch| {
            Expected {
                kmer_counts: sketch.kmer_counts.into_iter().collect(),
                c: sketch.c,
                k: sketch.k,
                file_name: sketch.file_name,
                sample_name: sketch.sample_name,
                paired: sketch.paired,
                mean_read_length_bits: sketch.mean_read_length.to_bits(),
            }
        }),
    )
}

fn generate() -> Result<Vec<Case>, Box<dyn Error>> {
    let mut cases = Vec::new();
    for (name, c, k, reads) in inputs() {
        assert!(reads.len() <= 8);
        verify_backend_neutral_inputs(&reads, c, k);
        for no_dedup in [false, true] {
            let mut prefixes = Vec::new();
            if let Some(empty) = oracle(&[], c, k, no_dedup)? {
                prefixes.push(empty);
            }
            for end in 1..=reads.len() {
                prefixes
                    .push(oracle(&reads[..end], c, k, no_dedup)?.ok_or("invalid fixture FASTQ")?);
            }
            // Empty FASTQ is currently rejected upstream. Fail fast if this
            // changes rather than silently changing prefix indexing below.
            assert_eq!(
                prefixes.len(),
                reads.len(),
                "upstream now accepts empty FASTQ"
            );
            cases.push(Case {
                name: format!("{name}_{}", if no_dedup { "no_dedup" } else { "dedup" }),
                c,
                k,
                no_dedup,
                reads: reads.clone(),
                prefixes,
            });
        }
    }
    verify_distinctions(&cases);
    Ok(cases)
}

// Pinned AVX2 extraction omits the final 1-3 windows when the window
// count is not divisible by four. This bounded input subset avoids that
// discrepancy; it is NOT a general cross-backend equivalence guarantee.
// Check omitted windows using upstream scalar extraction, never a custom
// hash/model. Expected sketches still come only from the public file oracle.
fn verify_backend_neutral_inputs(reads: &[String], c: usize, k: usize) {
    for read in reads {
        if read.len() < k {
            continue;
        }
        let windows = read.len() - k + 1;
        let retained = windows / 4 * 4;
        for offset in retained..windows {
            let mut sampled = Vec::new();
            fmh_seeds(&read.as_bytes()[offset..offset + k], &mut sampled, c, k);
            assert!(sampled.is_empty(), "AVX2-omitted window would be sampled");
        }
    }
}

fn verify_distinctions(cases: &[Case]) {
    let map = |case: usize, prefix: usize| &cases[case].prefixes[prefix].kmer_counts;
    assert!(!map(0, 0).is_empty(), "duplicate oracle must have markers");
    assert_eq!(map(0, 0), map(0, 1), "duplicate suppressed across boundary");
    for (marker, count) in map(0, 0) {
        assert_eq!(*count, 1);
        assert_eq!(map(1, 1)[marker], 2, "no-dedup counts duplicate");
    }
    assert!(map(2, 0).is_empty(), "shorter than k produces no markers");
    assert!(!map(2, 1).is_empty());
    assert_ne!(map(2, 1), map(2, 2), "65-base repeats are counted");
    assert_ne!(map(2, 2), map(2, 3), "66-base read contributes markers");
    assert_eq!(map(2, 3), map(2, 4), "66-base repeats are suppressed");
    assert_ne!(map(2, 4), map(3, 4));
    assert!(!map(4, 0).is_empty());
    assert_eq!(map(4, 0), map(4, 1), "400-base repeats are suppressed");
    assert_ne!(map(4, 1), map(4, 2));
    assert_ne!(map(4, 2), map(4, 3), "401-base repeats are counted");
    assert_ne!(map(4, 3), map(5, 3));
    let shared: Vec<_> = map(6, 0)
        .keys()
        .filter(|marker| {
            (0..5).all(|prefix| map(6, prefix).get(marker) == Some(&(prefix as u32 + 1)))
        })
        .collect();
    assert!(
        !shared.is_empty(),
        "shared marker must progress 1,2,3,4,5 despite final repeated signature"
    );
    assert_eq!(map(6, 4), map(7, 4), "count-four escape matches no-dedup");
    assert!(!map(8, 0).is_empty(), "k=31 oracle must have markers");
    assert_eq!(map(8, 0), map(8, 1), "k=31 duplicate suppressed");
    for (marker, count) in map(8, 0) {
        assert_eq!(*count, 1);
        assert_eq!(map(9, 1)[marker], 2, "k=31 no-dedup counts duplicate");
    }
    assert_ne!(map(0, 0), map(8, 0), "comparable DNA differs at k=21/31");
    eprintln!(
        "oracle distinctions verified; {} count-four shared markers",
        shared.len()
    );
}

struct TemporaryCwd {
    original: PathBuf,
    directory: PathBuf,
}

impl TemporaryCwd {
    fn enter() -> Result<Self, Box<dyn Error>> {
        let original = std::env::current_dir()?;
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)?
            .as_nanos();
        let directory =
            std::env::temp_dir().join(format!("sylph-reference-{}-{nonce}", std::process::id()));
        fs::create_dir(&directory)?;
        let temporary = Self {
            original,
            directory,
        };
        std::env::set_current_dir(&temporary.directory)?;
        Ok(temporary)
    }
}

impl Drop for TemporaryCwd {
    fn drop(&mut self) {
        std::env::set_current_dir(&self.original).expect("restore original working directory");
        fs::remove_dir_all(&self.directory).expect("remove oracle temporary directory");
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let backend = extraction_backend();
    let mut arguments = std::env::args_os().skip(1);
    let output = arguments.next().map(PathBuf::from).unwrap_or_else(|| {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../fixtures/read-sketch.yaml")
    });
    assert!(
        arguments.next().is_none(),
        "usage: sylph-read-sketch-reference [output.yaml]"
    );
    let output = std::path::absolute(output)?;
    let rustc = Command::new("rustc").arg("-vV").output()?;
    assert!(rustc.status.success(), "rustc metadata command failed");
    let metadata = String::from_utf8(rustc.stdout)?;
    let _temporary = TemporaryCwd::enter()?;
    let reference = Reference {
        revision: REVISION.into(),
        rustc: metadata
            .lines()
            .next()
            .ok_or("missing rustc version")?
            .into(),
        target: env!("REFERENCE_TARGET").into(),
        extraction_backend: backend.into(),
        cases: generate()?,
    };
    let yaml = serde_yaml::to_string(&reference)?;
    assert!(yaml.len() < 300_000, "fixture exceeds size budget");
    fs::create_dir_all(output.parent().ok_or("missing fixture directory")?)?;
    fs::write(&output, &yaml)?;
    eprintln!(
        "wrote {} cases, {} bytes to {}",
        reference.cases.len(),
        yaml.len(),
        output.display()
    );
    Ok(())
}

// Match the pinned upstream public extractor's runtime dispatch, without
// replacing its implementation or changing inputs for either backend.
fn extraction_backend() -> &'static str {
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") {
        return "avx2";
    }
    "scalar"
}
