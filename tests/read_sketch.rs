use std::collections::BTreeMap;

use serde::Deserialize;
use sylph::{
    read_sketch::{Deduplication, ReadSketcher, SketchError, SketchParams},
    types::SequencesSketch,
};

#[derive(Deserialize)]
struct Fixtures {
    revision: String,
    cases: Vec<Case>,
}

#[derive(Deserialize)]
struct Case {
    name: String,
    c: usize,
    k: usize,
    no_dedup: bool,
    reads: Vec<String>,
    prefixes: Vec<Expected>,
}

#[derive(Deserialize)]
struct Expected {
    kmer_counts: BTreeMap<u64, u32>,
    c: usize,
    k: usize,
    file_name: String,
    sample_name: Option<String>,
    paired: bool,
    mean_read_length_bits: u64,
}

fn fixtures() -> Fixtures {
    let fixtures: Fixtures =
        serde_yaml::from_str(include_str!("fixtures/read-sketch.yaml")).unwrap();
    assert_eq!(
        fixtures.revision,
        "1f7fcb53aa20b592d8b3b97f75d9971a60eb06e9"
    );
    for case in &fixtures.cases {
        assert_eq!(case.reads.len(), case.prefixes.len(), "{}", case.name);
    }
    fixtures
}

fn sketcher(case: &Case) -> ReadSketcher {
    ReadSketcher::new(
        SketchParams {
            c: case.c,
            k: case.k,
            deduplication: if case.no_dedup {
                Deduplication::Disabled
            } else {
                Deduplication::Exact
            },
        },
        "reads.fastq".to_owned(),
        Some("fixture".to_owned()),
    )
    .unwrap()
}

fn assert_sketch(actual: &SequencesSketch, expected: &Expected) {
    let counts: BTreeMap<_, _> = actual
        .kmer_counts
        .iter()
        .map(|(&key, &n)| (key, n))
        .collect();
    assert_eq!(counts, expected.kmer_counts);
    assert_eq!(actual.c, expected.c);
    assert_eq!(actual.k, expected.k);
    assert_eq!(actual.file_name, expected.file_name);
    assert_eq!(actual.sample_name, expected.sample_name);
    assert_eq!(actual.paired, expected.paired);
    assert_eq!(
        actual.mean_read_length.to_bits(),
        expected.mean_read_length_bits
    );
}

#[test]
fn duplicate_policy_survives_a_batch_boundary() {
    let fixtures = fixtures();
    let cases: Vec<_> = fixtures
        .cases
        .iter()
        .filter(|case| case.c == 1 && case.reads.len() == 2 && case.reads[0] == case.reads[1])
        .collect();
    assert!(cases.iter().any(|case| case.no_dedup));
    assert!(cases.iter().any(|case| !case.no_dedup));
    for case in cases {
        assert!(!case.prefixes[0].kmer_counts.is_empty());
        assert!(case.prefixes[0]
            .kmer_counts
            .values()
            .all(|&count| count == 1));
        if case.no_dedup {
            assert_eq!(
                case.prefixes[0].kmer_counts.len(),
                case.prefixes[1].kmer_counts.len()
            );
            assert!(case.prefixes[1]
                .kmer_counts
                .values()
                .all(|&count| count == 2));
        } else {
            assert_eq!(case.prefixes[0].kmer_counts, case.prefixes[1].kmer_counts);
        }
        let mut sketcher = sketcher(case);
        for (read, expected) in case.reads.iter().zip(&case.prefixes) {
            sketcher.extend([read.as_bytes()]).unwrap();
            assert_sketch(sketcher.sketch(), expected);
        }
    }
}

#[test]
fn every_partition_and_observed_prefix_matches_upstream() {
    for case in fixtures().cases {
        assert!(!case.reads.is_empty() && case.reads.len() <= 8);
        let partitions = 1usize << (case.reads.len() - 1);
        for cuts in 0..partitions {
            let mut sketcher = sketcher(&case);
            assert!(sketcher.sketch().kmer_counts.is_empty());
            assert_eq!(
                sketcher.sketch().mean_read_length.to_bits(),
                0.0f64.to_bits()
            );
            let mut start = 0;
            for end in 1..=case.reads.len() {
                if end == case.reads.len() || cuts & (1 << (end - 1)) != 0 {
                    sketcher
                        .extend(case.reads[start..end].iter().map(String::as_bytes))
                        .unwrap();
                    assert_sketch(sketcher.sketch(), &case.prefixes[end - 1]);
                    // Repeated observations and empty batches must not disturb
                    // the state used when the next nonempty batch arrives.
                    assert_sketch(sketcher.sketch(), &case.prefixes[end - 1]);
                    sketcher.extend(std::iter::empty::<&[u8]>()).unwrap();
                    assert_sketch(sketcher.sketch(), &case.prefixes[end - 1]);
                    start = end;
                }
            }
        }
    }
}

#[test]
fn invalid_input_preserves_the_accepted_prefix() {
    let fixtures = fixtures();
    let case = fixtures
        .cases
        .iter()
        .find(|case| !case.no_dedup && case.c == 1 && case.k == 21 && case.reads.len() == 2)
        .unwrap();
    let mut sketcher = sketcher(case);
    let error = sketcher
        .extend([case.reads[0].as_bytes(), b"ACZ".as_slice()])
        .unwrap_err();
    assert_eq!(error, SketchError::InvalidRead { read_index: 1 });
    assert_sketch(sketcher.sketch(), &case.prefixes[0]);
    assert_eq!(
        sketcher.extend([b"".as_slice()]).unwrap_err(),
        SketchError::InvalidRead { read_index: 0 }
    );
    sketcher.extend([case.reads[1].as_bytes()]).unwrap();
    assert_sketch(sketcher.sketch(), &case.prefixes[1]);

    for (c, k) in [(0, 21), (1, 0), (1, 20), (1, 32)] {
        let result = ReadSketcher::new(
            SketchParams {
                c,
                k,
                deduplication: Deduplication::Exact,
            },
            "reads.fastq".into(),
            None,
        );
        assert!(matches!(result, Err(SketchError::InvalidParameters { .. })));
    }
}
