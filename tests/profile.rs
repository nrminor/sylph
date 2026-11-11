use std::{collections::BTreeMap, path::Path};

use serde::Deserialize;
use sylph::{
    profile::{ProfileError, ProfileMode, ProfileOptions, ProfileRow, Profiler, SequenceEstimate},
    read_sketch::{Deduplication, ReadSketcher, SketchParams},
    types::{AdjustStatus, GenomeSketch, SequencesSketch},
};

#[derive(Deserialize)]
struct Fixture {
    revision: String,
    c: usize,
    k: usize,
    no_dedup: bool,
    reads: Vec<String>,
    prefixes: Vec<Prefix>,
}

#[derive(Deserialize)]
struct Prefix {
    end: usize,
    rows: Vec<BTreeMap<String, String>>,
}

#[derive(Deserialize)]
struct OptionsFixture {
    revision: String,
    fixed_cohort: String,
    prefix_end: usize,
    cases: Vec<OptionCase>,
}

#[derive(Deserialize)]
struct OptionCase {
    name: String,
    rows: Vec<BTreeMap<String, String>>,
}

fn assert_rounded(actual: f64, expected: &str, decimals: i32, context: &str) {
    let expected: f64 = expected.parse().expect("numeric pinned oracle value");
    let tolerance = 0.5 * 10.0_f64.powi(-decimals) + 1e-12;
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "{context}: actual {actual}, pinned CLI {expected}, tolerance {tolerance}"
    );
}

fn assert_interval(actual: Option<(f64, f64)>, expected: &str, scale: f64, context: &str) {
    if expected == "NA-NA" {
        assert_eq!(actual, None, "{context}");
        return;
    }
    let (low, high) = expected.split_once('-').expect("pinned interval endpoints");
    let (actual_low, actual_high) = actual.expect("pinned oracle has an interval");
    assert_rounded(actual_low * scale, low, 2, &format!("{context} lower"));
    assert_rounded(actual_high * scale, high, 2, &format!("{context} upper"));
}

fn assert_row(
    actual: &ProfileRow,
    expected: &BTreeMap<String, String>,
    mode: ProfileMode,
    context: &str,
) {
    let context = format!("{context}, {}", actual.genome_name);
    let coverage_field = if expected.contains_key("True_cov") {
        "True_cov"
    } else {
        "Eff_cov"
    };
    // The pinned CLI caps only this displayed ANI, not the scientific result.
    for (field, value, decimals) in [
        ("Adjusted_ANI", actual.adjusted_ani.min(1.0) * 100.0, 2),
        ("Naive_ANI", actual.naive_ani * 100.0, 2),
        (coverage_field, actual.effective_coverage, 3),
        ("Mean_cov_geq1", actual.mean_coverage, 3),
        ("Taxonomic_abundance", actual.taxonomic_abundance_percent, 4),
    ] {
        assert_rounded(
            value,
            &expected[field],
            decimals,
            &format!("{context}, {field}"),
        );
    }
    assert!(
        actual.adjusted_ani.is_finite(),
        "{context}: raw ANI is finite"
    );
    assert_eq!(
        actual.median_coverage,
        expected["Median_cov"].parse::<f64>().unwrap(),
        "{context}"
    );
    let (contained, total) = expected["Containment_ind"].split_once('/').unwrap();
    assert_eq!(
        actual.containment,
        (contained.parse().unwrap(), total.parse().unwrap()),
        "{context}"
    );
    assert_eq!(
        actual.kmers_reassigned,
        expected["kmers_reassigned"].parse::<usize>().unwrap(),
        "{context}"
    );
    let abundance = match (mode, &actual.sequence_estimate) {
        (ProfileMode::EstimateReadCounts, SequenceEstimate::EstimatedReadCount(value)) => *value,
        (
            ProfileMode::RelativeAbundance | ProfileMode::EstimateUnknown,
            SequenceEstimate::Percent(value),
        ) => *value,
        _ => panic!("{context}: incorrect sequence-estimate units for {mode:?}"),
    };
    assert_rounded(
        abundance,
        &expected["Sequence_abundance"],
        4,
        &format!("{context}, sequence abundance"),
    );
    match expected["Eff_lambda"].as_str() {
        "HIGH" => assert!(
            matches!(&actual.adjustment, AdjustStatus::High),
            "{context}"
        ),
        "LOW" => assert!(matches!(&actual.adjustment, AdjustStatus::Low), "{context}"),
        value => match &actual.adjustment {
            AdjustStatus::Lambda(lambda) => {
                assert_rounded(lambda.to_owned(), value, 3, &format!("{context}, lambda"))
            }
            other => panic!("{context}: expected numeric lambda, got {other:?}"),
        },
    }
    assert_interval(
        actual.ani_interval,
        &expected["ANI_5-95_percentile"],
        100.0,
        &format!("{context}, ANI interval"),
    );
    assert_interval(
        actual.lambda_interval,
        &expected["Lambda_5-95_percentile"],
        1.0,
        &format!("{context}, lambda interval"),
    );
}

fn assert_rows(
    actual: &[ProfileRow],
    expected: &[BTreeMap<String, String>],
    mode: ProfileMode,
    context: &str,
) {
    let mut actual_ids: Vec<_> = actual
        .iter()
        .map(|row| (row.genome_name.as_str(), row.contig_name.as_str()))
        .collect();
    let mut expected_ids: Vec<_> = expected
        .iter()
        .map(|row| (row["Genome_file"].as_str(), row["Contig_name"].as_str()))
        .collect();
    actual_ids.sort_unstable();
    expected_ids.sort_unstable();
    assert_eq!(actual_ids, expected_ids, "{context} final identities");
    for expected in expected {
        let actual = actual
            .iter()
            .find(|row| {
                row.genome_name == expected["Genome_file"]
                    && row.contig_name == expected["Contig_name"]
            })
            .unwrap();
        assert_row(actual, expected, mode, context);
    }
}

#[test]
fn default_profile_matches_pinned_raw_input_prefixes() {
    let fixture: Fixture =
        serde_yaml::from_str(include_str!("fixtures/profile/default.yaml")).unwrap();
    assert_eq!(fixture.revision, "1f7fcb53aa20b592d8b3b97f75d9971a60eb06e9");
    assert_eq!(
        fixture
            .prefixes
            .iter()
            .map(|prefix| prefix.end)
            .collect::<Vec<_>>(),
        [1, 23, 143]
    );
    let database =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/profile/reference.syldb");
    let profiler = Profiler::open_syldb(database, ProfileOptions::default(), 1).unwrap();
    let mut sketcher = ReadSketcher::new(
        SketchParams {
            c: fixture.c,
            k: fixture.k,
            deduplication: if fixture.no_dedup {
                Deduplication::Disabled
            } else {
                Deduplication::Exact
            },
        },
        "instrument-input".to_owned(),
        Some("synthetic-sample".to_owned()),
    )
    .unwrap();

    let mut start = 0;
    for prefix in fixture.prefixes {
        sketcher
            .extend(
                fixture.reads[start..prefix.end]
                    .iter()
                    .map(String::as_bytes),
            )
            .unwrap();
        let profile = profiler.profile(sketcher.sketch()).unwrap();
        // Logical metadata is ours; raw-prefix oracle filenames are not identities.
        assert_eq!(profile.sample_name, "synthetic-sample");
        assert_rows(
            &profile.rows,
            &prefix.rows,
            ProfileMode::RelativeAbundance,
            &format!("prefix {}", prefix.end),
        );
        start = prefix.end;
    }
}

fn options_for_case(name: &str) -> ProfileOptions {
    let mut options = ProfileOptions::default();
    match name {
        "minimum_ani" => options.minimum_ani = 99.9,
        "min_contain" => options.min_contain = 163,
        "min_number_kmers" => options.min_number_kmers = 190.0,
        "confidence_intervals_disabled" => options.confidence_intervals = false,
        "unknown_auto_identity" => options.mode = ProfileMode::EstimateUnknown,
        "estimated_read_counts" => options.mode = ProfileMode::EstimateReadCounts,
        "unknown_explicit_identity" => {
            options.mode = ProfileMode::EstimateUnknown;
            options.read_identity = Some(99.0);
        }
        _ => panic!("unexpected pinned option case: {name}"),
    }
    options
}

#[test]
fn supported_options_match_pinned_raw_input_profile() {
    let fixture: Fixture =
        serde_yaml::from_str(include_str!("fixtures/profile/default.yaml")).unwrap();
    let controls: OptionsFixture =
        serde_yaml::from_str(include_str!("fixtures/profile/options.yaml")).unwrap();
    assert_eq!(controls.revision, fixture.revision);
    assert_eq!(controls.fixed_cohort, "default.yaml");
    assert_eq!(controls.prefix_end, fixture.reads.len());
    assert_eq!(controls.cases.len(), 7);
    let mut sketcher = ReadSketcher::new(
        SketchParams {
            c: fixture.c,
            k: fixture.k,
            deduplication: Deduplication::Exact,
        },
        "instrument-input".to_owned(),
        Some("synthetic-sample".to_owned()),
    )
    .unwrap();
    sketcher
        .extend(fixture.reads.iter().map(String::as_bytes))
        .unwrap();
    let database =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/profile/reference.syldb");
    for case in controls.cases {
        let options = options_for_case(&case.name);
        let profiler = Profiler::open_syldb(&database, options, 1).unwrap();
        let profile = profiler.profile(sketcher.sketch()).unwrap();
        assert_eq!(profile.sample_name, "synthetic-sample", "{}", case.name);
        assert_eq!(
            profile.options, options,
            "{} fixed options metadata",
            case.name
        );
        assert_rows(&profile.rows, &case.rows, options.mode, &case.name);
    }
}

fn marker_reference() -> GenomeSketch {
    GenomeSketch {
        genome_kmers: (1..=10).collect(),
        pseudotax_tracked_nonused_kmers: Some(Vec::new()),
        file_name: "marker-reference".to_owned(),
        first_contig_name: "marker-contig".to_owned(),
        c: 1,
        k: 31,
        gn_size: 1000,
        min_spacing: 30,
    }
}

fn marker_sample(count: u32) -> SequencesSketch {
    let mut sample = SequencesSketch::new(
        "marker-input".to_owned(),
        1,
        31,
        false,
        Some("marker-sample".to_owned()),
        510.0,
    );
    sample
        .kmer_counts
        .extend((1..=10).map(|marker| (marker, count)));
    sample
}

#[test]
fn calculation_errors_preserve_reusable_profiler_and_owned_results() {
    let mut sample = marker_sample(2);
    let (retained, repeated) = {
        let profiler =
            Profiler::from_genomes(vec![marker_reference()], ProfileOptions::default(), 1).unwrap();
        let retained = profiler.profile(&sample).unwrap();
        assert_eq!(retained.rows.len(), 1);

        sample.kmer_counts.insert(1, 0);
        assert!(matches!(
            profiler.profile(&sample),
            Err(ProfileError::InvalidSketch(_))
        ));
        sample.kmer_counts.insert(1, 2);
        sample.k = 21;
        assert!(matches!(
            profiler.profile(&sample),
            Err(ProfileError::IncompatibleSketch {
                sample_k: 21,
                reference_k: 31,
                ..
            })
        ));
        sample.k = 31;
        sample.c = 2;
        assert!(matches!(
            profiler.profile(&sample),
            Err(ProfileError::IncompatibleSketch {
                sample_c: 2,
                reference_c: 1,
                ..
            })
        ));
        sample.c = 1;

        // With all counts maximal, the high median disables outlier trimming;
        // the retained u32 sum overflows without a huge marker allocation.
        sample
            .kmer_counts
            .values_mut()
            .for_each(|count| *count = u32::MAX);
        let overflow = profiler.profile(&sample);
        assert!(
            matches!(&overflow, Err(ProfileError::ArithmeticOverflow(_))),
            "expected retained-count overflow, got {overflow:?}"
        );
        sample.kmer_counts.values_mut().for_each(|count| *count = 2);
        let repeated = profiler.profile(&sample).unwrap();
        assert_eq!(
            repeated, retained,
            "errors must not poison the reusable reference"
        );
        (retained, repeated)
    };

    // Both returned profiles outlive the profiler and the source sample.
    sample.kmer_counts.clear();
    sample.sample_name = Some("changed-sample".to_owned());
    drop(sample);
    assert_eq!(retained, repeated);
    assert_eq!(retained.sample_name, "marker-sample");
    assert_eq!(retained.rows[0].genome_name, "marker-reference");
    assert_eq!(retained.rows[0].contig_name, "marker-contig");
    assert_eq!(retained.rows[0].containment, (10, 10));
    assert_eq!(retained.rows[0].effective_coverage, 2.0);
    assert_eq!(
        retained.rows[0].sequence_estimate,
        SequenceEstimate::Percent(100.0)
    );
}

#[test]
fn correction_requires_evidence_only_when_final_hits_remain() {
    for mode in [
        ProfileMode::EstimateUnknown,
        ProfileMode::EstimateReadCounts,
    ] {
        let options = ProfileOptions {
            mode,
            ..ProfileOptions::default()
        };
        let profiler = Profiler::from_genomes(vec![marker_reference()], options, 1).unwrap();
        let mut sample = marker_sample(1);
        sample.kmer_counts.clear();
        sample.mean_read_length = 0.0;
        let empty = profiler.profile(&sample).unwrap();
        assert_eq!(empty.rows.len(), 0, "{mode:?}: valid empty input");

        // A valid nonempty sample can have neither reference hits nor evidence
        // for automatic read identity. No correction is needed for zero rows.
        sample.mean_read_length = 510.0;
        sample.kmer_counts.insert(999, 1);
        let unmatched = profiler.profile(&sample).unwrap();
        assert_eq!(unmatched.rows.len(), 0, "{mode:?}: nonmatching markers");

        // Ten contained singleton markers give a positive hit but cannot
        // estimate automatic identity; long reads bypass short-read fallback.
        let positive = marker_sample(1);
        let undefined = profiler.profile(&positive);
        assert!(
            matches!(&undefined, Err(ProfileError::InsufficientEvidence(_))),
            "{mode:?}: positive undefined correction must fail, got {undefined:?}"
        );

        // The same matching evidence becomes a valid no-hit result when the
        // requested containment threshold removes every candidate.
        let filtered_profiler = Profiler::from_genomes(
            vec![marker_reference()],
            ProfileOptions {
                min_contain: 11,
                ..options
            },
            1,
        )
        .unwrap();
        let filtered = filtered_profiler.profile(&positive).unwrap();
        assert_eq!(filtered.rows.len(), 0, "{mode:?}: all candidates filtered");
    }
}

#[test]
fn constructor_rejects_invalid_options_and_references_without_affecting_valid_owners() {
    let first = marker_reference();
    let mut second = marker_reference();
    second.file_name = "other-reference".to_owned();
    second.first_contig_name = "other-contig".to_owned();
    second.genome_kmers = (11..=20).collect();
    let valid = vec![first.clone(), second.clone()];
    let options = ProfileOptions::default();
    let before = Profiler::from_genomes(valid.clone(), options, 2).unwrap();
    let sample = marker_sample(2);
    let expected = before.profile(&sample).unwrap();
    assert_eq!(expected.rows.len(), 1);

    let zero_workers = Profiler::from_genomes(valid.clone(), options, 0)
        .err()
        .expect("zero workers must fail");
    assert!(
        matches!(&zero_workers, ProfileError::InvalidOptions(_)),
        "{zero_workers:?}"
    );
    let invalid_options = ProfileOptions {
        minimum_ani: f64::NAN,
        ..options
    };
    let invalid_option = Profiler::from_genomes(valid.clone(), invalid_options, 2)
        .err()
        .expect("NaN ANI must fail");
    assert!(
        matches!(&invalid_option, ProfileError::InvalidOptions(_)),
        "{invalid_option:?}"
    );

    let mut missing = second.clone();
    missing.pseudotax_tracked_nonused_kmers = None;
    let mut incompatible = second.clone();
    incompatible.k = 21;
    let mut overlap = second;
    overlap.pseudotax_tracked_nonused_kmers = Some(vec![11]);
    let mut zero_size = first.clone();
    zero_size.gn_size = 0;
    for (label, genomes, expected_index) in [
        ("empty database", Vec::new(), 0),
        (
            "later missing profiling data",
            vec![first.clone(), missing],
            1,
        ),
        ("later inconsistent k", vec![first.clone(), incompatible], 1),
        ("duplicate identity", vec![first.clone(), first.clone()], 1),
        ("tracked/main marker overlap", vec![first, overlap], 1),
        ("zero genome size", vec![zero_size], 0),
    ] {
        let error = Profiler::from_genomes(genomes, options, 2)
            .err()
            .expect(label);
        assert!(
            matches!(&error, ProfileError::InvalidReference { index, .. } if *index == expected_index),
            "{label}: expected invalid reference {expected_index}, got {error:?}"
        );
    }

    let after = Profiler::from_genomes(valid, options, 2).unwrap();
    assert_eq!(after.profile(&sample).unwrap(), expected);
    assert_eq!(before.profile(&sample).unwrap(), expected);
}
