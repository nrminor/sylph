//! Reusable, synchronous profiling without CLI or process-wide initialization.
//!
//! References and scientific options are fixed when constructing a [`Profiler`].
//! Each call computes a complete profile, including marker reassignment and
//! abundance normalization, and returns an owned result. It does not incrementally
//! update a previous profile or retain the sample sketch.

use std::{fs::File, io, io::BufReader, path::Path};

use fxhash::FxHashSet;
use rayon::{ThreadPool, ThreadPoolBuilder};
use thiserror::Error;

use crate::{
    profile_core::{self, Parameters},
    types::{AdjustStatus, AniResult, GenomeSketch, SequencesSketch},
};

/// How sequence abundance is estimated. Taxonomic abundance remains a percentage.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum ProfileMode {
    /// Normalize sequence abundance among detected reference genomes.
    #[default]
    RelativeAbundance,
    /// Scale sequence abundance by the estimated fraction represented in the database.
    EstimateUnknown,
    /// Estimate read counts instead of sequence percentages, including unknown correction.
    /// This is sylph's experimental read-count estimator, not a read assignment count.
    EstimateReadCounts,
}

/// Scientific controls, with the same defaults as `sylph profile`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ProfileOptions {
    /// Minimum adjusted ANI in percent (0 through 100); default 95.
    /// Upstream warns that profiling below 95 percent can be inaccurate.
    pub minimum_ani: f64,
    /// Minimum contained markers per hit; default 7.
    pub min_contain: usize,
    /// Minimum sampled markers in a reference sketch; default 10.
    /// Fractional thresholds retain the existing CLI's comparison semantics.
    pub min_number_kmers: f64,
    /// Request confidence intervals; default true. Some estimates have no interval.
    pub confidence_intervals: bool,
    /// Sequence abundance mode; default relative abundance.
    pub mode: ProfileMode,
    /// Read identity in percent, greater than zero and at most 100.
    /// `None` uses upstream's automatic estimator. An override requires one of
    /// the unknown-correction modes. The automatic short-read fallback depends
    /// on marker-map iteration order; equivalent maps need not give equal estimates.
    pub read_identity: Option<f64>,
}

impl Default for ProfileOptions {
    fn default() -> Self {
        Self {
            minimum_ani: 95.0,
            min_contain: 7,
            min_number_kmers: 10.0,
            confidence_intervals: true,
            mode: ProfileMode::RelativeAbundance,
            read_identity: None,
        }
    }
}

/// A sequence-abundance estimate with explicit units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SequenceEstimate {
    /// Sequence abundance as a percentage.
    Percent(f64),
    /// Rounded experimental read-count estimate, not an observed integer count.
    EstimatedReadCount(f64),
}

/// An owned, final profile. Options distinguish disabled from unavailable intervals.
#[derive(Debug, PartialEq)]
pub struct Profile {
    /// Logical sample name, falling back to the sketch's input name.
    pub sample_name: String,
    /// The fixed scientific options used for this calculation.
    pub options: ProfileOptions,
    /// Final rows in descending taxonomic abundance order; tied order is unspecified.
    pub rows: Vec<ProfileRow>,
}

/// Final scientific estimates for one reference, without retaining its sketch.
#[derive(Debug, PartialEq)]
pub struct ProfileRow {
    /// Reference input name stored in the database.
    pub genome_name: String,
    /// First contig name stored in the database.
    pub contig_name: String,
    /// Naive ANI as a fraction, not a percentage.
    pub naive_ani: f64,
    /// Unrounded adjusted ANI as a fraction. Upstream estimates can exceed one;
    /// only the CLI's displayed adjusted ANI is capped at 100 percent.
    pub adjusted_ani: f64,
    /// Final effective coverage, including unknown correction when requested.
    pub effective_coverage: f64,
    /// The upstream `Mean_cov_geq1` estimate.
    pub mean_coverage: f64,
    /// Median matched-marker coverage.
    pub median_coverage: f64,
    /// Contained markers and total sampled reference markers, after reassignment.
    pub containment: (usize, usize),
    /// Coverage-adjustment status; `Low` and `High` are not numeric lambda estimates.
    pub adjustment: AdjustStatus,
    /// ANI interval in fractional units, when available and requested.
    pub ani_interval: Option<(f64, f64)>,
    /// Lambda interval before unknown correction, when available and requested.
    pub lambda_interval: Option<(f64, f64)>,
    /// Taxonomic abundance in percent among final detected references.
    pub taxonomic_abundance_percent: f64,
    /// Sequence percentage or experimental estimated read count.
    pub sequence_estimate: SequenceEstimate,
    /// Number of contained markers assigned to competing references.
    pub kmers_reassigned: usize,
}

/// A configuration, reference, sample, or calculation failure. No partial profile is returned.
#[derive(Debug, Error)]
pub enum ProfileError {
    /// Scientific options or worker count are invalid.
    #[error("invalid profile configuration: {0}")]
    InvalidOptions(&'static str),
    /// A reference fails the profiling admission contract.
    #[error("invalid reference at index {index}: {reason}")]
    InvalidReference {
        /// Zero-based position in the loaded reference vector.
        index: usize,
        /// Violated invariant.
        reason: &'static str,
    },
    /// Sample metadata or marker counts are invalid.
    #[error("invalid sample sketch: {0}")]
    InvalidSketch(&'static str),
    /// The sample must have matching k and be at least as dense as each reference.
    #[error("incompatible sample (c={sample_c}, k={sample_k}) and reference (c={reference_c}, k={reference_k})")]
    IncompatibleSketch {
        /// Sample sampling factor.
        sample_c: usize,
        /// Sample k-mer length.
        sample_k: usize,
        /// Reference sampling factor.
        reference_c: usize,
        /// Reference k-mer length.
        reference_k: usize,
    },
    /// A calculation exceeds an upstream integer accumulator's range.
    #[error("{0} exceeds the supported arithmetic range")]
    ArithmeticOverflow(&'static str),
    /// The requested correction cannot be estimated from the available evidence.
    #[error("insufficient evidence: {0}")]
    InsufficientEvidence(&'static str),
    /// A calculation produced an undefined or non-finite value.
    #[error("undefined profile estimate: {0}")]
    NonFiniteEstimate(&'static str),
    /// Reference file I/O failed.
    #[error("could not read reference database: {0}")]
    Io(#[from] io::Error),
    /// Reference bytes are not a standard sylph database.
    #[error("could not decode standard sylph database: {0}")]
    Decode(#[from] bincode::Error),
    /// The dedicated worker pool could not be created.
    #[error("could not create profiling workers: {0}")]
    ThreadPool(#[from] rayon::ThreadPoolBuildError),
}

/// One validated reference collection, fixed options, and a dedicated local worker pool.
///
/// Construction does not initialize logging or Rayon's global pool. Calculations
/// are synchronous and are not cancelled by dropping an async caller's future.
/// Hosts should run them outside async executor threads and control concurrency.
pub struct Profiler {
    genomes: Vec<GenomeSketch>,
    options: ProfileOptions,
    parameters: Parameters,
    pool: ThreadPool,
    k: usize,
    minimum_c: usize,
}

impl Profiler {
    /// Load a standard `.syldb` once and retain its reference sketches in memory.
    /// Two-stage databases, raw FASTA, and compressed inputs are not accepted.
    /// This loader is for operator-provided local databases, not untrusted uploads.
    ///
    /// # Errors
    /// Returns I/O, decoding, validation, or worker-pool errors.
    pub fn open_syldb(
        path: impl AsRef<Path>,
        options: ProfileOptions,
        workers: usize,
    ) -> Result<Self, ProfileError> {
        let reader = BufReader::new(File::open(path)?);
        let genomes = bincode::deserialize_from(reader)?;
        Self::from_genomes(genomes, options, workers)
    }

    /// Adopt already-loaded reference sketches without copying them.
    ///
    /// # Errors
    /// Rejects an empty database, incompatible reference parameters, missing
    /// profiling markers, duplicate reference identities or markers, invalid
    /// options, zero workers, or failure to create the local worker pool.
    pub fn from_genomes(
        genomes: Vec<GenomeSketch>,
        options: ProfileOptions,
        workers: usize,
    ) -> Result<Self, ProfileError> {
        options.validate()?;
        if workers == 0 {
            return Err(ProfileError::InvalidOptions("workers must be positive"));
        }
        let (k, minimum_c) = validate_references(&genomes)?;
        let pool = ThreadPoolBuilder::new().num_threads(workers).build()?;
        Ok(Self {
            genomes,
            options,
            parameters: options.parameters(),
            pool,
            k,
            minimum_c,
        })
    }

    /// Calculate a complete profile of this sketch without retaining or mutating it.
    /// Valid empty/no-hit samples return an empty result. Detections can change
    /// or disappear when later evidence is profiled; this is not incremental inference.
    ///
    /// # Errors
    /// Rejects invalid or incompatible sketches, unrepresentable arithmetic,
    /// and undefined estimates in the requested mode. No partial rows are returned.
    pub fn profile(&self, sketch: &SequencesSketch) -> Result<Profile, ProfileError> {
        if sketch.c == 0 || !matches!(sketch.k, 21 | 31) {
            return Err(ProfileError::InvalidSketch(
                "c must be positive and k must be 21 or 31",
            ));
        }
        if sketch.k != self.k || sketch.c > self.minimum_c {
            return Err(ProfileError::IncompatibleSketch {
                sample_c: sketch.c,
                sample_k: sketch.k,
                reference_c: self.minimum_c,
                reference_k: self.k,
            });
        }
        if !sketch.mean_read_length.is_finite()
            || sketch.mean_read_length < 0.0
            || (!sketch.kmer_counts.is_empty() && sketch.mean_read_length == 0.0)
        {
            return Err(ProfileError::InvalidSketch(
                "mean read length must be finite and positive for a nonempty sketch",
            ));
        }
        if sketch.kmer_counts.values().any(|count| *count == 0) {
            return Err(ProfileError::InvalidSketch(
                "stored marker counts must be positive",
            ));
        }
        let references: Vec<_> = self.genomes.iter().collect();
        let results = self
            .pool
            .install(|| profile_core::calculate(&self.parameters, &references, sketch))?;
        Ok(Profile {
            sample_name: sketch
                .sample_name
                .as_ref()
                .unwrap_or(&sketch.file_name)
                .clone(),
            options: self.options,
            rows: results
                .into_iter()
                .map(|result| owned_row(result, self.options.mode))
                .collect(),
        })
    }
}

impl ProfileOptions {
    fn validate(&self) -> Result<(), ProfileError> {
        if !self.minimum_ani.is_finite() || !(0.0..=100.0).contains(&self.minimum_ani) {
            return Err(ProfileError::InvalidOptions(
                "minimum ANI must be between 0 and 100 percent",
            ));
        }
        if !self.min_number_kmers.is_finite() || self.min_number_kmers < 0.0 {
            return Err(ProfileError::InvalidOptions(
                "minimum reference markers must be finite and nonnegative",
            ));
        }
        if let Some(identity) = self.read_identity {
            if !identity.is_finite() || identity <= 0.0 || identity > 100.0 {
                return Err(ProfileError::InvalidOptions(
                    "read identity must be greater than 0 and at most 100 percent",
                ));
            }
            if self.mode == ProfileMode::RelativeAbundance {
                return Err(ProfileError::InvalidOptions(
                    "read identity requires unknown correction",
                ));
            }
        }
        Ok(())
    }

    fn parameters(self) -> Parameters {
        Parameters {
            minimum_ani: Some(self.minimum_ani),
            pseudotax: true,
            min_number_kmers: self.min_number_kmers,
            min_contain: self.min_contain,
            min_count_correct: 3.0,
            ratio: false,
            mme: false,
            nb: false,
            mle: false,
            no_adj: false,
            no_ci: !self.confidence_intervals,
            mean_coverage: false,
            seq_id: self.read_identity,
            estimate_unknown: self.mode != ProfileMode::RelativeAbundance,
            estimate_read_counts: self.mode == ProfileMode::EstimateReadCounts,
            redundant_ani: 99.0,
            log_reassignments: false,
        }
    }
}

fn validate_references(genomes: &[GenomeSketch]) -> Result<(usize, usize), ProfileError> {
    let first = genomes.first().ok_or(ProfileError::InvalidReference {
        index: 0,
        reason: "database is empty",
    })?;
    let mut minimum_c = first.c;
    let mut identities = FxHashSet::default();
    for (index, genome) in genomes.iter().enumerate() {
        let invalid = |reason| ProfileError::InvalidReference { index, reason };
        if genome.c == 0 || !matches!(genome.k, 21 | 31) || genome.k != first.k {
            return Err(invalid(
                "c must be positive and all references must share k=21 or k=31",
            ));
        }
        if genome.gn_size == 0 || genome.genome_kmers.len() > genome.gn_size {
            return Err(invalid(
                "genome size must be positive and cover its sampled markers",
            ));
        }
        let tracked = genome
            .pseudotax_tracked_nonused_kmers
            .as_ref()
            .ok_or_else(|| invalid("database lacks required profiling markers"))?;
        if !identities.insert((&genome.file_name, &genome.first_contig_name)) {
            return Err(invalid("duplicate reference input and contig name"));
        }
        let mut markers = FxHashSet::default();
        for marker in genome.genome_kmers.iter().chain(tracked) {
            if !markers.insert(*marker) {
                return Err(invalid("duplicate marker within a reference"));
            }
        }
        minimum_c = minimum_c.min(genome.c);
    }
    Ok((first.k, minimum_c))
}

fn owned_row(result: AniResult<'_>, mode: ProfileMode) -> ProfileRow {
    let sequence = result
        .seq_abund
        .expect("complete profiling assigns sequence abundance");
    ProfileRow {
        genome_name: result.gn_name.to_owned(),
        contig_name: result.contig_name.to_owned(),
        naive_ani: result.naive_ani,
        adjusted_ani: result.final_est_ani,
        effective_coverage: result.final_est_cov,
        mean_coverage: result.mean_cov,
        median_coverage: result.median_cov,
        containment: result.containment_index,
        adjustment: result.lambda,
        ani_interval: result.ani_ci.0.zip(result.ani_ci.1),
        lambda_interval: result.lambda_ci.0.zip(result.lambda_ci.1),
        taxonomic_abundance_percent: result
            .rel_abund
            .expect("complete profiling assigns taxonomic abundance"),
        sequence_estimate: match mode {
            ProfileMode::EstimateReadCounts => SequenceEstimate::EstimatedReadCount(sequence),
            ProfileMode::RelativeAbundance | ProfileMode::EstimateUnknown => {
                SequenceEstimate::Percent(sequence)
            }
        },
        kmers_reassigned: result
            .kmers_lost
            .expect("complete profiling includes reassignment"),
    }
}
