//! Ordered, single-end sketch accumulation without file or process operations.

use fxhash::FxHashSet;
use thiserror::Error;

use crate::{
    constants::MAX_DEDUP_COUNT,
    sketch_core::{dup_removal_lsh_full_exact, extract_markers, pair_kmer_single, Marker},
    types::SequencesSketch,
};

/// The single-end duplicate policy, fixed for an accumulator's lifetime.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Deduplication {
    /// Count every sampled marker occurrence.
    Disabled,
    /// Use sylph's exact signature sets, including its count-four escape.
    Exact,
}

/// Scientific parameters checked when constructing a read sketcher.
#[derive(Clone, Copy, Debug)]
pub struct SketchParams {
    /// Sampling factor. Must be greater than zero.
    pub c: usize,
    /// K-mer length. Must be 21 or 31, supported by both extraction backends.
    pub k: usize,
    /// Single-end duplicate policy.
    pub deduplication: Deduplication,
}

/// A rejected configuration or read. Read indices are local to an `extend` call.
#[derive(Debug, Error, Eq, PartialEq)]
pub enum SketchError {
    /// The extraction backend cannot accept these parameters.
    #[error("c must be positive and k must be 21 or 31; got c={c}, k={k}")]
    InvalidParameters {
        /// Requested sampling factor.
        c: usize,
        /// Requested k-mer length.
        k: usize,
    },
    /// An empty sequence or a byte outside the IUPAC nucleotide alphabet.
    #[error("read {read_index} must contain nonempty IUPAC nucleotide letters")]
    InvalidRead {
        /// Number of complete reads accepted earlier in this call.
        read_index: usize,
    },
    /// A marker count would exceed the sketch format's `u32` range.
    #[error("read {read_index} would overflow a marker count")]
    MarkerCountOverflow {
        /// Number of complete reads accepted earlier in this call.
        read_index: usize,
    },
    /// Another read would exceed the sequential mean's exact integer counter.
    #[error("read {read_index} would exceed the exact read-count limit of 2^53")]
    ReadCountLimit {
        /// Number of complete reads accepted earlier in this call.
        read_index: usize,
    },
}

/// Persistent counts, duplicate history, and mean length for one ordered sample.
///
/// Batch boundaries and sketch observations do not reset this state. A sketch
/// alone does not contain the history needed to resume accumulation.
pub struct ReadSketcher {
    sketch: SequencesSketch,
    deduplication: Deduplication,
    signatures: FxHashSet<(u64, [Marker; 2])>,
    markers: Vec<u64>,
    reads: f64,
}

impl ReadSketcher {
    /// Start an empty sketch with a logical input name; no file is opened.
    ///
    /// # Errors
    /// Rejects a zero sampling factor or a k-mer length other than 21 or 31.
    pub fn new(
        params: SketchParams,
        file_name: String,
        sample_name: Option<String>,
    ) -> Result<Self, SketchError> {
        if params.c == 0 || !matches!(params.k, 21 | 31) {
            return Err(SketchError::InvalidParameters {
                c: params.c,
                k: params.k,
            });
        }
        Ok(Self {
            sketch: SequencesSketch::new(file_name, params.c, params.k, false, sample_name, 0.0),
            deduplication: params.deduplication,
            signatures: FxHashSet::default(),
            markers: Vec::new(),
            reads: 0.0,
        })
    }

    /// Add complete, borrowed sequences in order, without retaining the reads.
    ///
    /// Letter case and IUPAC ambiguity use the existing sylph byte mapping.
    /// Marker extraction uses the same CPU dispatch as whole-file sketching.
    /// The sequential running mean includes every accepted read, even one with
    /// no markers or one whose markers are all suppressed as duplicates.
    ///
    /// # Errors
    /// Rejects invalid reads and unrepresentable counts. Earlier complete reads
    /// remain accepted; the rejected read changes neither counts, duplicate
    /// history, nor mean. The error's index identifies that accepted prefix.
    pub fn extend<'a>(
        &mut self,
        reads: impl IntoIterator<Item = &'a [u8]>,
    ) -> Result<(), SketchError> {
        for (read_index, sequence) in reads.into_iter().enumerate() {
            if sequence.is_empty()
                || !sequence
                    .iter()
                    .all(|base| b"ACGTURYSWKMBDHVNacgturyswkmbdhvn".contains(base))
            {
                return Err(SketchError::InvalidRead { read_index });
            }
            if self.reads >= 9_007_199_254_740_992.0 {
                return Err(SketchError::ReadCountLimit { read_index });
            }
            self.markers.clear();
            extract_markers(sequence, &mut self.markers, self.sketch.c, self.sketch.k);
            let signature = if sequence.len() > 400 {
                None
            } else {
                pair_kmer_single(sequence)
            };
            if self.would_overflow(signature) {
                return Err(SketchError::MarkerCountOverflow { read_index });
            }
            let mut removed = 0;
            for marker in &self.markers {
                dup_removal_lsh_full_exact(
                    &mut self.sketch.kmer_counts,
                    &mut self.signatures,
                    marker,
                    signature,
                    &mut removed,
                    self.deduplication == Deduplication::Disabled,
                    Some(MAX_DEDUP_COUNT),
                );
            }
            self.reads += 1.0;
            #[allow(
                clippy::cast_precision_loss,
                reason = "Preserve the sequential reference's read-length conversion"
            )]
            let length = sequence.len() as f64;
            self.sketch.mean_read_length =
                self.sketch.mean_read_length + (length - self.sketch.mean_read_length) / self.reads;
        }
        Ok(())
    }

    /// Borrow the current sketch without copying counts or resetting history.
    #[must_use]
    pub const fn sketch(&self) -> &SequencesSketch {
        &self.sketch
    }

    fn would_overflow(&self, signature: Option<([Marker; 2], [Marker; 2])>) -> bool {
        let bound = u64::try_from(self.markers.len()).unwrap_or(u64::MAX);
        self.markers.iter().any(|marker| {
            let count = self.sketch.kmer_counts.get(marker).copied().unwrap_or(0);
            let remaining = u64::from(u32::MAX - count);
            if bound <= remaining {
                return false;
            }
            // Only near capacity: count occurrences without allocating or
            // changing the extraction order. All state is still untouched.
            let occurrences = self.markers.iter().filter(|other| *other == marker).count();
            if u64::try_from(occurrences).unwrap_or(u64::MAX) <= remaining {
                return false;
            }
            if self.deduplication == Deduplication::Exact && count < MAX_DEDUP_COUNT {
                if let Some((first, second)) = signature {
                    if count > 0
                        && (first == second
                            || self.signatures.contains(&(*marker, first))
                            || self.signatures.contains(&(*marker, second)))
                    {
                        return false;
                    }
                    // One new observation can raise the count. Further copies
                    // are suppressed unless that first one reaches the escape.
                    return count + 1 >= MAX_DEDUP_COUNT;
                }
            }
            true
        })
    }
}
