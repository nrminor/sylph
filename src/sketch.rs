//! Marker extraction and, with `file-io`, legacy file sketching entry points.

pub use crate::sketch_core::{extract_markers, extract_markers_positions};

#[cfg(feature = "file-io")]
pub use crate::sketch_files::{
    check_vram_and_block, is_fasta, is_fastq, sketch, sketch_genome, sketch_genome_individual,
    sketch_pair_sequences, sketch_sequences_needle,
};
