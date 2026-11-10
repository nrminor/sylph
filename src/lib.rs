//! Read sketch accumulation is available without default features.
//!
//! Capabilities are nested: `cli` (the default) enables `profiling`, which
//! enables `file-io`. File and profiling entry points retain their legacy
//! process and global-state behavior; use [`read_sketch`] for embedding.

pub mod constants;
pub mod read_sketch;
pub mod seeding;
pub mod sketch;
mod sketch_core;
pub mod types;

#[cfg(feature = "file-io")]
pub mod cmdline;
#[cfg(feature = "profiling")]
pub mod contain;
#[cfg(feature = "profiling")]
pub mod inference;
#[cfg(feature = "file-io")]
pub mod inspect;
#[cfg(feature = "file-io")]
pub mod parallel_sketch;
#[cfg(feature = "file-io")]
mod sketch_files;
#[cfg(feature = "file-io")]
pub mod twostage_db;

#[cfg(target_arch = "x86_64")]
pub mod avx2_seeding;
