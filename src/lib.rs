//! Read sketch accumulation is available without default features.
//!
//! `cli` (the default) enables legacy `profiling`, which enables `file-io` and
//! `profile-core`. Enable `profile-core` alone for reusable standard-database
//! profiling through `profile::Profiler`, without FASTX or native compression.
//! Legacy command entry points retain their process and global-state behavior.

pub mod constants;
pub mod read_sketch;
pub mod seeding;
pub mod sketch;
mod sketch_core;
pub mod types;

#[cfg(feature = "profile-core")]
pub mod profile;
#[cfg(feature = "profile-core")]
mod profile_core;

#[cfg(feature = "file-io")]
pub mod cmdline;
#[cfg(feature = "profiling")]
pub mod contain;
#[cfg(feature = "profiling")]
pub mod inference;
#[cfg(all(feature = "profile-core", not(feature = "profiling")))]
mod inference;
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
