# sylph - fast and precise species-level metagenomic + taxonomic profiling for shotgun sequencing 

<p align="center"><img src="assets/jay_logo.png" width = 450 /></p>
 
> [!IMPORTANT]
> Documentation for sylph has moved to https://sylph-docs.github.io/. All GitHub documentation (e.g., Wikis) are out of date. 

**Sylph** is a tool for ultrafast taxonomic profiling of metagenomic shotgun reads. Sylph can detect species and their abundances against massive databases such as GTDB (~200k genomes) in seconds. 

### Why sylph?

1. **Ultrafast, multithreaded, multi-sample**: sylph can be > 50x faster than other methods. Sylph only takes ~5 GB of RAM for profiling against the entire GTDB-R232 database (200k genomes).

2. **Precise species-level profiling**: sylph has less false positives than Kraken and is about as precise and sensitive as marker gene methods (MetaPhlAn, mOTUs). 

3. **Accurate (containment) ANI information**: sylph can give accurate **ANI estimates** between reference genomes and your metagenome sample down to 0.1x coverage.

4. **Customizable databases and pre-built databases**: We offer pre-built databases of [prokaryotes, viruses, eukaryotes](pre‐built-databases.md). Custom databases (e.g. using your own MAGs) are easy to build.  

5. **Short or long reads**: Sylph works for nanopore, PacBio, and short reads. Sylph is the most accurate method [on Oxford Nanopore's independent benchmarks](https://nanoporetech.com/resource-centre/genomic-and-epigenomic-insights-into-microbial-biology-with-nanopore-metagenomic-and-isolate-sequencing).

Documentation, installation, and usage information is available at https://sylph-docs.github.io/.

## Incremental single-end sketches

This fork adds `sylph::read_sketch::ReadSketcher`. It accepts ordered batches of complete sequences without opening files, initializing a logger or thread pool, or retaining the reads. Parameters are fixed at construction. The current interface supports `k = 21` or `31` and a positive sampling factor `c`.

For an embedded accumulator, select this fork with `default-features = false` in the Cargo dependency declaration. The core retains the existing scientific types and their Serde implementations, but does not compile file readers, compression backends, CLI parsing, or profiling dependencies. Calls to `ReadSketcher` are unchanged.

Legacy features include their dependencies; `profile-core` can also be selected without file I/O:

```text
default -> cli -> profiling --+--> file-io
                             +--> profile-core

core is always available
```

| Selection | Additional capabilities |
| --- | --- |
| No default features | `read_sketch`, `types`, constants, seeding, and `sketch` marker extraction |
| `file-io` | File sketching and its parallel pipeline, database conversion, inspection, and plain configuration structs in `cmdline` |
| `profile-core` | Reusable profiling and standard `.syldb` loading, without FASTX or native compression |
| `profiling` | Existing query/profile routines and inference, including `file-io` and `profile-core` |
| `cli` (default) | Complete `sylph` executable and Clap parsing for the configuration structs |

The default build preserves the existing CLI and public module paths. The binary requires `cli`; there is no reduced-command CLI. `file-io` retains the existing gzip/bzip2/xz readers and native `zlib-ng` selection. A core-only consumer can choose its own compression backend without sylph enabling one. Cargo features are additive: another dependency enabling sylph's defaults will also enable them for that build.

Legacy command entry points retain their existing global logger/thread-pool initialization and process-exit behavior. Do not call `contain::contain` from a service; use `profile::Profiler` instead. Configuration structs without Clap are plain data; Rust's `SketchArgs::default()` still supplies field defaults, not the CLI's scientific defaults.

```rust
use sylph::read_sketch::{Deduplication, ReadSketcher, SketchParams};

let mut sketcher = ReadSketcher::new(
    SketchParams { c: 200, k: 31, deduplication: Deduplication::Exact },
    "instrument-input".to_owned(),
    Some("sample-1".to_owned()),
)?;
sketcher.extend(batch.iter().map(Vec::as_slice))?;
let current = sketcher.sketch();
```

`sketch()` borrows the current result without copying the count map or resetting duplicate history. It is not a recovery checkpoint: `SequencesSketch` does not contain the state needed to resume deduplication. Batch boundaries must preserve read order. Paired reads, profile computation, transport, and restart recovery are not part of this interface.

`Deduplication::Exact` preserves sylph's single-end rules, including signature eligibility and the count-four escape; it does not suppress duplicates forever. `Disabled` counts every sampled occurrence. The running mean includes all accepted reads, including reads with no markers. On an ingestion error, earlier complete reads remain accepted and the rejected read leaves the sketch unchanged. The error reports its zero-based index within that call. Marker-count overflow and a read count beyond the exact floating counter range are errors rather than wrapping arithmetic.

Compatibility tests compare every prefix at every contiguous partition of small fixtures with expected results generated from unmodified upstream revision `1f7fcb53aa20b592d8b3b97f75d9971a60eb06e9`. They compare marker maps and metadata semantically, and mean lengths by their exact floating-point bits. Normal tests read checked-in expectations and do not build or run the reference generator:

```sh
cargo test --locked --test read_sketch
cargo test --locked --no-default-features --lib --test read_sketch
```

To regenerate the expectations explicitly:

```sh
cargo run --locked --manifest-path tests/reference/Cargo.toml
```

The generator has a separate lockfile and a Git dependency on the unchanged upstream revision. Its native compression dependencies require a C compiler and CMake. Fixture metadata records the revision, compiler, target, and extraction backend. Generation checks that the inputs actually distinguish the duplicate policies and boundary cases.

This patch preserves upstream CPU dispatch. At the pinned revision, AVX2 extraction omits up to three trailing windows that scalar extraction includes. The fixtures avoid sampled omitted windows; passing them is not a general claim of scalar/AVX2 equivalence. The accumulator still follows whole-file extraction on the same backend. Overflow-limit behavior is guarded in code but is not exercised by the small scientific fixtures.

### Reusable profiling

Select `default-features = false, features = ["profile-core"]` to calculate a complete profile from an existing sketch. `Profiler` loads a standard `.syldb` once, validates every reference, and owns an explicitly sized local Rayon pool. It neither configures global logging nor initializes Rayon's global pool. `from_genomes` adopts already-loaded sketches through the same validation path.

```rust
use sylph::profile::{ProfileOptions, Profiler};

let profiler = Profiler::open_syldb("reference.syldb", ProfileOptions::default(), 3)?;
let profile = profiler.profile(sketcher.sketch())?;
```

The calculation includes filtering, shared-marker reassignment, coverage correction, and abundance normalization. Results own their names and numeric estimates, not the reference or sample sketches. Options remain fixed for the profiler's lifetime and accompany each result. Percentages and experimental estimated read counts have distinct `SequenceEstimate` variants. ANI values and intervals use fractional units; raw adjusted ANI can exceed one, as in upstream's internal estimates. The legacy CLI alone caps its displayed adjusted ANI at 100 percent.

`ProfileOptions::default()` matches the CLI's scientific profiling defaults. Controls include minimum ANI, contained/reference-marker thresholds, confidence intervals, unknown-fraction correction, experimental read-count estimation, and a read-identity override for corrected modes. Missing intervals remain absent; the result's options distinguish unavailable intervals from disabled calculation. Unknown correction retains upstream's automatic identity heuristic, whose short-read fallback depends on marker-map iteration order. Equal marker maps alone do not establish equal corrected profiles. Tied detections also retain upstream's unspecified parallel ordering.

Invalid configurations, reference/sample incompatibility, unsupported arithmetic ranges, and undefined estimates return errors without partial rows. Valid empty/no-hit samples return empty profiles. Standard databases are operator-provided local inputs; this is not an untrusted-upload parser. References must contain profiling markers and use a common supported k; a sample must have the same k and must not be sparser than any reference. Two-stage databases and raw/compressed input loading remain outside this interface.

Coverage sums retain upstream's `u32` limit, including each bootstrap resample. Inference vectors are conservatively limited to `i32::MAX` entries because legacy estimators use signed histogram counters. Unknown correction also checks the sample-wide count sum and covered-base multiplication for overflow. These checks reject unsupported arithmetic instead of changing the estimators' integer widths.

Calls are synchronous, retain no sample state, and are not cancellable. A service should run them outside its async executor threads and control concurrent work. Do not hold an ingestion lock throughout profiling unless pausing ingestion is intentional. This interface adds no refresh scheduler, snapshot copying, publication endpoint, or restart recovery.

Profile compatibility uses a separate synthetic cohort and a driver linked to the unchanged upstream revision above. It exercises ordered prefixes, strict shared-marker competition, filtering, multiple survivors, abundance normalization, and confidence intervals. Expectations are compared semantically at the original CLI's printed precision, not as proof of identical unrounded results. Normal tests consume the frozen data:

```sh
cargo test --locked --no-default-features --features profile-core --test profile
```

Explicit regeneration uses the separately locked reference crate. The driver calls upstream's public parser and command routines; it is not the original CLI executable and contains no adapted scientific implementation. Generation requires upstream's native build dependencies. Use a disposable directory for intermediate files:

```sh
cargo run --locked --manifest-path tests/reference/Cargo.toml --bin profile-reference -- generate "$ORACLE_DIR"
cargo run --locked --manifest-path tests/reference/Cargo.toml --bin profile-reference -- freeze "$ORACLE_DIR" tests/fixtures/profile
cargo run --locked --manifest-path tests/reference/Cargo.toml --bin profile-reference -- generate-options "$ORACLE_DIR"
cargo run --locked --manifest-path tests/reference/Cargo.toml --bin profile-reference -- freeze-options "$ORACLE_DIR" tests/fixtures/profile/options.yaml
```

## Citing sylph

Jim Shaw and Yun William Yu. Rapid species-level metagenome profiling and containment estimation with sylph (2024). Nature Biotechnology.
