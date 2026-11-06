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
```

To regenerate the expectations explicitly:

```sh
cargo run --locked --manifest-path tests/reference/Cargo.toml
```

The generator has a separate lockfile and a Git dependency on the unchanged upstream revision. Its native compression dependencies require a C compiler and CMake. Fixture metadata records the revision, compiler, target, and extraction backend. Generation checks that the inputs actually distinguish the duplicate policies and boundary cases.

This patch preserves upstream CPU dispatch. At the pinned revision, AVX2 extraction omits up to three trailing windows that scalar extraction includes. The fixtures avoid sampled omitted windows; passing them is not a general claim of scalar/AVX2 equivalence. The accumulator still follows whole-file extraction on the same backend. Overflow-limit behavior is guarded in code but is not exercised by the small scientific fixtures.

## Citing sylph

Jim Shaw and Yun William Yu. Rapid species-level metagenome profiling and containment estimation with sylph (2024). Nature Biotechnology.
