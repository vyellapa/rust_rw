# delly-rs

A parallel Rust port of [Delly](https://github.com/dellytools/delly) v2.8.1 short-read SV
calling (`delly sr`). For the same inputs and options, its BCF output (and `.csi` index,
and the optional `-d` gzipped read dump) is **byte-for-byte identical** to the C++ Delly
binary built on the same machine.

```
delly sr -g ref.fa -o calls.bcf sample1.bam sample2.bam ...
```

All `delly sr` options are supported: discovery (`-t -q -r -s -c -z -m -n -p`), genotyping
(`-v -u -d -a --sex`), `-x` exclude files (plain or gzipped), `-h` threads, CRAM input and
VCF-to-stdout when `-o` is omitted.

## Build

Requirements: Rust (stable), a C/C++ compiler, and htslib (the same build the reference
Delly links against; Homebrew's is used by default).

```
cargo build --release              # binary: target/release/delly
HTSLIB_PREFIX=/usr/local cargo build --release   # other htslib location
HTSLIB_PREFIX=/usr cargo build --release         # Linux system htslib
```

htslib (BAM/CRAM/BCF I/O) is called through FFI, and edlib is compiled from Delly's own
bundled `edlib.cpp`. Everything Delly itself does is reimplemented in Rust.

## What "byte-identical" required

Delly's output depends on behaviour that comes from the C++ standard library and Boost, not
from Delly's own code. These pieces are reproduced exactly:

| C++ behaviour | Where it matters | Rust |
|---|---|---|
| libc++ `std::sort` (introsort) on comparators with ties | SV, paired-end, breakpoint and span sorting | `src/cxx_sort.rs`, a faithful port, verified on 23,994 cases |
| `std::unordered_set<std::string>` iteration order + `std::hash<std::string>` (CityHash64) | Order of reads fed into the consensus MSA | `src/cxx_hash.rs`, verified on 20k insertion sequences / 2M hashes |
| Boost 1.92 `hash_combine` / `hash_range` | Read/pair identifiers, which feed sort order and fragment de-duplication | `src/cxx_hash.rs` |
| clang FMA contraction (`-ffp-contract=on`) | MSA profile score, genotype likelihoods, consensus entropy | `f32/f64::mul_add` in exactly those spots |
| compile-time folded `log10(2)` | genotype likelihoods | exact bit pattern |
| `boost::iostreams::gzip_compressor` | `-d` read dump | system zlib with Boost's header/footer |
| signed `char`, `unsigned`/`size_t` wrap-around and narrowing | hashes, coordinates | explicit wrapping arithmetic |

The output also contains `##fileDate`, so outputs only match when produced on the same
calendar day.

## Parallelism

C++ Delly parallelises per sample for scanning, and per SV for assembly and alignment batches.
delly-rs adds:

* **Scanning:** one task per (sample, chromosome). Mate pairs that span chromosomes are
  recorded as ordered events and replayed in chromosome order afterwards.
* **Split-read collection:** one task per chromosome. Reads are then inserted into the
  libc++-ordered sets in the original order.
* **Genotyping:** one task per (sample, chromosome). Each task emits an ordered event log:
  junction-read candidates, spanning pairs, and translocation mates. Per-sample replay admits
  candidates in the original C++ batches, then aligns independent SVs in parallel. When
  there are few SVs, read waves within an SV also run in parallel. Each wave is at most the
  remaining read capacity, so a reached cap skips the remaining alignments. Ordered merging
  retains reference subsampling, cross-chromosome state and SR/PE dump interleaving.
* **Probe generation and paired-end clustering:** parallel per SV / per SV type.
* **Profile alignment on x86:** FMA-capable CPUs use a hardware-specialized kernel selected
  once per alignment. Other CPUs use the portable kernel. Both preserve fused rounding;
  no global native-CPU flags or fast-math are required.

The output is independent of `-h`. Benchmark on representative data: chromosome count,
read caps, coverage and number of SVs affect which phases can use all workers.

## Verifying

`scripts/compare.sh <cpp-delly> <rust-delly> <outdir> <name> [delly sr args]` runs both
binaries and compares the BCF, `.csi` and dump byte for byte. In development the following
were all identical: the bundled example, two simulated data sets (1–3 samples, every SV
type including BND and inline-sequence INS/DEL), thread counts 1/4/8/10/16, `-t` filters,
exclude files, parameter sweeps, read caps (`-a 5/20/30`), `--sex none/male/file/auto`
with a haploid chrX, `-v` genotyping, `-d` dumps, CRAM input and stdout VCF.

The comparison script uses fresh temporary outputs and returns nonzero for failed callers,
missing/invalid outputs or byte differences. Final logs remain available after a failed run.

```
HTSLIB_PREFIX=/usr cargo test --release
python3 scripts/test_compare.py
python3 scripts/validate.py /path/to/reference-delly
# Compare calls across builds with different compression/HTSlib versions:
python3 scripts/validate.py /path/to/reference-delly --comparison records
# Optional large-input regression (processes over 4 GiB and requires gzip):
HTSLIB_PREFIX=/usr cargo test --release input_over_four_gib_is_not_truncated -- --ignored
```

`validate.py` requires samtools and bcftools. It generates tied sites, a two-chromosome
fixture with real translocation mates, CRAM and multi-sample inputs; checks read caps,
forced batch boundaries, dumps and stdout; and retains outputs/logs for inspection.
It can also compare a saved Rust binary against the optimized build. For C++ equivalence,
use the documented reference toolchain and the bundled reference ordering fix below.
`--comparison records` instead checks row counts, sample names/order and every field of
the rendered VCF records exactly, while ignoring compression, header metadata, index
bytes and diagnostic read-dump ordering. This mode can validate the downloaded Linux
DELLY container without requiring its Boost and HTSlib builds to match the Rust build.

## Scope and caveats

* Only `sr` is ported. `lr`, `asm`, `cnv`, `merge` and `filter` are not.
* Identity holds against a reference built with **Apple clang/libc++ and Boost 1.92**. A
  Delly built with GCC/libstdc++, or with another Boost, sorts ties and iterates hash sets
  differently, so it can produce different (equally valid) output on some inputs.
* The bundled C++ reference restores serial probe collection order before its original
  position-only sort (`delly/src/probeorder.h`). This preserves its original single-threaded
  tie permutation. An unpatched multi-threaded reference can produce different read-dump
  ordering at tied breakpoints; rebuilding it with this fix makes the reference deterministic.
* For speed, delly-rs opens one BAM handle per (sample, chromosome) task.
* `DELLY_TIMING=1` prints per-phase timings. `DELLY_BATCH_BASE` exists only to test the
  genotyping batch-boundary emulation.
