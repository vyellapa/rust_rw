# rust_rw

A parallel Rust rewrite of DELLY 2.8.1 short-read structural-variant discovery and genotyping (`delly sr`).

- `delly-rs/`: Rust implementation, FFI shims, validation scripts and regression tests.
- `delly/`: bundled C++ reference source and example inputs. Upstream snapshot: `dellytools/delly` commit `9e568441b9627003291c07b0135dfaf6eb1bdf73`, with deterministic probe ordering added locally.

The acceptance criterion is identical row counts, sample order and complete VCF record fields, including coordinates, alleles, genotypes, quality and INFO/FORMAT values. Compressed BCF and CSI bytes may differ across HTSlib/compiler builds.

## Build

Install stable Rust, a C/C++ compiler and HTSlib development libraries. On Linux with system HTSlib:

```bash
cd delly-rs
HTSLIB_PREFIX=/usr cargo build --release
./target/release/delly --version
```

On macOS, the build defaults to the Homebrew HTSlib prefix `/opt/homebrew`. See [the Rust README](delly-rs/README.md) for compatibility details and options.

## Validate

```bash
cd delly-rs
HTSLIB_PREFIX=/usr cargo test --release
python3 scripts/test_compare.py
python3 scripts/validate.py /path/to/cpp-delly --comparison records
```

The Python checks require samtools, bcftools and a C++ compiler. `validate.py` retains outputs and logs and reports nonzero on mismatches. Its default `--comparison bytes` also compares compressed BCF, CSI and optional dump bytes when matching toolchains are available.

See [VALIDATION.md](VALIDATION.md) for the completed comparison against the downloaded official C++ executable, plus implementation and test details in [delly-rs/README.md](delly-rs/README.md).

## C++ reference build

The original HTSlib submodule pin is preserved under `delly/src/htslib`. Initialize it with `git submodule update --init --recursive`, or build against installed HTSlib with the Makefile overrides documented in `delly/Makefile`. The C++ reference also requires Boost development libraries.

The local probe-ordering change preserves the original single-threaded collection/sort order while making parallel collection deterministic. Licensing and upstream notices are retained in `delly/LICENSE`, `delly/CREDITS`, and the original source files.
