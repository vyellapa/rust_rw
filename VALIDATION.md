# Validation results

Validated on Linux x86_64 on 2026-11-07.

## Official C++ reference

The C++ executable was extracted from the downloaded `delly-v2.8.1.sif` image. It reports DELLY 2.8.1, Boost 1.74.0 and HTSlib 1.21. The Rust executable uses HTSlib 1.19.

- SIF SHA-256: `b8e87db9769d0991ca4ed391618ee1b5c9e66d611ebf8e6de37f32614a007196`
- Extracted executable SHA-256: `5e6ea67e2d108f9bcb0ffa06099b13507e00e4d3a826dfee6cde007b3af88ce3`

All **38 comparison cases** had identical sample names/order, row counts, and complete rendered VCF records against this C++ reference. The suite checks every VCF record field without numeric tolerances. Header metadata, compressed-file bytes, index offsets and diagnostic dump ordering are excluded in records mode.

Coverage includes bundled-example discovery at 1/4/8 threads; 256 tied genotyping sites with read caps 0/1/5/250; a two-chromosome fixture with actual translocation mate pairs; additional caps 1/20/250; forced batch boundaries; dump mode; CRAM; multiple samples; and stdout VCF. The discovery fixture emits one row, tied-site fixtures emit 256 rows, and translocation fixtures emit 33 rows.

On the bundled example, even the decompressed BCF header and binary records matched exactly; compressed BCF and CSI bytes differed between builds. These are fixture results, not a guarantee for all datasets.

## Regression checks

- 10 Rust unit tests passed, including cap/parity behavior, ordered dump replay, cross-chromosome mates and FMA/fallback alignment equality.
- The optional over-4-GiB gzip regression passed full decompression length and CRC validation.
- Five Python harness tests passed, including failed callers with stale files, missing/invalid outputs, differing outputs, and deterministic C++ probe ordering under shuffled collection.
- The same 38 cases also matched BCF/CSI/dump bytes against the saved pre-optimization Rust binary.

The Apple clang/libc++ and Boost 1.92 byte-identity target documented by the original port was not revalidated against an Apple C++ build in this Linux environment.
