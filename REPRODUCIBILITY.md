# Reproducibility scope

This repository supports source verification, LaTeX compilation, figure regeneration, source traceability, and selected deterministic helper operations. It is not described as a complete end-to-end reproduction of historical calculations for which the original analysis object is unavailable.

## Included

- Main and supplementary LaTeX sources
- Deterministic figure-generation code for all manuscript and supplementary figures
- Consolidated numerical-results CSV
- Executable GSE152250 external-validation script
- Derived GSE152250 sample-level CTHRC1 values, time-point summary, exact sample mapping, schema record, and official-input SHA-256 hashes
- Saudi-center target-status CSV
- Saudi target-extraction helper
- Reference audit and reproducibility note

## Explicitly not included / not claimed

- Original publisher-hosted PLOS S1 XLS binary
- Original GSE97358 raw-count object required for a new donor-blocked negative-binomial sensitivity
- Original four-family timing-model specification file
- Third-party raw GEO/ArrayExpress datasets, including the GSE152250 FPKM and SOFT source files (downloaded from NCBI at execution time and hash-locked, but not redistributed)

## Consequence for claims

The manuscript does not claim independent re-extraction of the Saudi TGFB1/CTHRC1 rows, a newly rerun donor-blocked negative-binomial coefficient, or a newly reconstructed four-family timing-model comparison. GSE152250 supports independent within-study directional/time-course replication of the 24-h CTHRC1 response; its three culture replicates per arm are not treated as donor-level evidence and its FPKM-scale magnitude is not pooled numerically with GSE97358.
