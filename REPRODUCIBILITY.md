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


## BioSystems ED3R1 source-locked GEO replay (17 October 2026)

The companion `external_validation/gse_fourcase_r2/` includes `download_sources.py`, `replay.py`, and `additional_tests.py`.
The corrected reproducibility contract limits the GSE245825 Benjamini--Hochberg family to the nine unique planned target/comparator combinations (3 targets CTHRC1/ASPN/WNT5A by 3 comparators COL1A1/POSTN/FSTL3). The prior ED3 replay enumerated 13 comparisons including redundant reversed pairs. The manuscript already reported the nine-comparison family, so ED3R1 brings code and output into agreement without changing Spearman rho or its raw P-value. The nine-comparison CTHRC1--COL1A1 q value is 3.455467599044716e-05.

For an internet-enabled run: `python external_validation/gse_fourcase_r2/download_sources.py --out source_data`; `python external_validation/gse_fourcase_r2/replay.py --sources source_data --out replay_results`; then copy `additional_tests.py` into `replay_results` and execute it there. The seven source gzip files are downloaded from official GEO only, SHA-256 checked, and excluded from manuscript release archives.

The ED3R1 reanalysis scripts and curated derived summaries were published with author authorization in the public project repository at `external_validation/gse_fourcase_r2/` (main commit `2fb9b47dd1a4d7c8e55e583ac5279e2ba4dd5eb8`). All three authors approved the ED3R1 BioSystems manuscript and submission, as reported by the corresponding author on 17 October 2026. The publication remains an author-controlled portal action; no submission tracking number has been issued.

## BioSystems pre-submission source synchronization (18 October 2026)

`main.tex` and `supplement.tex` match the author-approved ED3R1 scientific submission source; the added `external_validation/external_temporal_sensitivity.py` reproduces the post-hoc GSE152250 leave-one-culture-out and exact conditional treatment-by-time sensitivity from the already-released sample-value CSV. Its interaction reference value P=0.065 is **not** a verified biological delay.

For GSE298023, COL1A1/POSTN histology associations are previously published author calculations extracted from a selected supplement; no CTHRC1-histology association can be recovered from the absence of that gene in the selected list. The 58-biopsy GSE245825 correlations are gene coexpression, not histological fibrosis associations. This repository does not contain patient-identifiable data or original third-party GEO raw/normalized matrices.
