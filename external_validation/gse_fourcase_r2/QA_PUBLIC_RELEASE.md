# Project 30 — ED3R1 public release QA (17 October 2026)

The code and summary outputs were sourced from the author-confirmed ED3R1 Overleaf/code package. No main manuscript or Supplementary TeX source has been changed during this public release.

**Validated local source files:** 7 official NCBI GEO gzips. Their filenames, exact byte sizes, official HTTPS source URLs and SHA-256 digests appear in `Source_SHA256_Audit.csv` and `replay.py`. No source gzips are in GitHub.

**Reproduction check:** the original `replay.py` ran with `--sources` pointing to a private/local SHA-matched source directory and an isolated output directory, followed by `additional_tests.py`. Twelve deterministic calculated files were compared against the validated ED3R1 package (including the gene-level working outputs kept private): **12/12 exact byte matches**. No external Colab session or data-sharing account was required.

**Public content filter:** Only summary-level model-free gene comparisons, source annotations and HOLD documents are published. Source binary files (`*.gz`, `*.xlsx`, `*.xls`), full gene matrices and per-sample targets (`GSE202714_Gene_Values.csv`, `GSE245825_Gene_Values.csv`, `GSE234715_GeneValues_NO_GROUP_COMPARISON.csv`) are deliberately absent from the public tree. The scripts may regenerate these omitted tables locally after downloading official inputs.

**Correction:** The ED3R1 code tests nine nonredundant DCM gene–comparator pairs, not 13. The BH correction for CTHRC1/COL1A1 yields q=3.455467599044716e-05. No claims of gene–histology correlation, CTHRC1-specific treatment efficacy, myocardial transport/delay, biological onset, or identified reaction topology arise from the four GEO analyses.

**Author confirmation:** The corresponding author communicated author approval for the BioSystems ED3R1 manuscript and subsequently authorized public posting of code and *derived* results in this existing GitHub repository. Journal portal submission is a separate action.

**Limit:** The prior repository `main.tex` and `supplement.tex` do not become the BioSystems ED3R1 manuscript by publishing this subdirectory. No claim of all-manuscript computational reproduction is made.
