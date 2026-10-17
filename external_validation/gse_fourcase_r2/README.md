# BioSystems ED3R1: four GEO datasets (public code and derived summaries)

**Project 30 — cardiac omics, CTHRC1, and evidence boundaries for delayed Schnakenberg-type models.**
The authors authorized public release of the computational code and derived summary outputs on 17 October 2026. **This is not a record of journal submission or acceptance.**

This directory contains the **independently replayed, ED3R1-corrected** external-data analysis from four official NCBI GEO accessions. Only aggregate/summary results, code, schema annotations, and hold records are distributed here. Original third-party GEO matrices (including `.gz` and `.xlsx`) and sample-level expression/count matrices are **not included**.

The repository's root `main.tex`/`supplement.tex` may be an earlier manuscript build; the **ED3R1 BioSystems manuscript and supplement are not replaced by this code release**.

## Included datasets and limits

| GEO | Analysis | Status / scientific limit |
|---|---|---|
| [GSE246298](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246298) | Three cardiac organoid lines, control vs fibrotic vehicle and FGFR inhibitor | **Exploratory descriptive** within-line expression contrasts, not proven drug efficacy or a specific causal CTHRC1 effect |
| [GSE202714](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202714) | 3 fibroblast sources × 3 cultures per arm; TGF-β1 vs control | **Exploratory** source-stratified culture responses; no donor-independent replication claim |
| [GSE245825](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245825) | 58 DCM biopsies, 9 declared coexpression tests | **Coexpression only**; cannot measure CTHRC1 vs tissue collagen without sample-keyed histology |
| [GSE234715](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234715) | 12 engineered-heart-tissue samples | **HOLD** injury/control contrasts: native matrix columns not conclusively mapped to source GSM identifiers |

### Checked figures (not clinical validation)

- GSE246298 CTHRC1 mean contrasts on `log2(source-provided normalized counts + 1)`: fibrotic vehicle vs control **+1.1886978289578805**, inhibitor vs fibrotic vehicle **−1.2330769881611643** across three source lines (descriptive only).
- GSE202714 CTHRC1 source-specific contrasts: **+0.8525371140882392** (hiPSC-CF), plus positive contrasts in aHCF and fHCF; see source-stratified exact-permutation summary. Cultures, not donors, are the units analyzed.
- GSE245825 CTHRC1–COL1A1: **Spearman rho 0.5412065103915417** (58 DCM biopsies); BH-adjusted q across the **9 predeclared comparisons** is **3.455467599044716×10⁻⁵**. It is **not a quantitative collagen histology association**.

The delayed Schnakenberg model's cardiac state mapping, delay, diffusion matrix, Turing/Hopf classification, and clinical decision utility **remain unvalidated** by these data.

## Reproduce from original GEO files (not included)

Use Python 3.11+ with NumPy and SciPy. Run from repository root:

```bash
python -m pip install -r external_validation/gse_fourcase_r2/requirements.txt
python external_validation/gse_fourcase_r2/download_sources.py --out source_data
python external_validation/gse_fourcase_r2/replay.py --sources source_data --out replay_results
cp external_validation/gse_fourcase_r2/additional_tests.py replay_results/
(cd replay_results && python additional_tests.py)
```

Each of the seven official input `.gz` files must pass the exact SHA-256 checks in `replay.py` before the scientific script runs. Downloaded raw data and generated sample-level values are intended to stay **local** and are ignored by this subdirectory's `.gitignore`. The download script can also use a previously downloaded, author-held folder (`--from-dir`). The public code can generate detailed sample-level audit tables locally, but only the summary result files are hosted here.

`replay.py`'s default `--out` refers to an older analysis-environment path; **always pass `--out replay_results`** as above. `additional_tests.py` must run from a directory populated by `replay.py`; its sample-level input tables are deliberately not republished here.

## ED3R1 correction, reproducibility and provenance

The earlier exploratory ED3 helper enumerated 13 coexpression tests because it included redundant reverse directions. **ED3R1 corrected it to 9 unique target/comparator pairs** (targets CTHRC1, ASPN, WNT5A × comparators COL1A1, POSTN, FSTL3), consistent with the manuscript. Rho and unadjusted p were unchanged; the BH correction is tied to these 9 tests. Code and results were replayed on the seven hash-checked original GEO source files and **12/12 deterministic output files matched byte-for-byte**. An independent sensitivity script computes source-stratified exact-label permutations and 3,000-bootstrap within-DCM coexpression intervals; see `QA_PUBLIC_RELEASE.md` for the exact audit scope and limits.

The public release includes `Source_SHA256_Audit.csv`, `Ensembl_NCBI_Target_Map.csv`, a SHA manifest for the published files, and status/hold documents. It does **not** contain raw GEO records, source `.gz` files, any 58-biopsy fibrosis measurement table, any patient clinical record, or unpublished original instrument data.

**Rights and citation:** original data belong to their respective depositors and must be obtained from NCBI GEO with original attribution. No third-party data license is asserted by this repository. Public availability of code does not itself grant a particular open-source license; check repository license/author statements before reuse. Cite the original source studies and accession IDs in publications.
