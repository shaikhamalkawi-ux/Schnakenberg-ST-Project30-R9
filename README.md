# From Cardiac Omics to Identifiable Dynamics

This repository accompanies the manuscript **From Cardiac Omics to Identifiable Dynamics: Evidence Boundaries and Measurement Decisions for a Delayed Schnakenberg Model**.

It contains the manuscript and Supplementary Information LaTeX sources, deterministic figure-generation code, source-traceability tables, and limited helper code. Third-party raw GEO/ArrayExpress data and the publisher-hosted PLOS S1 XLS are not redistributed; obtain them from their original sources.

## Primary empirical anchor

The primary human analysis uses 84 donor-linked primary human cardiac-fibroblast basal/TGF-beta1 pairs from GSE97358. CTHRC1 increased in 82/84 pairs, with a mean paired effect of +0.882 log2 units and donor-bootstrap 95% interval [0.801, 0.962].

## Independent external validation

GSE152250 provides an independent human primary cardiac-fibroblast TGF-beta experiment with 3 control and 3 treated culture replicates at 6, 24, and 48 h. The 24-h comparison was prespecified as the primary external-validation point because it matches GSE97358. CTHRC1 mean expression was 108.605 versus 171.000 FPKM at 24 h (log2 ratio of means +0.655), with all three treated replicates above all three controls. The direction remained positive at 48 h and was near-null at 6 h. This is within-study directional/time-course replication, not donor-level evidence and not a cross-study pooled effect.

## Repository contents

- `main.tex` - main manuscript source
- `supplement.tex` - Supplementary Information source
- `make_figures.py` - deterministic figure-generation code; run it before LaTeX compilation to create the figure assets
- `support/consolidated_numerical_results.csv` - machine-readable counterpart of the consolidated numerical-results table
- `external_validation/gse152250_validate_v2.py` - executable GSE152250 external-validation analysis
- `external_validation/results/` - derived GSE152250 sample values, summaries, mapping, input hashes, and audit report
- `support/saudi_target_status.csv` - Saudi-center source-status table
- `support/extract_saudi_targets.py` - deterministic helper for target extraction when the original PLOS S1 XLS is supplied
- `support/reference_audit.csv` - reference audit
- `support/reproducibility_note.md` - scope and limitations
- `DATA_SOURCES.md` and `REPRODUCIBILITY.md` - data provenance and reproducibility boundaries
- `BUILD.md` and `requirements.txt` - build instructions and Python dependencies

## Reproducibility boundaries

The repository does not claim independent re-extraction of the Saudi TGFB1/CTHRC1 rows from the original publisher XLS, a newly rerun donor-blocked negative-binomial coefficient for GSE97358, or reconstruction of the unavailable four-family timing-model specification. No myocardial delay, physical diffusion coefficient, or cardiac Turing/Hopf/Turing-Hopf certificate is asserted from the current evidence.


## BioSystems ED3R1 source-locked GEO replay (17 October 2026)

The companion `external_validation/gse_fourcase_r2/` includes `download_sources.py`, `replay.py`, and `additional_tests.py`.
The corrected reproducibility contract limits the GSE245825 Benjamini--Hochberg family to the nine unique planned target/comparator combinations (3 targets CTHRC1/ASPN/WNT5A by 3 comparators COL1A1/POSTN/FSTL3). The prior ED3 replay enumerated 13 comparisons including redundant reversed pairs. The manuscript already reported the nine-comparison family, so ED3R1 brings code and output into agreement without changing Spearman rho or its raw P-value. The nine-comparison CTHRC1--COL1A1 q value is 3.455467599044716e-05.

For an internet-enabled run: `python external_validation/gse_fourcase_r2/download_sources.py --out source_data`; `python external_validation/gse_fourcase_r2/replay.py --sources source_data --out replay_results`; then copy `additional_tests.py` into `replay_results` and execute it there. The seven source gzip files are downloaded from official GEO only, SHA-256 checked, and excluded from manuscript release archives.

The ED3R1 reanalysis scripts and curated derived summaries were published with author authorization in the public project repository at `external_validation/gse_fourcase_r2/` (main commit `2fb9b47dd1a4d7c8e55e583ac5279e2ba4dd5eb8`). All three authors approved the ED3R1 BioSystems manuscript and submission, as reported by the corresponding author on 17 October 2026. The publication remains an author-controlled portal action; no submission tracking number has been issued.

## BioSystems ED3R1 manuscript synchronized (18 October 2026)

The named-author `main.tex` and `supplement.tex` at the default branch are the author-confirmed **BioSystems** pre-submission scientific source (22-page manuscript; 12-page Supplement in the verified local build). Ghassan Malkawi is first and corresponding author; Ahmed Elsayed and Hassan Y. Alfifi follow. This is **not** an accepted or already submitted manuscript.

Curated source-replayed GEO extension code and selected **derived** summary tables: `external_validation/gse_fourcase_r2/`. Three contexts were admitted: GSE202714 human fibroblast TGF-beta1 response across cell sources, GSE246298 experimental organoid intervention, and GSE245825 within-DCM transcript coexpression. GSE234715 remains a sample-mapping HOLD; no quantitative injury conclusion is drawn.

The GSE298023 tissue-histology comparison reads selected **published supplementary tables only**; it does **not** produce a new CTHRC1–histological fibrosis coefficient. Do not confuse transcript coexpression with collagen-fraction correlation or clinical prediction. Original GEO source files and publisher supplements must be fetched from their original hosts.

The original reference audit is historical; the current key/DOI presence audit for the BioSystems source is `support/reference_audit_ED3R1.csv`. This audit does not independently validate each DOI. This repository is not a Zenodo DOI record; cite a project Zenodo DOI only if an existing, matching published record is independently confirmed.
