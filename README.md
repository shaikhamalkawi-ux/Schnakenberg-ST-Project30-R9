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
