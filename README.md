# From Cardiac Omics to Identifiable Dynamics

This repository accompanies the manuscript **From Cardiac Omics to Identifiable Dynamics: Evidence Boundaries and Measurement Decisions for a Delayed Schnakenberg Model**.

It contains the manuscript and Supplementary Information LaTeX sources, deterministic figure-generation code, final figure assets, source-traceability tables, and limited helper code. Third-party raw GEO/ArrayExpress data and the publisher-hosted PLOS S1 XLS are not redistributed; obtain them from their original sources.

## Primary empirical anchor

The primary human analysis uses 84 donor-linked primary human cardiac-fibroblast basal/TGF-beta1 pairs from GSE97358. CTHRC1 increased in 82/84 pairs, with a mean paired effect of +0.882 log2 units and donor-bootstrap 95% interval [0.801, 0.962].

## Repository contents

- `main.tex` - main manuscript source
- `supplement.tex` - Supplementary Information source
- `figures/` - final PDF figure assets used by LaTeX
- `make_figures.py` - deterministic figure-generation code
- `support/consolidated_numerical_results.csv` - machine-readable counterpart of the consolidated numerical-results table
- `support/saudi_target_status.csv` - Saudi-center source-status table
- `support/extract_saudi_targets.py` - deterministic helper for target extraction when the original PLOS S1 XLS is supplied
- `support/reference_audit.csv` - reference audit
- `support/reproducibility_note.md` - scope and limitations
- `DATA_SOURCES.md` and `REPRODUCIBILITY.md` - data provenance and reproducibility boundaries
- `BUILD.md` and `requirements.txt` - build instructions and Python dependencies

## Reproducibility boundaries

The repository does not claim independent re-extraction of the Saudi TGFB1/CTHRC1 rows from the original publisher XLS, a newly rerun donor-blocked negative-binomial coefficient for GSE97358, or reconstruction of the unavailable four-family timing-model specification. No myocardial delay, physical diffusion coefficient, or cardiac Turing/Hopf/Turing-Hopf certificate is asserted from the current evidence.
