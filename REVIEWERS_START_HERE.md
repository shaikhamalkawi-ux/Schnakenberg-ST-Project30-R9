# Reviewer access and reproducibility — BioSystems ED3R1

**Manuscript:** From Cardiac Omics to Identifiable Dynamics: Evidence Boundaries and Measurement Decisions for a Delayed Schnakenberg Model  
**Submission preparation:** 18 October 2026  
**Authors (in order):** Ghassan Malkawi (first and corresponding), Ahmed Elsayed, Hassan Y. Alfifi  
**Public source:** https://github.com/shaikhamalkawi-ux/Schnakenberg-ST-Project30-R9

This publicly readable repository is linked directly in the manuscript's **Code availability** section (page 18 of the 22-page BioSystems submission PDF). It is an author-approved pre-submission research record, **not** a statement that the manuscript has been accepted or submitted.

## Where to find things

| Item | Location | What it supports |
| --- | --- | --- |
| Main manuscript and Supplement source | [main.tex](main.tex), [supplement.tex](supplement.tex) | Exact BioSystems ED3R1 scientific narrative, results and inference limits |
| Figures | [make_figures.py](make_figures.py), generated `figures/` at build time | Recreate the publication figures |
| Primary human evidence | [support/consolidated_numerical_results.csv](support/consolidated_numerical_results.csv) | Previously established, donor-linked 84-pair summary |
| Independent 6/24/48-h validation | [external_validation/](external_validation/) | GSE152250 mapping, source checks, code and derived evidence |
| Post-hoc temporal and culture-deletion audit | [external_validation/external_temporal_sensitivity.py](external_validation/external_temporal_sensitivity.py) and [derived outputs](external_validation/results/) | Descriptive sensitivity; **not** identification of onset, hard delay or donor effects |
| Three GEO exploratory extensions | [external_validation/gse_fourcase_r2/](external_validation/gse_fourcase_r2/) | GSE202714, GSE246298 and GSE245825 independent clean replay, statistical tables, source hashes |
| Contextual Saudi-center evidence | [support/saudi_target_status.csv](support/saudi_target_status.csv) | Source-status and transferability limits, **not** Saudi-population validation |
| Reproduction boundaries | [REPRODUCIBILITY.md](REPRODUCIBILITY.md), [DATA_SOURCES.md](DATA_SOURCES.md) | Exactly what was, and was not, independently re-estimated |
| Authoritative publication metadata | [README.md](README.md), [BUILD.md](BUILD.md) | Dependencies, source build instructions and paper identity |

## Build the manuscript

Use Python 3 and a full LaTeX installation with `latexmk` and `elsarticle`.

```bash
python -m pip install -r requirements.txt
python make_figures.py
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error supplement.tex
```

The GitHub Actions [clean-build workflow](.github/workflows/clean-build.yml) automates the source build and produces verification artifacts. The pre-submission manuscript is 22 pages and the Supplement is 12 pages in the verified build.

## Replay the GEO exploratory extensions

The code retrieves the **seven published GEO gzip source files** from the original NCBI archive and validates their expected SHA-256 digests; source binaries are deliberately **not redistributed** in this repository.

```bash
python -m pip install -r external_validation/gse_fourcase_r2/requirements.txt
python external_validation/gse_fourcase_r2/download_sources.py --out source_data
python external_validation/gse_fourcase_r2/replay.py --sources source_data --out replay_results
cp external_validation/gse_fourcase_r2/additional_tests.py replay_results/
(cd replay_results && python additional_tests.py)
```

Consult the [four-case README](external_validation/gse_fourcase_r2/README.md) and [source SHA register](external_validation/gse_fourcase_r2/Source_SHA256_Audit.csv). The manuscript's exploratory GSE245825 Benjamini-Hochberg family contains **nine unique predeclared comparisons**; the earlier 13-row duplicate-containing draft was corrected before this release. The submitted GSE245825 association is **within-cohort gene coexpression**, not histological fibrosis prediction.

## GSE298023 published human-tissue context

The limited audit extracts published selected-list entries from Youness et al. (2025), Additional file 2, Tables S5/S6/S9 (DOI [10.1186/s13073-025-01580-z](https://doi.org/10.1186/s13073-025-01580-z)). The original publisher-hosted supplementary workbook is not redistributed. The audit's exact extraction code and derived cross-checks are also included in the editable source bundle supplied to the journal. Absence from the selected list does not mean zero association.

## Interpretation and availability limits

- GSE97358: 84 donor-linked fibroblast basal/TGF-β1 pairs are the primary human perturbation anchor; this repository is not a re-execution of any unavailable historical raw-count object.
- GSE152250: independent **culture-level** directional/temporal confirmation, not donor-level pooling.
- GSE202714: source-stratified **culture** responses, not independent-donor validation.
- GSE246298: experimental **organoid transcript response**, not clinical drug efficacy.
- GSE245825: correlation of gene expression within DCM biopsies, **not** gene-to-collagen histology linkage.
- GSE234715: sample mapping remains **HOLD**; no injury contrast is claimed.
- GSE298023: selected published histology-related tables only; no independently computed CTHRC1 histology coefficient.
- No direct cardiac identification of the Schnakenberg states, myocardial transport/diffusion, hard delay `tau`, or Turing/Hopf/Turing–Hopf stability certificate is claimed.

This public repository provides the paper's source code, traceable derived summaries and reproducibility guidance. Original third-party datasets remain with their respective GEO/publisher hosts. No project-specific Zenodo DOI is asserted.
