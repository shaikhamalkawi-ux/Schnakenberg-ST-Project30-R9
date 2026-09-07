# From Cardiac Omics to Identifiable Dynamics

This private repository accompanies the manuscript **“From Cardiac Omics to Identifiable Dynamics: Evidence Boundaries and Measurement Decisions for a Delayed Schnakenberg Model.”**

The repository is organized to support transparent review of the manuscript, figures, source-traceability materials, and limited deterministic helper code. It does **not** redistribute third-party raw datasets. Public data should be obtained from the original repositories and publications listed below.

## Repository status

This repository is currently **PRIVATE** and is being used as the working deposit for Project 30 / R9. It is not yet the public archival record.

The manuscript currently contains temporary repository placeholders. The final public GitHub URL and Zenodo DOI/record URL should be inserted only after the author approves public release.

## Primary result retained in the manuscript

The primary human analysis uses 84 donor-linked primary human cardiac-fibroblast basal/TGF-beta1 pairs from GSE97358. CTHRC1 increased in 82/84 pairs; the mean paired effect is +0.882 log2 units with donor-bootstrap 95% interval [0.801, 0.962]. The repository preserves the manuscript-level evidence boundaries around this result.

## Planned repository structure

- `main.tex` - main manuscript source
- `supplement.tex` - Supplementary Material source
- `figures/` - final manuscript figure assets
- `make_figures.py` - deterministic figure-generation code used for the included figures
- `support/` - source-traceability and reference-audit materials
- `manuscript/` - main, supplement, and combined PDFs
- `DATA_SOURCES.md` - public accession and source list
- `REPRODUCIBILITY.md` - reproducibility scope and explicit open evidence items
- `requirements.txt` - Python dependencies for figure generation
- `CITATION.cff.template` - citation metadata template to complete when author identities are released
- `LICENSE_SELECTION_REQUIRED.md` - reminder that repository licensing remains an author decision
- `MANIFEST_SHA256.txt` - SHA-256 manifest for repository contents

## Public data sources

The manuscript uses the following GEO accessions: GSE265828, GSE267256, GSE261428, GSE123018, GSE96991, GSE96975, GSE97358, GSE225336, and GSE132143. Saudi-center transcriptome context comes from Colak et al. (2009, 2016) and ArrayExpress E-TABM-480. The PLOS ONE 2016 S1 table is identified by DOI `10.1371/journal.pone.0162669.s005`.

## Important reproducibility boundaries

1. The original publisher-hosted PLOS S1 XLS is not redistributed here; Saudi TGFB1/CTHRC1 target-level statements remain published-source contextual reports unless the XLS is independently supplied and re-extracted.
2. The original GSE97358 raw-count object is not bundled; no new donor-blocked negative-binomial sensitivity is claimed in this release.
3. The original four-family timing-model specification file is unavailable; that comparison remains a supporting diagnostic only.
4. No myocardial delay, physical diffusion coefficient, or cardiac Turing/Hopf/Turing-Hopf certificate is asserted from the current public evidence.

## Release status

Keep this repository private until the manuscript/repository release decision is finalized. Zenodo should likewise remain a draft until publication is explicitly approved.
