# Build instructions

Third-party raw datasets are not redistributed.

## Regenerate figures

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python make_figures.py
```

## Re-run the GSE152250 external validation

Download the two official NCBI GEO files without redistributing them:

```bash
mkdir -p external_validation/input external_validation/output
curl -L --fail -o external_validation/input/GSE152250_Fibroblasts_FPKM.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152250/suppl/GSE152250_Fibroblasts_FPKM.txt.gz
curl -L --fail -o external_validation/input/GSE152250_family.soft.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152250/soft/GSE152250_family.soft.gz
sha256sum external_validation/input/*
python external_validation/gse152250_validate_v2.py \
  --fpkm external_validation/input/GSE152250_Fibroblasts_FPKM.txt.gz \
  --soft external_validation/input/GSE152250_family.soft.gz \
  --outdir external_validation/output
```

The committed `external_validation/results/` files record the accepted derived result and official-input hashes. The raw third-party GEO files are not committed.

## Compile the main manuscript

```bash
pdflatex main.tex
pdflatex main.tex
```

## Compile the Supplementary Information

```bash
pdflatex supplement.tex
pdflatex supplement.tex
```

The manuscript uses explicit `thebibliography` sections; BibTeX is not required.
