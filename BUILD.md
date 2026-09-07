# Build instructions

The repository is intentionally source-first. Third-party raw datasets are not redistributed.

## 1. Create the manuscript figures

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python make_figures.py
```

This creates the PDF and PNG figure assets under `figures/`.

## 2. Compile the main manuscript

```bash
pdflatex main.tex
pdflatex main.tex
```

## 3. Compile the Supplementary Information

```bash
pdflatex supplement.tex
pdflatex supplement.tex
```

The manuscript uses an explicit `thebibliography` section; BibTeX is not required.

## Reproducibility boundary

Figure generation and LaTeX compilation are deterministic from the files in this private repository. Historical analyses for which the original raw-count object, PLOS S1 XLS binary, or original timing-model specification are not present are not represented as independently rerun calculations; see `REPRODUCIBILITY.md` and `support/reproducibility_note.md`.
