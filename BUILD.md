# Build instructions

Third-party raw datasets are not redistributed.

## Regenerate figures

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python make_figures.py
```

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
