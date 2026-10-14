#!/usr/bin/env python3
"""
Independent external validation of the manuscript's primary CTHRC1 direction
using GEO GSE152250 human cardiac fibroblasts.

Prespecified validation target:
- Primary external replication time point: 24 h (matches GSE97358).
- Secondary temporal context: 6 h and 48 h.
- Unit: independent culture replicate within GSE152250; no donor pairing is assumed.
- Primary validation criterion: mean CTHRC1 expression under TGF-beta1 is greater
  than control at 24 h. Small-n exact permutation statistics are reported as
  descriptive support, not used to redefine the criterion post hoc.
- Cross-study effect sizes are not pooled with GSE97358.
"""
from __future__ import annotations

import argparse
import gzip
import itertools
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


TARGET_GENES = ["CTHRC1"]\nGENE_IDS = {"CTHRC1": "ENSG00000164932"}
TIMEPOINTS = [6, 24, 48]
PSEUDOCOUNT = 0.5


def parse_soft(path: Path) -> pd.DataFrame:
    rows = []
    current = None
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current:
                    rows.append(current)
                current = {"gsm": line.split("=", 1)[1].strip(), "descriptions": []}
            elif current is not None and line.startswith("!Sample_title = "):
                current["title"] = line.split("=", 1)[1].strip()
            elif current is not None and line.startswith("!Sample_description = "):
                current["descriptions"].append(line.split("=", 1)[1].strip())
    if current:
        rows.append(current)

    out = []
    for r in rows:
        title = r.get("title", "")
        if not title.lower().startswith("cardiac "):
            continue
        m = re.search(r"cardiac\s+(\d+)h\s+(CTRL|TGFb1)\s*(\d+)?", title, flags=re.I)
        if not m:
            continue
        out.append(
            {
                "gsm": r["gsm"],
                "title": title,
                "time_h": int(m.group(1)),
                "group": "TGFb1" if m.group(2).lower().startswith("tgf") else "CTRL",
                "replicate": int(m.group(3)) if m.group(3) else None,
                "description": " | ".join(r["descriptions"]),
            }
        )
    df = pd.DataFrame(out)
    if len(df) != 18:
        raise RuntimeError(f"Expected 18 cardiac samples in SOFT metadata, found {len(df)}")
    return df.sort_values(["time_h", "group", "replicate", "gsm"]).reset_index(drop=True)


def norm_token(x: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(x).lower())


def load_expression(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    if df.shape[1] < 19:
        raise RuntimeError(f"Unexpected FPKM table shape: {df.shape}")
    return df


def find_symbol_column(df: pd.DataFrame) -> str:
    """Locate either a gene-symbol column or an Ensembl-gene-ID column."""
    ensembl = GENE_IDS["CTHRC1"].upper()
    exact = []
    contains = []
    for col in df.columns:
        vals = df[col].astype(str).str.strip()
        upper = vals.str.upper()
        base = upper.str.replace(r"\.\d+$", "", regex=True)
        if upper.eq("CTHRC1").any() or base.eq(ensembl).any():
            exact.append(col)
        elif upper.str.contains(r"\bCTHRC1\b", regex=True).any():
            contains.append(col)
    if exact:
        return exact[0]
    if contains:
        return contains[0]
    raise RuntimeError(
        "Could not locate CTHRC1 symbol or Ensembl ID ENSG00000164932 in annotation columns"
    )

