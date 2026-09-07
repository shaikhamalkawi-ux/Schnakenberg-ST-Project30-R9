#!/usr/bin/env python3
"""Deterministic target-row extractor for the Colak et al. 2016 PLOS S1 XLS.

This helper deliberately refuses to guess column meanings. Supply the original
publisher XLS plus the exact gene and fold-change column names. Optional
significance columns can also be supplied. It writes all matching source rows
for the predefined target panel without imputing missing targets.

Example:
  python extract_saudi_targets.py pone.0162669.s005.xls \
      --gene-column GeneSymbol --fold-column FoldChange \
      --fdr-column FDR --output saudi_target_rows.csv

The script has NOT been executed in this study because the original XLS
binary was not available in the active analysis materials.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

TARGETS = ["CTHRC1","TGFB1","TGFBR1","TGFBR2","POSTN","FN1","LOX","COL1A1","COL3A1","ASPN"]

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("xls", type=Path)
    ap.add_argument("--sheet", default=0, help="Sheet name or zero-based sheet index")
    ap.add_argument("--gene-column", required=True)
    ap.add_argument("--fold-column", required=True)
    ap.add_argument("--p-column")
    ap.add_argument("--fdr-column")
    ap.add_argument("--output", type=Path, default=Path("saudi_target_rows.csv"))
    args = ap.parse_args()

    if not args.xls.exists():
        raise SystemExit(f"Input file not found: {args.xls}")
    sheet = int(args.sheet) if str(args.sheet).isdigit() else args.sheet
    df = pd.read_excel(args.xls, sheet_name=sheet)
    required = [args.gene_column, args.fold_column]
    optional = [c for c in [args.p_column, args.fdr_column] if c]
    missing = [c for c in required + optional if c not in df.columns]
    if missing:
        raise SystemExit(
            "Requested columns are not present: " + ", ".join(missing) +
            "\nAvailable columns: " + ", ".join(map(str, df.columns))
        )

    genes = df[args.gene_column].astype(str).str.strip().str.upper()
    hit = df.loc[genes.isin(TARGETS), required + optional].copy()
    hit.insert(0, "target", genes.loc[hit.index])
    hit.insert(1, "source_row_index", hit.index.astype(int))
    hit = hit.sort_values(["target", "source_row_index"])
    hit.to_csv(args.output, index=False)

    observed = set(hit["target"])
    status = pd.DataFrame({
        "target": TARGETS,
        "matched_source_row": [t in observed for t in TARGETS],
        "interpretation": ["source row extracted" if t in observed else "no matching row extracted; do not interpret as biological absence" for t in TARGETS],
    })
    status.to_csv(args.output.with_name(args.output.stem + "_status.csv"), index=False)
    print(f"Wrote {len(hit)} matching source rows to {args.output}")
    print(f"Wrote target presence status to {args.output.with_name(args.output.stem + '_status.csv')}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
