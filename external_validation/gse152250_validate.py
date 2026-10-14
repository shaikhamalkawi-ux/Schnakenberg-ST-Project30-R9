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


TARGET_GENES = ["CTHRC1", "POSTN", "FN1", "LOX", "COL1A1", "COL3A1", "ASPN"]
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
    exact = []
    contains = []
    for c in df.columns:
        vals = df[c].astype(str)
        if vals.str.upper().eq("CTHRC1").any():
            exact.append(c)
        elif vals.str.upper().str.contains(r"\bCTHRC1\b", regex=True).any():
            contains.append(c)
    if exact:
        return exact[0]
    if contains:
        return contains[0]
    raise RuntimeError("Could not locate CTHRC1 in any annotation column")


def candidate_ids(row: pd.Series) -> list[str]:
    ids = [row["gsm"], row["title"], row["description"]]
    desc = str(row["description"])
    ids += re.findall(r"[A-Za-z]+-\d+", desc)
    return [x for x in ids if x and str(x) != "nan"]


def map_sample_columns(meta: pd.DataFrame, expr: pd.DataFrame, symbol_col: str) -> pd.DataFrame:
    # Expression columns are inferred by numeric convertibility; annotations are excluded.
    numeric_cols = []
    for c in expr.columns:
        if c == symbol_col:
            continue
        s = pd.to_numeric(expr[c], errors="coerce")
        if s.notna().mean() > 0.95:
            numeric_cols.append(c)

    if len(numeric_cols) < 18:
        raise RuntimeError(f"Found only {len(numeric_cols)} numeric expression columns")

    col_norm = {c: norm_token(c) for c in numeric_cols}
    mapped = []
    used = set()

    for _, r in meta.iterrows():
        ids = candidate_ids(r)
        choices = []
        for c, nc in col_norm.items():
            if c in used:
                continue
            score = 0
            for ident in ids:
                ni = norm_token(ident)
                if not ni:
                    continue
                if nc == ni:
                    score = max(score, 100)
                elif ni in nc or nc in ni:
                    score = max(score, 80 if len(ni) >= 6 else 20)
            if score:
                choices.append((score, c))
        choices.sort(reverse=True)
        if choices:
            c = choices[0][1]
            used.add(c)
            mapped.append({**r.to_dict(), "expr_column": c, "mapping_score": choices[0][0]})
        else:
            mapped.append({**r.to_dict(), "expr_column": None, "mapping_score": 0})

    out = pd.DataFrame(mapped)
    missing = out["expr_column"].isna().sum()
    if missing:
        # Last-resort: if there are exactly 18 cardiac-expression columns and all metadata
        # are unmapped, use table order only when column names themselves encode time/group.
        raise RuntimeError(
            f"Could not map {missing} cardiac samples to FPKM columns. "
            f"Expression columns: {numeric_cols[:30]}"
        )
    if out["expr_column"].duplicated().any():
        raise RuntimeError("Duplicate expression-column mapping detected")
    return out


def extract_gene_rows(expr: pd.DataFrame, symbol_col: str, genes: list[str]) -> pd.DataFrame:
    symbols = expr[symbol_col].astype(str).str.upper().str.strip()
    sub = expr.loc[symbols.isin([g.upper() for g in genes])].copy()
    if sub.empty:
        raise RuntimeError("No target-gene rows found")
    sub["_SYMBOL_"] = symbols.loc[sub.index]
    # Collapse duplicate symbol rows conservatively by summing FPKM columns.
    return sub


def exact_permutation(values: np.ndarray, labels: np.ndarray) -> tuple[float, float]:
    """Exact one-sided and two-sided permutation p-values for difference in means."""
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels, dtype=int)
    n_t = int(labels.sum())
    obs = values[labels == 1].mean() - values[labels == 0].mean()
    stats = []
    for idx in itertools.combinations(range(len(values)), n_t):
        lab = np.zeros(len(values), dtype=int)
        lab[list(idx)] = 1
        stats.append(values[lab == 1].mean() - values[lab == 0].mean())
    stats = np.asarray(stats)
    tol = 1e-12
    p_one = float(np.mean(stats >= obs - tol))
    p_two = float(np.mean(np.abs(stats) >= abs(obs) - tol))
    return p_one, p_two


def summarize_gene(gene: str, row: pd.Series, mapping: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_rows = []
    for _, m in mapping.iterrows():
        val = pd.to_numeric(row[m["expr_column"]], errors="coerce")
        sample_rows.append(
            {
                "gene": gene,
                "gsm": m["gsm"],
                "title": m["title"],
                "time_h": int(m["time_h"]),
                "group": m["group"],
                "replicate": m["replicate"],
                "expr_column": m["expr_column"],
                "FPKM": float(val),
                "log2_FPKM_plus_0.5": float(np.log2(float(val) + PSEUDOCOUNT)),
            }
        )
    sample_df = pd.DataFrame(sample_rows)

    summaries = []
    for t in TIMEPOINTS:
        d = sample_df[sample_df.time_h == t].copy()
        ctrl = d[d.group == "CTRL"]
        tgf = d[d.group == "TGFb1"]
        if len(ctrl) != 3 or len(tgf) != 3:
            raise RuntimeError(f"{gene} {t} h does not have 3 CTRL + 3 TGF samples")

        mean_ctrl = ctrl.FPKM.mean()
        mean_tgf = tgf.FPKM.mean()
        if mean_ctrl > 0 and mean_tgf > 0:
            log2_ratio_means = math.log2(mean_tgf / mean_ctrl)
        else:
            log2_ratio_means = math.log2((mean_tgf + 1e-9) / (mean_ctrl + 1e-9))

        vals = d["log2_FPKM_plus_0.5"].to_numpy()
        labs = (d["group"] == "TGFb1").astype(int).to_numpy()
        p_one, p_two = exact_permutation(vals, labs)

        summaries.append(
            {
                "gene": gene,
                "time_h": t,
                "n_CTRL": len(ctrl),
                "n_TGFb1": len(tgf),
                "CTRL_mean_FPKM": float(mean_ctrl),
                "TGFb1_mean_FPKM": float(mean_tgf),
                "CTRL_median_FPKM": float(ctrl.FPKM.median()),
                "TGFb1_median_FPKM": float(tgf.FPKM.median()),
                "log2_ratio_of_mean_FPKM": float(log2_ratio_means),
                "mean_log2_FPKM_plus_0.5_difference": float(
                    tgf.log2_FPKM_plus_0.5.mean() - ctrl.log2_FPKM_plus_0.5.mean()
                ),
                "all_TGF_above_all_CTRL": bool(tgf.FPKM.min() > ctrl.FPKM.max()),
                "exact_perm_p_one_sided_directional": p_one,
                "exact_perm_p_two_sided": p_two,
                "direction": "UP" if mean_tgf > mean_ctrl else ("DOWN" if mean_tgf < mean_ctrl else "TIE"),
            }
        )
    return sample_df, pd.DataFrame(summaries)


def choose_gene_row(sub: pd.DataFrame, symbol_col: str, gene: str, sample_cols: list[str]) -> pd.Series:
    rows = sub[sub["_SYMBOL_"] == gene.upper()].copy()
    if len(rows) == 1:
        return rows.iloc[0]
    # If duplicate annotations exist, choose the row with largest total expression,
    # and record duplicates in the schema report.
    totals = rows[sample_cols].apply(pd.to_numeric, errors="coerce").sum(axis=1)
    return rows.loc[totals.idxmax()]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fpkm", required=True)
    ap.add_argument("--soft", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    fpkm = Path(args.fpkm)
    soft = Path(args.soft)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    meta = parse_soft(soft)
    expr = load_expression(fpkm)
    symbol_col = find_symbol_column(expr)
    mapping = map_sample_columns(meta, expr, symbol_col)
    mapping.to_csv(outdir / "GSE152250_cardiac_sample_mapping.csv", index=False)

    sample_cols = mapping["expr_column"].tolist()
    sub = extract_gene_rows(expr, symbol_col, TARGET_GENES)

    all_samples = []
    all_summary = []
    duplicate_counts = {}
    for gene in TARGET_GENES:
        gene_rows = sub[sub["_SYMBOL_"] == gene.upper()]
        duplicate_counts[gene] = int(len(gene_rows))
        if len(gene_rows) == 0:
            continue
        row = choose_gene_row(sub, symbol_col, gene, sample_cols)
        sdf, sumdf = summarize_gene(gene, row, mapping)
        all_samples.append(sdf)
        all_summary.append(sumdf)

    samples = pd.concat(all_samples, ignore_index=True)
    summary = pd.concat(all_summary, ignore_index=True)
    samples.to_csv(outdir / "GSE152250_target_gene_sample_values.csv", index=False)
    summary.to_csv(outdir / "GSE152250_target_gene_timepoint_summary.csv", index=False)

    c = summary[summary.gene == "CTHRC1"].set_index("time_h")
    if 24 not in c.index:
        raise RuntimeError("CTHRC1 24-h result missing")
    primary = c.loc[24]
    go = bool(primary["direction"] == "UP")

    result = {
        "dataset": "GSE152250",
        "primary_validation_time_h": 24,
        "prespecified_primary_criterion": "CTHRC1 mean FPKM under TGFb1 > CTRL at 24 h",
        "primary_criterion_pass": go,
        "decision": "GO" if go else "NO-GO",
        "interpretation": (
            "Independent human cardiac-fibroblast directional external replication"
            if go
            else "Independent dataset does not replicate the prespecified 24-h CTHRC1 direction"
        ),
        "cross_study_pooling": False,
        "donor_pairing_assumed": False,
        "pseudocount_for_log_transform": PSEUDOCOUNT,
        "symbol_column": symbol_col,
        "expression_table_shape": list(expr.shape),
        "cardiac_sample_count": int(len(mapping)),
        "target_gene_duplicate_row_counts": duplicate_counts,
        "CTHRC1": {
            str(int(t)): {
                k: (bool(v) if isinstance(v, (np.bool_, bool)) else float(v) if isinstance(v, (np.floating, float)) else int(v) if isinstance(v, (np.integer, int)) else v)
                for k, v in row.items()
                if k not in {"gene", "time_h"}
            }
            for t, row in c.reset_index().set_index("time_h").iterrows()
        },
    }
    (outdir / "GSE152250_external_validation_result.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )

    lines = [
        "# GSE152250 external validation",
        "",
        "Prespecified primary target: CTHRC1 at 24 h, matching the GSE97358 primary time point.",
        "Independent culture replicates are treated as unpaired. Cross-study effect sizes are not pooled.",
        "",
        f"Decision: **{result['decision']}**",
        "",
        "| Time | CTRL mean FPKM | TGFb1 mean FPKM | log2 ratio of means | exact one-sided p | exact two-sided p | Direction |",
        "|---:|---:|---:|---:|---:|---:|:---|",
    ]
    for t in TIMEPOINTS:
        r = c.loc[t]
        lines.append(
            f"| {t} h | {r['CTRL_mean_FPKM']:.6g} | {r['TGFb1_mean_FPKM']:.6g} | "
            f"{r['log2_ratio_of_mean_FPKM']:.4f} | {r['exact_perm_p_one_sided_directional']:.4f} | "
            f"{r['exact_perm_p_two_sided']:.4f} | {r['direction']} |"
        )
    lines += [
        "",
        "Interpretation boundary: this is a within-study directional/time-course replication in an independent human cardiac-fibroblast dataset. It does not create donor-level inference, and its FPKM-scale effect is not pooled numerically with the GSE97358 donor-paired estimand.",
    ]
    (outdir / "GSE152250_external_validation_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
