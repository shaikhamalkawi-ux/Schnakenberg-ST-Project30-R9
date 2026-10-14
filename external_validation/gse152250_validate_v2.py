#!/usr/bin/env python3
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

CTHRC1_ENSEMBL = "ENSG00000164932"
TIMEPOINTS = (6, 24, 48)
PSEUDOCOUNT = 0.5


def parse_soft(path: Path) -> pd.DataFrame:
    rows = []
    current = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    rows.append(current)
                current = {"gsm": line.split("=", 1)[1].strip(), "descriptions": []}
            elif current is not None and line.startswith("!Sample_title = "):
                current["title"] = line.split("=", 1)[1].strip()
            elif current is not None and line.startswith("!Sample_description = "):
                current["descriptions"].append(line.split("=", 1)[1].strip())
    if current is not None:
        rows.append(current)

    cardiac = []
    pat = re.compile(r"^cardiac\s+(6|24|48)h\s+(CTRL|TGFb1)\s*(\d+)?$", re.I)
    for r in rows:
        title = r.get("title", "")
        m = pat.match(title.strip())
        if not m:
            continue
        cardiac.append({
            "gsm": r["gsm"],
            "title": title,
            "time_h": int(m.group(1)),
            "group": "TGFb1" if m.group(2).lower() == "tgfb1" else "CTRL",
            "replicate": int(m.group(3)) if m.group(3) else None,
            "description": " | ".join(r.get("descriptions", [])),
        })

    df = pd.DataFrame(cardiac)
    if len(df) != 18:
        raise RuntimeError(f"Expected 18 cardiac samples, found {len(df)}")
    return df.sort_values(["time_h", "group", "replicate", "gsm"]).reset_index(drop=True)


def norm(x: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(x).lower())


def load_expr(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)


def find_cthrc1_row(expr: pd.DataFrame) -> tuple[str, int]:
    for col in expr.columns:
        s = expr[col].astype(str).str.strip().str.upper()
        base = s.str.replace(r"\.\d+$", "", regex=True)
        hit = s.eq("CTHRC1") | base.eq(CTHRC1_ENSEMBL)
        if hit.any():
            return col, int(np.flatnonzero(hit.to_numpy())[0])
    raise RuntimeError("CTHRC1/ENSG00000164932 not found in any annotation column")


def expression_columns(expr: pd.DataFrame, annotation_col: str) -> list[str]:
    cols = []
    for col in expr.columns:
        if col == annotation_col:
            continue
        numeric = pd.to_numeric(expr[col], errors="coerce")
        if numeric.notna().mean() >= 0.95:
            cols.append(col)
    return cols


def sample_identifiers(r: pd.Series) -> list[str]:
    ids = [r["gsm"], r["title"], r["description"]]
    ids.extend(re.findall(r"[A-Za-z]+-\d+", str(r["description"])))
    return [str(x) for x in ids if str(x) and str(x) != "nan"]


def map_columns(meta: pd.DataFrame, expr_cols: list[str]) -> pd.DataFrame:
    norm_cols = {c: norm(c) for c in expr_cols}
    used = set()
    out = []
    for _, r in meta.iterrows():
        scored = []
        for col, ncol in norm_cols.items():
            if col in used:
                continue
            best = 0
            for ident in sample_identifiers(r):
                ni = norm(ident)
                if not ni:
                    continue
                if ncol == ni:
                    best = max(best, 100)
                elif len(ni) >= 6 and (ni in ncol or ncol in ni):
                    best = max(best, 80)
            if best:
                scored.append((best, col))
        scored.sort(reverse=True)
        if not scored:
            out.append({**r.to_dict(), "expr_column": None, "mapping_score": 0})
        else:
            score, col = scored[0]
            used.add(col)
            out.append({**r.to_dict(), "expr_column": col, "mapping_score": score})
    mapped = pd.DataFrame(out)
    if mapped["expr_column"].isna().any():
        bad = mapped.loc[mapped["expr_column"].isna(), ["gsm", "title", "description"]]
        raise RuntimeError("Unmapped cardiac samples: " + bad.to_json(orient="records"))
    return mapped


def exact_perm(values: np.ndarray, labels: np.ndarray) -> tuple[float, float]:
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
    p_one = float(np.mean(stats >= obs - 1e-12))
    p_two = float(np.mean(np.abs(stats) >= abs(obs) - 1e-12))
    return p_one, p_two


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fpkm", required=True)
    ap.add_argument("--soft", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    meta = parse_soft(Path(args.soft))
    expr = load_expr(Path(args.fpkm))

    schema = {
        "shape": [int(expr.shape[0]), int(expr.shape[1])],
        "columns": [str(c) for c in expr.columns],
    }
    (outdir / "GSE152250_expression_schema.json").write_text(json.dumps(schema, indent=2))

    annotation_col, row_pos = find_cthrc1_row(expr)
    expr_cols = expression_columns(expr, annotation_col)
    mapped = map_columns(meta, expr_cols)
    mapped.to_csv(outdir / "GSE152250_cardiac_sample_mapping.csv", index=False)

    row = expr.iloc[row_pos]
    sample_rows = []
    for _, m in mapped.iterrows():
        val = float(pd.to_numeric(row[m["expr_column"]], errors="raise"))
        sample_rows.append({
            "gsm": m["gsm"],
            "title": m["title"],
            "time_h": int(m["time_h"]),
            "group": m["group"],
            "replicate": int(m["replicate"]) if pd.notna(m["replicate"]) else None,
            "expr_column": m["expr_column"],
            "CTHRC1_FPKM": val,
            "log2_FPKM_plus_0.5": float(np.log2(val + PSEUDOCOUNT)),
        })
    samples = pd.DataFrame(sample_rows)
    samples.to_csv(outdir / "GSE152250_CTHRC1_sample_values.csv", index=False)

    summary_rows = []
    for t in TIMEPOINTS:
        d = samples[samples["time_h"] == t].copy()
        ctrl = d[d["group"] == "CTRL"]
        tgf = d[d["group"] == "TGFb1"]
        if len(ctrl) != 3 or len(tgf) != 3:
            raise RuntimeError(f"{t} h expected 3 CTRL and 3 TGFb1, got {len(ctrl)} and {len(tgf)}")

        mc = float(ctrl["CTHRC1_FPKM"].mean())
        mt = float(tgf["CTHRC1_FPKM"].mean())
        log2_ratio = math.log2((mt + 1e-12) / (mc + 1e-12))
        values = d["log2_FPKM_plus_0.5"].to_numpy()
        labels = (d["group"] == "TGFb1").astype(int).to_numpy()
        p1, p2 = exact_perm(values, labels)

        summary_rows.append({
            "time_h": t,
            "CTRL_mean_FPKM": mc,
            "TGFb1_mean_FPKM": mt,
            "CTRL_median_FPKM": float(ctrl["CTHRC1_FPKM"].median()),
            "TGFb1_median_FPKM": float(tgf["CTHRC1_FPKM"].median()),
            "log2_ratio_of_mean_FPKM": float(log2_ratio),
            "mean_log2_FPKM_plus_0.5_difference": float(
                tgf["log2_FPKM_plus_0.5"].mean() - ctrl["log2_FPKM_plus_0.5"].mean()
            ),
            "all_TGF_above_all_CTRL": bool(tgf["CTHRC1_FPKM"].min() > ctrl["CTHRC1_FPKM"].max()),
            "exact_perm_p_one_sided_directional": p1,
            "exact_perm_p_two_sided": p2,
            "direction": "UP" if mt > mc else ("DOWN" if mt < mc else "TIE"),
        })

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(outdir / "GSE152250_CTHRC1_timepoint_summary.csv", index=False)

    primary = summary.loc[summary["time_h"] == 24].iloc[0]
    passed = bool(primary["direction"] == "UP")
    result = {
        "dataset": "GSE152250",
        "source_file": "GSE152250_Fibroblasts_FPKM.txt.gz",
        "annotation_column": annotation_col,
        "CTHRC1_identifier": CTHRC1_ENSEMBL,
        "primary_time_h": 24,
        "prespecified_primary_criterion": "mean CTHRC1 FPKM in TGFb1 > CTRL at 24 h",
        "primary_criterion_pass": passed,
        "decision": "GO" if passed else "NO-GO",
        "donor_pairing_assumed": False,
        "cross_study_numeric_pooling": False,
        "timepoints": summary.to_dict(orient="records"),
    }
    (outdir / "GSE152250_external_validation_result.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )

    lines = [
        "# GSE152250 external validation",
        "",
        "Primary validation time point was prespecified as 24 h to match GSE97358.",
        "GSE152250 culture replicates are analyzed as independent/unpaired; no donor pairing is assumed.",
        "Cross-study effect sizes are not pooled.",
        "",
        f"Decision: **{result['decision']}**",
        "",
        "| Time | CTRL mean FPKM | TGFb1 mean FPKM | log2 ratio of means | one-sided exact p | two-sided exact p | all TGF > all CTRL |",
        "|---:|---:|---:|---:|---:|---:|:---:|",
    ]
    for _, r in summary.iterrows():
        lines.append(
            f"| {int(r['time_h'])} h | {r['CTRL_mean_FPKM']:.6g} | {r['TGFb1_mean_FPKM']:.6g} | "
            f"{r['log2_ratio_of_mean_FPKM']:.4f} | {r['exact_perm_p_one_sided_directional']:.4f} | "
            f"{r['exact_perm_p_two_sided']:.4f} | {bool(r['all_TGF_above_all_CTRL'])} |"
        )
    lines += [
        "",
        "Boundary: this tests independent within-study directional replication of the CTHRC1 response. "
        "It does not convert three culture replicates into donor-level evidence and does not pool FPKM-scale "
        "magnitudes with the donor-paired GSE97358 estimand.",
    ]
    (outdir / "GSE152250_external_validation_report.md").write_text("\n".join(lines) + "\n")

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
