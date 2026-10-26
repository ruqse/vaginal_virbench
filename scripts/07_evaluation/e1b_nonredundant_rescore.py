#!/usr/bin/env python3
"""Rescore the pooled 13-sample multi-evidence benchmark under the
non-redundant E1b coverage estimator, from existing per-contig tool calls.

No tool is re-run. Benchmark eligibility and labels depend only on the tier
(Tier 1-2 viral, Tier 0 negative, Tier 3 and -1 excluded), so the Table S26
tier changes only change which contigs are scored and with which label.

Labels: primary benchmark-ready contigs (>= 1,500 bp) with the Table S26 tier
changes applied. Calls: primary pooled predictions. Contigs that only enter the
benchmark under the corrected labels (Tier 3 -> 0) take their calls from the
no-E1b pooled predictions, which match the primary file on every shared row.

Before rescoring, the script reproduces the published primary confusion counts
and MCC (overall_metrics.tsv) and the published paired contig-bootstrap
intervals (paired_mcc_intervals.json, seed 42, 2,000 draws).

Usage:
    python3 scripts/07_evaluation/e1b_nonredundant_rescore.py
"""
import csv
import gzip
import json
import shutil
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/07_evaluation"))
from evaluate_metagenome import _auprc  # noqa: E402

TR = ROOT / "results/test_real"
OUT = TR / "secondary_benchmark_pooled13_e1b_nonredundant"
TABLES = ROOT / "results/tables"
EXPORT = {"overall_metrics.tsv": "table_s26_rescored_metrics.tsv",
          "pairwise_mcc.tsv": "table_s26_rescored_pairwise.tsv"}
MIN_LEN = 1500  # benchmark length gate (pooled_multievidence_benchmark.py --min-length)
PAIR_TOOLS = ["geNomad", "ViraLM", "VirSorter2", "VIBRANT", "Jaeger"]  # as revision_pairwise_bootstrap.py
LABEL = {"0": "negative", "1": "viral", "2": "viral"}  # Tier 3 and -1 are excluded


def load_pred(d):
    with gzip.open(TR / d / "pooled_tool_predictions.tsv.gz", "rt") as fh:
        return {(r["sample"], r["orig_contig_id"], r["tool"]): r
                for r in csv.DictReader(fh, delimiter="\t")}


def counts(y, p):
    return (int(((y == 1) & (p == 1)).sum()), int(((y == 0) & (p == 1)).sum()),
            int(((y == 0) & (p == 0)).sum()), int(((y == 1) & (p == 0)).sum()))


def mcc_bench(tp, fp, tn, fn):  # guard as pooled_multievidence_benchmark._mcc
    d = np.sqrt(max(tp + fp, 1) * max(tp + fn, 1) * max(tn + fp, 1) * max(tn + fn, 1))
    return (tp * tn - fp * fn) / d if d > 0 else 0.0


def mcc_boot(y, p):  # guard as revision_pairwise_bootstrap.mcc
    tp, fp, tn, fn = counts(y, p)
    d = np.sqrt(float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
    return float((tp * tn - fp * fn) / d) if d else 0.0


def metrics(rows):
    out = []
    for tool in sorted({r["tool"] for r in rows}):
        rs = [r for r in rows if r["tool"] == tool]
        y = np.array([r["label"] == "viral" for r in rs], int)
        p = np.array([int(r["y_pred"]) for r in rs])
        s = np.array([float(r["score"]) for r in rs])
        tp, fp, tn, fn = counts(y, p)
        prec = tp / (tp + fp) if tp + fp else 0.0
        rec = tp / (tp + fn) if tp + fn else 0.0
        out.append({"tool": tool, "TP": tp, "FP": fp, "TN": tn, "FN": fn,
                    "precision": round(prec, 4), "recall": round(rec, 4),
                    "F1": round(2 * prec * rec / (prec + rec), 4) if prec + rec else 0.0,
                    "MCC": round(mcc_bench(tp, fp, tn, fn), 4),
                    "AUPRC": round(_auprc(y, s), 4) if len(np.unique(s)) > 1 else "NA",
                    "n_total": len(rs)})
    out.sort(key=lambda r: r["MCC"], reverse=True)
    for i, r in enumerate(out, 1):
        r["rank"] = i
    return out


def paired(rows):
    """Contig bootstrap, 2,000 draws, seed 42, as revision_pairwise_bootstrap.py."""
    data = {}
    for r in rows:
        if r["tool"] in PAIR_TOOLS:
            d = data.setdefault(f'{r["sample"]}|{r["orig_contig_id"]}', {"y": int(r["label"] == "viral")})
            d[r["tool"]] = int(r["y_pred"])
    ids = sorted(data)
    y = np.array([data[c]["y"] for c in ids])
    P = {t: np.array([data[c][t] for c in ids]) for t in PAIR_TOOLS}
    rng = np.random.default_rng(42)
    pairs = list(combinations(PAIR_TOOLS, 2))
    point = {t: mcc_boot(y, P[t]) for t in PAIR_TOOLS}
    diffs = {q: [] for q in pairs}
    for _ in range(2000):
        draw = rng.integers(0, len(y), len(y))
        v = {t: mcc_boot(y[draw], P[t][draw]) for t in PAIR_TOOLS}
        for a, b in pairs:
            diffs[a, b].append(v[a] - v[b])
    res = []
    for a, b in pairs:
        lo, hi = np.percentile(diffs[a, b], [2.5, 97.5])
        res.append({"tool_a": a, "tool_b": b, "delta_a_minus_b": point[a] - point[b],
                    "ci_low": float(lo), "ci_high": float(hi)})
    return res


def main():
    prim = load_pred("secondary_benchmark_pooled13")
    noe1b = load_pred("secondary_benchmark_pooled13_no_e1b")
    for k in prim.keys() & noe1b.keys():
        assert (prim[k]["y_pred"], prim[k]["score"]) == (noe1b[k]["y_pred"], noe1b[k]["score"]), k
    prim_rows = list(prim.values())
    tools = {r["tool"] for r in prim_rows}

    # 1. Reproduce the published primary benchmark.
    with open(TR / "secondary_benchmark_pooled13/overall_metrics.tsv") as fh:
        pub = {r["tool"]: r for r in csv.DictReader(fh, delimiter="\t")}
    old = {r["tool"]: r for r in metrics(prim_rows)}
    for t, r in old.items():
        q = pub[t]
        assert (r["TP"], r["FP"], r["TN"], r["FN"]) == tuple(int(q[k]) for k in ("TP", "FP", "TN", "FN")), t
        assert abs(r["MCC"] - float(q["MCC"])) < 1e-4, t
        # Stored scores have 6 significant figures, so AUPRC may differ in the 4th decimal.
        assert q["AUPRC"] == "NA" or abs(r["AUPRC"] - float(q["AUPRC"])) <= 1.5e-4, t
    pub_pairs = json.loads((ROOT / "results/statistics/paired_mcc_intervals.json").read_text())["multi_evidence"]
    for a, b in zip(paired(prim_rows), pub_pairs):
        assert (a["tool_a"], a["tool_b"]) == (b["tool_a"], b["tool_b"])
        assert abs(a["ci_low"] - b["ci_low"]) < 1e-12 and abs(a["ci_high"] - b["ci_high"]) < 1e-12, a

    # 2. Apply the Table S26 tier changes at the benchmark length gate.
    drop, add = set(), set()
    with open(TABLES / "table_s26_e1b_parser_audit_tier_changes.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if int(r["contig_id"].split("_length_")[1].split("_")[0]) < MIN_LEN:
                continue
            k = (r["sample"], r["contig_id"])
            was, now = LABEL.get(r["tier_frozen"]), LABEL.get(r["tier_corrected"])
            if was and not now:
                drop.add(k)
            elif now and not was:
                add.add(k)
            elif was != now:
                raise SystemExit(f"viral/negative flip not handled: {k}")
    corr = [r for k, r in prim.items() if k[:2] not in drop]
    for k, r in noe1b.items():
        if k[:2] in add:
            assert r["label"] == "negative", k
            corr.append(r)
    assert len(corr) == len(prim_rows) + len(tools) * (len(add) - len(drop)), "missing calls for added contigs"

    # 3. Rescore and write.
    OUT.mkdir(parents=True, exist_ok=True)
    cols = ["tool", "rank", "MCC", "TP", "FP", "TN", "FN", "precision", "recall", "F1", "AUPRC", "n_total"]
    with open(OUT / "overall_metrics.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols + ["MCC_primary", "rank_primary", "delta_MCC", "AUPRC_primary"])
        for r in metrics(corr):
            o = old[r["tool"]]
            w.writerow([r[c] for c in cols] + [o["MCC"], o["rank"], round(r["MCC"] - o["MCC"], 4),
                                               pub[r["tool"]]["AUPRC"]])
    with open(OUT / "pairwise_mcc.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["context", "tool_a", "tool_b", "delta_a_minus_b", "ci_low", "ci_high",
                    "ci_excludes_zero", "bootstrap_draws", "unit", "seed"])
        for q in paired(corr):
            w.writerow(["multi_evidence_e1b_nonredundant", q["tool_a"], q["tool_b"], q["delta_a_minus_b"],
                        q["ci_low"], q["ci_high"], q["ci_low"] > 0 or q["ci_high"] < 0, 2000, "contig", 42])
    for f, dest in EXPORT.items():
        shutil.copyfile(OUT / f, TABLES / dest)
    print(f"positives dropped {len(drop)}, negatives added {len(add)}, contigs scored {len(corr) // len(tools)}")
    print((OUT / "overall_metrics.tsv").read_text())


if __name__ == "__main__":
    main()
