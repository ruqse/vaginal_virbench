#!/usr/bin/env python3
"""Sensitivity of the MiTCH rank-replication correlation (Table S22, Panel D).

The headline replication statistic is a Spearman correlation between per-tool
pooled MCC in MiTCH and in the matched four-evidence (no-E5) UChoose benchmark,
over all 14 tools (rho = 0.943). Six of those tools sit at MCC <= 0.07 in MiTCH,
so this script asks how much of the agreement is carried by that weak cluster,
and separates rank agreement from agreement in level:

  * rho recomputed over three tool sets (all 14; MiTCH MCC >= 0.10; the five
    decision-relevant tools of Table 1),
  * a sample-block bootstrap CI for each (resample the 30 MiTCH samples,
    recompute pooled MCC, re-rank, correlate against the fixed UChoose ranking),
  * a permutation p-value (shuffle the UChoose values against the MiTCH ones),
  * Lin's concordance correlation coefficient, which unlike rho penalises the
    systematic level shift between cohorts.

Ranks use average ranks for ties. The production path
(scripts/10_vmgc_crosscheck/scaleup_analysis.py) uses argsort().argsort(),
i.e. ordinal ranks; the two agree to three decimals on these data because no
consequential ties occur, but average ranks are the correct general choice.

Writes the Panel D rows of results/tables/table_s22_cross_track_spearman.tsv.
"""
import csv
import gzip
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PRED = ROOT / "results/test_real/viralm_investigation/mitch_scaleup30/pooled_tool_predictions.tsv.gz"
S27 = ROOT / "results/tables/table_s27_mitch_ranking_ci.tsv"
OUT = ROOT / "results/tables/table_s22d_replication_sensitivity.tsv"

B_BOOT = 2000
N_PERM = 20000
SEED = 12345
DECISION_TOOLS = ["geNomad", "ViraLM", "VirSorter2", "VIBRANT", "Jaeger"]
WEAK_CUTOFF = 0.10


def mcc(tp, fp, tn, fn):
    d = math.sqrt(max(tp + fp, 1) * max(tp + fn, 1) * max(tn + fp, 1) * max(tn + fn, 1))
    return (tp * tn - fp * fn) / d if d > 0 else 0.0


def load_confusion(path):
    conf, tools, samples = {}, set(), set()
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tool, sample = r["tool"], r["sample"]
            tools.add(tool)
            samples.add(sample)
            v = conf.setdefault((tool, sample), [0, 0, 0, 0])
            viral = r["label"] == "viral"
            pred = r["y_pred"] == "1"
            if viral:
                v[0] += pred
                v[3] += not pred
            else:
                v[1] += pred
                v[2] += not pred
    return conf, sorted(tools), sorted(samples)


def pooled_mcc(conf, tool, samples):
    tp = fp = tn = fn = 0
    for s in samples:
        a = conf.get((tool, s))
        if a:
            tp += a[0]
            fp += a[1]
            tn += a[2]
            fn += a[3]
    return mcc(tp, fp, tn, fn)


def average_ranks(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return ranks


def spearman(x, y):
    return float(np.corrcoef(average_ranks(list(x)), average_ranks(list(y)))[0, 1])


def lin_ccc(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    cov = np.cov(x, y, bias=True)[0, 1]
    return float(2 * cov / (x.var() + y.var() + (x.mean() - y.mean()) ** 2))


def main():
    conf, tools, samples = load_confusion(PRED)
    point = {t: pooled_mcc(conf, t, samples) for t in tools}

    published = {}
    with open(S27) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            published[r["tool"]] = (float(r["MCC"]), float(r["disc_no_e5_MCC"]))
    drift = max(abs(point[t] - published[t][0]) for t in tools)
    if drift > 5e-4:
        raise SystemExit(f"pooled MCC disagrees with Table S27 (max diff {drift:.5f})")
    uchoose = {t: published[t][1] for t in tools}

    tool_sets = [
        ("all 14 tools", tools),
        (f"MiTCH MCC >= {WEAK_CUTOFF:.2f}", [t for t in tools if point[t] >= WEAK_CUTOFF]),
        ("five decision-relevant tools", [t for t in DECISION_TOOLS if t in point]),
    ]

    rng = np.random.default_rng(SEED)
    rows = [[
        "tool_set", "n_tools", "spearman_rho", "boot_CI_low", "boot_CI_high",
        "permutation_p", "lin_ccc", "mean_MCC_MiTCH", "mean_MCC_UChoose_noE5",
    ]]
    for name, subset in tool_sets:
        x = [point[t] for t in subset]
        y = [uchoose[t] for t in subset]
        rho = spearman(x, y)

        boots = []
        for _ in range(B_BOOT):
            draw = [samples[i] for i in rng.integers(0, len(samples), len(samples))]
            m = {t: pooled_mcc(conf, t, draw) for t in subset}
            boots.append(spearman([m[t] for t in subset], y))
        lo, hi = np.percentile(boots, [2.5, 97.5])

        perm = [spearman(x, list(rng.permutation(y))) for _ in range(N_PERM)]
        p = (np.sum(np.abs(perm) >= abs(rho) - 1e-12) + 1) / (N_PERM + 1)

        rows.append([
            name, len(subset), round(rho, 3), round(float(lo), 3), round(float(hi), 3),
            # p is (hits + 1) / (N + 1), so the floor is 1/(N+1); report it as a bound
            # rather than printing a rounded 0.0000 that reads as an exact zero.
            f"{p:.4f}" if p >= 1e-4 else f"<{1.0 / (N_PERM + 1):.5f}",
            round(lin_ccc(x, y), 3),
            round(float(np.mean(x)), 3), round(float(np.mean(y)), 3),
        ])

    with open(OUT, "w", newline="") as fh:
        csv.writer(fh, delimiter="\t").writerows(rows)
    for r in rows:
        print("\t".join(str(c) for c in r))
    print(f"\nwrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
