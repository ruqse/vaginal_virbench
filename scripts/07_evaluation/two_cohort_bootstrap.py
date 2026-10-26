#!/usr/bin/env python3
"""Primary-cohort MCC CIs and a two-cohort bootstrap for the MiTCH replication.

Table S27 / Fig. 7A give MiTCH pooled MCC with a sample-block bootstrap CI but
the matched primary (UChoose, four-evidence, no-E5) MCC as a bare point value,
and the Table S22D Spearman CI resamples MiTCH only, holding the UChoose values
fixed. This script adds the missing half:

  * primary no-E5 pooled MCC with a block-bootstrap 95% CI. Blocks are
    participants (11), because UC055 and UC093 each contribute two visits
    (supplementary_methods.md, line 91); sample blocks (13) are reported as a
    sensitivity check;
  * per-tool delta MCC (MiTCH - primary) with a CI from resampling both
    cohorts independently;
  * Spearman rho between cohorts with three CIs: MiTCH resampled only (the
    published S22D scheme), primary resampled only, and both resampled.

MiTCH blocks are samples (30), as in Table S27. Percentile intervals, B = 2,000.

It also writes the Table S27 rates sheets (results/tables/table_s27_rates_*.tsv):
per-cohort TPR and FPR, MiTCH MCC re-scored at the primary viral fraction, the
rank correlation of TPR - FPR (independent of viral fraction), and the viral
fraction at which geNomad's MCC reaches VIBRANT's under each cohort's rates.
"""
import csv
import gzip
import math
import shutil
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[2]
MITCH_PRED = ROOT / "results/test_real/viralm_investigation/mitch_scaleup30/pooled_tool_predictions.tsv.gz"
PRIM_DIR = ROOT / "results/test_real/secondary_benchmark_pooled13_no_e5"
PRIM_PRED = PRIM_DIR / "pooled_tool_predictions.tsv.gz"
PRIM_METRICS = PRIM_DIR / "overall_metrics.tsv"
S27 = ROOT / "results/tables/table_s27_mitch_ranking_ci.tsv"
S22D = ROOT / "results/tables/table_s22d_replication_sensitivity.tsv"
OUT = ROOT / "results/test_real/two_cohort_bootstrap"

TABLES = ROOT / "results/tables"
EXPORT = {"cross_cohort_rates.tsv": "table_s27_rates_cross_cohort.tsv",
          "cross_cohort_summary.tsv": "table_s27_rates_summary.tsv"}
GRID = np.round(np.arange(0.001, 0.5001, 0.0005), 4)   # viral fractions searched for the crossover

B = 2000
SEED = 12345
DECISION_TOOLS = ["geNomad", "ViraLM", "VirSorter2", "VIBRANT", "Jaeger"]
WEAK_CUTOFF = 0.10


def load_counts(path, block_of):
    """Return tools, blocks, and counts[tool, block, (tp, fp, tn, fn)]."""
    acc = {}
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            v = acc.setdefault((r["tool"], block_of(r["sample"])), [0, 0, 0, 0])
            viral, pred = r["label"] == "viral", r["y_pred"] == "1"
            v[0 if viral and pred else 3 if viral else 1 if pred else 2] += 1
    tools = sorted({t for t, _ in acc})
    blocks = sorted({b for _, b in acc})
    arr = np.zeros((len(tools), len(blocks), 4))
    for (t, b), v in acc.items():
        arr[tools.index(t), blocks.index(b)] = v
    return tools, blocks, arr


def mcc(c):
    """MCC over the last axis (tp, fp, tn, fn); same max(., 1) guard as scaleup_analysis.py."""
    tp, fp, tn, fn = (c[..., i] for i in range(4))
    d = np.sqrt(np.maximum(tp + fp, 1) * np.maximum(tp + fn, 1)
                * np.maximum(tn + fp, 1) * np.maximum(tn + fn, 1))
    return (tp * tn - fp * fn) / d


def boot_counts(rng, counts):
    """B x tools x 4 pooled counts under block resampling with replacement."""
    n = counts.shape[1]
    w = np.stack([np.bincount(rng.integers(0, n, n), minlength=n) for _ in range(B)])
    return np.einsum("bk,tkc->btc", w, counts)


def rates(c):
    """TPR, FPR and labelled viral fraction over the last axis (tp, fp, tn, fn)."""
    tp, fp, tn, fn = (c[..., i] for i in range(4))
    with np.errstate(invalid="ignore", divide="ignore"):
        return tp / (tp + fn), fp / (fp + tn), (tp + fn) / c.sum(axis=-1)


def mcc_at(tpr, fpr, pi):
    """Expected MCC at viral fraction pi with TPR and FPR held fixed; 0 for a tool with no calls."""
    q = pi * tpr + (1 - pi) * fpr
    with np.errstate(invalid="ignore", divide="ignore"):
        v = (tpr - fpr) * np.sqrt(pi * (1 - pi) / (q * (1 - q)))
    return np.where(q > 0, v, 0.0)


def crossover(tpr_a, fpr_a, tpr_b, fpr_b):
    """Lowest viral fraction at which tool a's MCC reaches tool b's (NaN if none in GRID).

    Inputs are length-B vectors; returns a length-B vector."""
    d = mcc_at(tpr_a[:, None], fpr_a[:, None], GRID) - mcc_at(tpr_b[:, None], fpr_b[:, None], GRID)
    hit = d >= 0
    first = hit.argmax(axis=1)
    return np.where(hit.any(axis=1), GRID[first], np.nan)


def spearman_rows(x, y):
    """Row-wise Spearman (average ranks) between two B x k matrices."""
    rx, ry = rankdata(x, axis=1), rankdata(y, axis=1)
    rx -= rx.mean(axis=1, keepdims=True)
    ry -= ry.mean(axis=1, keepdims=True)
    den = np.sqrt((rx ** 2).sum(axis=1) * (ry ** 2).sum(axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return (rx * ry).sum(axis=1) / den


def ci(v):
    v = v[np.isfinite(v)]
    lo, hi = np.percentile(v, [2.5, 97.5])
    return float(lo), float(hi)


def r4(x):
    return round(float(x), 4)


def main():
    m_tools, m_blocks, m_counts = load_counts(MITCH_PRED, lambda s: s)
    p_tools, p_part, p_counts = load_counts(PRIM_PRED, lambda s: s.split("_")[0])
    _, p_samp, p_counts_s = load_counts(PRIM_PRED, lambda s: s)
    assert m_tools == p_tools, "tool sets differ between cohorts"
    tools = m_tools
    print(f"MiTCH: {len(m_blocks)} samples | primary: {len(p_samp)} samples, "
          f"{len(p_part)} participants | {len(tools)} tools")

    m_pt = mcc(m_counts.sum(axis=1))
    p_pt = mcc(p_counts.sum(axis=1))

    # Point estimates must reproduce the published values.
    with open(S27) as fh:
        s27 = {r["tool"]: r for r in csv.DictReader(fh, delimiter="\t")}
    with open(PRIM_METRICS) as fh:
        om = {r["tool"]: float(r["MCC"]) for r in csv.DictReader(fh, delimiter="\t")}
    for i, t in enumerate(tools):
        for got, want, what in ((m_pt[i], float(s27[t]["MCC"]), "MiTCH vs S27"),
                                (p_pt[i], om[t], "primary vs overall_metrics"),
                                (p_pt[i], float(s27[t]["disc_no_e5_MCC"]), "primary vs S27")):
            if abs(got - want) > 5e-4:
                raise SystemExit(f"{t}: {what} mismatch ({got:.4f} vs {want:.4f})")

    rng = np.random.default_rng(SEED)
    m_bc = boot_counts(rng, m_counts)    # MiTCH, sample blocks
    p_bc = boot_counts(rng, p_counts)    # primary, participant blocks
    p_b, m_b = mcc(p_bc), mcc(m_bc)
    p_bs = mcc(boot_counts(rng, p_counts_s))   # primary, sample blocks (sensitivity)

    # ---- per-tool table ----
    OUT.mkdir(parents=True, exist_ok=True)
    order = np.argsort(-m_pt)
    hdr = ["tool", "MCC_MiTCH", "MiTCH_CI_low", "MiTCH_CI_high",
           "MCC_primary_noE5", "primary_CI_low", "primary_CI_high",
           "primary_CI_low_sampleblocks", "primary_CI_high_sampleblocks",
           "delta_MiTCH_minus_primary", "delta_CI_low", "delta_CI_high", "delta_CI_excludes_0"]
    rows = [hdr]
    deltas, excl_neg = [], 0
    for i in order:
        t = tools[i]
        d_b = m_b[:, i] - p_b[:, i]
        dlo, dhi = ci(d_b)
        excl = dlo > 0 or dhi < 0
        deltas.append(m_pt[i] - p_pt[i])
        excl_neg += excl and dhi < 0
        rows.append([t, r4(m_pt[i]), *map(r4, ci(m_b[:, i])),
                     r4(p_pt[i]), *map(r4, ci(p_b[:, i])), *map(r4, ci(p_bs[:, i])),
                     r4(m_pt[i] - p_pt[i]), r4(dlo), r4(dhi), excl])
    with open(OUT / "primary_mcc_ci_and_delta.tsv", "w", newline="") as fh:
        csv.writer(fh, delimiter="\t").writerows(rows)

    # ---- rho table ----
    with open(S22D) as fh:
        s22d = {r["tool_set"]: r for r in csv.DictReader(fh, delimiter="\t")}
    sets = [("all 14 tools", tools),
            (f"MiTCH MCC >= {WEAK_CUTOFF:.2f}", [t for t in tools if m_pt[tools.index(t)] >= WEAK_CUTOFF]),
            ("five decision-relevant tools", [t for t in DECISION_TOOLS if t in tools])]
    rho_rows = [["tool_set", "n_tools", "spearman_rho",
                 "S22D_published_CI_low", "S22D_published_CI_high",
                 "MiTCH_only_CI_low", "MiTCH_only_CI_high",
                 "primary_only_CI_low", "primary_only_CI_high",
                 "two_cohort_CI_low", "two_cohort_CI_high",
                 "two_cohort_CI_low_primary_sampleblocks", "two_cohort_CI_high_primary_sampleblocks",
                 "n_degenerate_draws"]]
    for name, subset in sets:
        idx = [tools.index(t) for t in subset]
        mx, px = m_pt[idx][None, :], p_pt[idx][None, :]
        rho = float(spearman_rows(mx, px)[0])
        mo = spearman_rows(m_b[:, idx], np.repeat(px, B, 0))
        po = spearman_rows(np.repeat(mx, B, 0), p_b[:, idx])
        both = spearman_rows(m_b[:, idx], p_b[:, idx])
        both_s = spearman_rows(m_b[:, idx], p_bs[:, idx])
        nan = int(sum((~np.isfinite(a)).sum() for a in (mo, po, both, both_s)))
        pub = s22d.get(name, {})
        rho_rows.append([name, len(subset), round(rho, 3),
                         pub.get("boot_CI_low", "NA"), pub.get("boot_CI_high", "NA"),
                         *[round(x, 3) for a in (mo, po, both, both_s) for x in ci(a)], nan])
    with open(OUT / "two_cohort_rho.tsv", "w", newline="") as fh:
        csv.writer(fh, delimiter="\t").writerows(rho_rows)

    for table in (rows, rho_rows):
        for r in table:
            print("\t".join(str(c) for c in r))
        print()
    n_lower = sum(d < 0 for d in deltas)
    print(f"MiTCH lower than primary: {n_lower}/{len(deltas)} tools, "
          f"median delta {np.median(deltas):.3f}; "
          f"{excl_neg} with delta CI entirely below 0")
    cross_cohort_rates(tools, m_counts.sum(axis=1), p_counts.sum(axis=1), m_bc, p_bc, m_pt, p_pt, m_b, p_b)
    print(f"wrote {OUT.relative_to(ROOT)}/")


def cross_cohort_rates(tools, m_c, p_c, m_bc, p_bc, m_pt, p_pt, m_b, p_b):
    """Table S27 rates sheets: how much of the cross-cohort MCC difference is viral fraction.

    MiTCH TPR and FPR are re-scored at the primary viral fraction (each draw at
    that draw's primary fraction), and the geNomad-VIBRANT crossover is located
    under each cohort's error rates. Same draws as the MCC intervals above.
    """
    (mT, mF, mpi), (pT, pF, ppi) = rates(m_c), rates(p_c)
    (mTb, mFb, mpib), (pTb, pFb, ppib) = rates(m_bc), rates(p_bc)
    # self-check: re-scoring at a cohort's own fraction must return its MCC
    called = (m_c[:, 0] + m_c[:, 1]) > 0
    assert np.allclose(mcc_at(mT, mF, mpi)[called], m_pt[called], atol=1e-9), "mcc_at != MCC at native fraction"
    m_at = mcc_at(mT, mF, ppi[0])
    m_atb = mcc_at(mTb, mFb, ppib[:, :1])
    d_at, d_atb = m_at - p_pt, m_atb - p_b

    hdr = ["tool", "TPR_primary", "TPR_primary_CI_low", "TPR_primary_CI_high",
           "FPR_primary", "FPR_primary_CI_low", "FPR_primary_CI_high",
           "TPR_MiTCH", "TPR_MiTCH_CI_low", "TPR_MiTCH_CI_high",
           "FPR_MiTCH", "FPR_MiTCH_CI_low", "FPR_MiTCH_CI_high",
           "MCC_primary_noE5", "MCC_MiTCH", "MCC_MiTCH_at_primary_fraction",
           "MCC_MiTCH_at_primary_fraction_CI_low", "MCC_MiTCH_at_primary_fraction_CI_high",
           "delta_at_primary_fraction", "delta_at_primary_fraction_CI_low",
           "delta_at_primary_fraction_CI_high", "delta_at_primary_fraction_CI_excludes_0"]
    rows = [hdr]
    for i in np.argsort(-m_pt):
        lo, hi = ci(d_atb[:, i])
        rows.append([tools[i], r4(pT[i]), *map(r4, ci(pTb[:, i])), r4(pF[i]), *map(r4, ci(pFb[:, i])),
                     r4(mT[i]), *map(r4, ci(mTb[:, i])), r4(mF[i]), *map(r4, ci(mFb[:, i])),
                     r4(p_pt[i]), r4(m_pt[i]), r4(m_at[i]), *map(r4, ci(m_atb[:, i])),
                     r4(d_at[i]), r4(lo), r4(hi), lo > 0 or hi < 0])

    g, v = tools.index("geNomad"), tools.index("VIBRANT")
    rho = lambda a, b: spearman_rows(np.atleast_2d(a), np.atleast_2d(b))
    summary = [["quantity", "value", "CI_low", "CI_high", "note"],
               ["viral_fraction_primary_noE5", r4(ppi[0]), *map(r4, ci(ppib[:, 0])), "participant blocks"],
               ["viral_fraction_MiTCH", r4(mpi[0]), *map(r4, ci(mpib[:, 0])), "sample blocks"],
               ["rho_MCC", r4(rho(m_pt, p_pt)[0]), *map(r4, ci(rho(m_b, p_b))), "both cohorts resampled"],
               ["rho_TPR_minus_FPR", r4(rho(mT - mF, pT - pF)[0]),
                *map(r4, ci(rho(mTb - mFb, pTb - pFb))), "does not depend on viral fraction"],
               ["rho_MCC_MiTCH_at_primary_fraction", r4(rho(m_at, p_pt)[0]),
                *map(r4, ci(rho(m_atb, p_b))), "MiTCH TPR/FPR re-scored at the primary viral fraction"],
               ["n_tools_MiTCH_lower_observed", int((m_pt < p_pt).sum()), "", "", f"of {len(tools)}"],
               ["n_tools_MiTCH_lower_at_primary_fraction", int((m_at < p_pt - 1e-12).sum()), "", "", f"of {len(tools)}"],
               ["median_delta_MCC_observed", r4(np.median(m_pt - p_pt)),
                *map(r4, ci(np.median(m_b - p_b, axis=1))), "MiTCH minus primary"],
               ["median_delta_MCC_at_primary_fraction", r4(np.median(d_at)),
                *map(r4, ci(np.median(d_atb, axis=1))), "MiTCH re-scored at primary fraction minus primary"]]
    for name, T, F, Tb, Fb in (("primary", pT, pF, pTb, pFb), ("MiTCH", mT, mF, mTb, mFb)):
        pt = crossover(T[[g]], F[[g]], T[[v]], F[[v]])[0]
        b = crossover(Tb[:, g], Fb[:, g], Tb[:, v], Fb[:, v])
        summary.append([f"crossover_fraction_geNomad_reaches_VIBRANT_{name}_rates", r4(pt),
                        *map(r4, ci(b)), f"grid {GRID[0]}-{GRID[-1]}; draws without crossover: {int(np.isnan(b).sum())}"])

    for fname, table in (("cross_cohort_rates.tsv", rows), ("cross_cohort_summary.tsv", summary)):
        with open(OUT / fname, "w", newline="") as fh:
            csv.writer(fh, delimiter="\t").writerows(table)
        shutil.copyfile(OUT / fname, TABLES / EXPORT[fname])
    for r in summary:
        print("\t".join(str(c) for c in r))


if __name__ == "__main__":
    main()
