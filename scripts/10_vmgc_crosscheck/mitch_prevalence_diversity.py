#!/usr/bin/env python3
"""MiTCH CST-IV versus CST-I/III/V contrast: estimand, prevalence and host diversity.

Table S31 reports the high- (CST-IV, n = 19) minus low-diversity (CST-I/III/V,
n = 11) difference in POOLED MCC (TP/FP/TN/FN summed across samples, then MCC).
The two arms differ in the labelled viral fraction of their evaluated contig
sets, and MCC depends on class prevalence. This script asks how much of the S31
difference survives other summaries, and whether the viral fraction tracks
bacterial diversity.

Outputs (results/test_real/viralm_investigation/mitch_scaleup30/prevalence_diversity/,
copied verbatim to results/tables/):

  Table S39  contrast_estimands.tsv
      Per tool, high-minus-low delta MCC under four estimands, each with a
      sample-block bootstrap 95% CI:
        pooled            reproduces Table S31 exactly (same draws, same MCC)
        mean_per_sample   arithmetic mean of per-sample MCC in each arm
        prev_std_low      pooled MCC with the high arm's TPR and FPR evaluated
                          at the low arm's labelled viral fraction
        prev_std_cohort   both arms' pooled TPR and FPR evaluated at the fixed
                          whole-cohort labelled viral fraction
  Table S40a per_sample_diversity.tsv
      Per sample: CST, arm, bacterial Shannon diversity, richness, Lactobacillus
      fraction, contig counts by label/category, labelled viral fraction.
  Table S40b diversity_correlations.tsv
      Spearman correlations (all 30 samples; CST-IV-B only, n = 18) among
      bacterial Shannon diversity, labelled viral contigs, viral fraction and
      per-sample MCC, with permutation p-values and bootstrap 95% CIs.

Prevalence standardisation: for an arm with sensitivity TPR and false-positive
rate FPR, expected MCC at prevalence pi uses tp = pi*TPR, fn = pi*(1-TPR),
fp = (1-pi)*FPR, tn = (1-pi)*(1-FPR). This holds each arm's TPR and FPR fixed,
i.e. it assumes they would not change with prevalence. It is a sensitivity
analysis of the metric, not a biological adjustment.

Labels are the frozen four-evidence MiTCH benchmark labels on contigs >= 1,500 bp;
the viral fraction is a property of the evaluated contig set, not a measured
biological viral abundance. Bacterial diversity is computed from VIRGO2
taxon-level normalised counts. Correlations are exploratory and unadjusted for
multiple testing. Tools are not rerun; no labels or thresholds change.

    python3 scripts/10_vmgc_crosscheck/mitch_prevalence_diversity.py
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/test_real/viralm_investigation/mitch_scaleup30"
PRED = BASE / "pooled_tool_predictions.tsv.gz"
TRUTH = BASE / "pooled_benchmark_ready_ground_truth.tsv"
CST = ROOT / "results/mgcst_mitch30/mgCSTs_with_CSTs_mitch30.csv"
TAXA = ROOT / "results/mgcst_mitch30/norm_counts_taxa.csv"
S30 = ROOT / "results/tables/table_s30_mitch_per_sample_mcc.tsv"
S31 = ROOT / "results/tables/table_s31_mitch_diversity_contrast.tsv"
OUT = BASE / "prevalence_diversity"
TABLES = ROOT / "results/tables"
EXPORT = {
    "contrast_estimands.tsv": "table_s39_mitch_contrast_estimands.tsv",
    "per_sample_diversity.tsv": "table_s40a_mitch_per_sample_diversity.tsv",
    "diversity_correlations.tsv": "table_s40b_mitch_diversity_correlations.tsv",
}

B = 2000
SEED = 12345           # Table S31 bootstrap seed; draws are replayed in S31 order
PERM = 10000
PERM_SEED = 20260915
RICHNESS_MIN_RELABUND = 0.001


def require(cond, msg):
    if not cond:
        raise ValueError(msg)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


# ------------------------------------------------------------------ MCC
def mcc_s31(tp, fp, tn, fn):
    """Scalar MCC exactly as scaleup_diversity_contrast.py (Table S31)."""
    d = math.sqrt(max(tp + fp, 1) * max(tp + fn, 1) * max(tn + fp, 1) * max(tn + fn, 1))
    return (tp * tn - fp * fn) / d if d > 0 else 0.0


def mcc_at_prevalence(tp, fp, tn, fn, pi):
    """Expected MCC with this arm's TPR and FPR evaluated at viral prevalence pi."""
    pos, neg = tp + fn, fp + tn
    if pos == 0 or neg == 0:
        return 0.0
    tpr, fpr = tp / pos, fp / neg
    etp, efn = pi * tpr, pi * (1 - tpr)
    efp, etn = (1 - pi) * fpr, (1 - pi) * (1 - fpr)
    d = math.sqrt((etp + efp) * (etp + efn) * (etn + efp) * (etn + efn))
    return (etp * etn - efp * efn) / d if d > 0 else 0.0


# ------------------------------------------------------------------ inputs
def load_arms():
    """Sample -> arm in CSV row order, exactly as the S31 script builds its lists."""
    arm = {}
    with open(CST, newline="") as fh:
        for r in csv.DictReader(fh):
            c = (r.get("CST") or "").strip()
            if c:
                arm[r["sampleID"]] = "high" if c.startswith("IV") else "low"
    low = [s for s, a in arm.items() if a == "low"]
    high = [s for s, a in arm.items() if a == "high"]
    return arm, low, high


def load_confusion():
    conf, tools = {}, set()
    with gzip.open(PRED, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            t, s = r["tool"], r["sample"]
            tools.add(t)
            v = conf.setdefault((t, s), [0, 0, 0, 0])
            viral, pred = r["label"] == "viral", r["y_pred"] == "1"
            if viral:
                v[0] += pred
                v[3] += not pred
            else:
                v[1] += pred
                v[2] += not pred
    return conf, sorted(tools)


def summed(conf, tool, samples):
    tp = fp = tn = fn = 0
    for s in samples:
        a = conf.get((tool, s))
        if a:
            tp += a[0]; fp += a[1]; tn += a[2]; fn += a[3]
    return tp, fp, tn, fn


def prevalence(conf, tool, samples):
    tp, fp, tn, fn = summed(conf, tool, samples)
    return (tp + fn) / (tp + fp + tn + fn)


# ------------------------------------------------------------------ contrast
def contrast(conf, tools, low, high):
    all_samples = low + high
    pi_cohort = prevalence(conf, tools[0], all_samples)
    per_sample = {(t, s): mcc_s31(*conf[(t, s)]) for t in tools for s in all_samples}

    def estimands(t, lo, hi):
        cl, ch = summed(conf, t, lo), summed(conf, t, hi)
        pi_low = (cl[0] + cl[3]) / sum(cl)
        return {
            "pooled": (mcc_s31(*cl), mcc_s31(*ch)),
            "mean_per_sample": (float(np.mean([per_sample[(t, s)] for s in lo])),
                                float(np.mean([per_sample[(t, s)] for s in hi]))),
            "prev_std_low": (mcc_s31(*cl), mcc_at_prevalence(*ch, pi_low)),
            "prev_std_cohort": (mcc_at_prevalence(*cl, pi_cohort),
                                mcc_at_prevalence(*ch, pi_cohort)),
        }

    names = ("pooled", "mean_per_sample", "prev_std_low", "prev_std_cohort")
    boot = {t: {n: [] for n in names} for t in tools}
    rng = np.random.default_rng(SEED)
    for _ in range(B):   # identical draw order to scaleup_diversity_contrast.py
        rs_low = [low[i] for i in rng.integers(0, len(low), len(low))]
        rs_high = [high[i] for i in rng.integers(0, len(high), len(high))]
        for t in tools:
            est = estimands(t, rs_low, rs_high)
            for n in names:
                boot[t][n].append(est[n][1] - est[n][0])

    rows = []
    for t in tools:
        est = estimands(t, low, high)
        cl, ch = summed(conf, t, low), summed(conf, t, high)
        row = {"tool": t, "n_low": len(low), "n_high": len(high),
               "viral_frac_low": round((cl[0] + cl[3]) / sum(cl), 4),
               "viral_frac_high": round((ch[0] + ch[3]) / sum(ch), 4),
               "viral_frac_cohort": round(pi_cohort, 4)}
        for n in names:
            lo_ci, hi_ci = np.percentile(boot[t][n], [2.5, 97.5])
            row[f"{n}_MCC_low"] = round(est[n][0], 4)
            row[f"{n}_MCC_high"] = round(est[n][1], 4)
            row[f"{n}_delta"] = round(est[n][1] - est[n][0], 4)
            row[f"{n}_CIl"] = round(float(lo_ci), 4)
            row[f"{n}_CIh"] = round(float(hi_ci), 4)
            row[f"{n}_excludes_0"] = bool(lo_ci > 0 or hi_ci < 0)
        rows.append(row)
    df = pd.DataFrame(rows)
    df["_order"] = df[["pooled_MCC_low", "pooled_MCC_high"]].max(axis=1)
    return df.sort_values("_order", ascending=False).drop(columns="_order"), per_sample


# ------------------------------------------------------------------ diversity
def per_sample_table(conf, tools, arm, per_sample):
    cst = pd.read_csv(CST)[["sampleID", "CST"]].rename(columns={"sampleID": "sample"})
    taxa = pd.read_csv(TAXA).set_index("Sample").fillna(0)
    rel = taxa.div(taxa.sum(axis=1), axis=0)
    logp = np.log(rel.where(rel > 0, 1.0))
    div = pd.DataFrame({
        "bacterial_shannon": -(rel * logp).sum(axis=1),
        "bacterial_richness_ge_0.1pct": (rel >= RICHNESS_MIN_RELABUND).sum(axis=1),
        "lactobacillus_fraction": rel[[c for c in rel.columns if c.startswith("Lactobacillus")]].sum(axis=1),
    }).rename_axis("sample").reset_index()

    gt = pd.read_csv(TRUTH, sep="\t", usecols=["sample", "label", "category", "length"])
    counts = gt.groupby("sample").agg(
        n_contigs=("label", "size"),
        n_viral=("label", lambda x: int((x == "viral").sum())),
        n_negative=("label", lambda x: int((x == "negative").sum())),
    )
    cats = (gt[gt.label == "viral"].groupby(["sample", "category"]).size()
            .unstack(fill_value=0).add_prefix("n_viral_"))
    counts = counts.join(cats).fillna(0).reset_index()
    counts["viral_fraction"] = counts.n_viral / counts.n_contigs

    df = cst.merge(div, on="sample").merge(counts, on="sample")
    df.insert(2, "arm", df["sample"].map({s: ("CST-IV" if a == "high" else "CST-I/III/V")
                                          for s, a in arm.items()}))
    for t in tools:
        df[f"MCC_{t}"] = df["sample"].map(lambda s, t=t: per_sample[(t, s)])
    require(len(df) == 30, f"expected 30 samples, got {len(df)}")
    return df


def avg_rank(v):
    """Average ranks (ties share the mean rank), as in Spearman's rho."""
    _, inv, cnt = np.unique(v, return_inverse=True, return_counts=True)
    upper = np.cumsum(cnt)
    return ((upper - cnt + 1 + upper) / 2.0)[inv]


def pearson(a, b):
    if a.std() == 0 or b.std() == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def spearman(x, y):
    return pearson(avg_rank(x), avg_rank(y))


def correlations(df, tools):
    rng_p = np.random.default_rng(PERM_SEED)
    rng_b = np.random.default_rng(SEED)
    pairs = [("bacterial_shannon", "n_viral"),
             ("bacterial_shannon", "viral_fraction"),
             ("bacterial_shannon", "n_contigs")]
    pairs += [("viral_fraction", f"MCC_{t}") for t in tools]
    pairs += [("bacterial_shannon", f"MCC_{t}") for t in tools]
    rows = []
    for scope, sub in (("all_samples", df), ("CST-IV-B_only", df[df.CST == "IV-B"])):
        n = len(sub)
        boot_idx = rng_b.integers(0, n, (B, n))
        for x, y in pairs:
            xv, yv = sub[x].to_numpy(float), sub[y].to_numpy(float)
            rho = spearman(xv, yv)
            if math.isnan(rho):
                p = lo = hi = float("nan")
            else:
                rx, ry = avg_rank(xv), avg_rank(yv)   # permuting y permutes its ranks
                perm = np.array([pearson(rx, rng_p.permutation(ry)) for _ in range(PERM)])
                p = (np.sum(np.abs(perm) >= abs(rho) - 1e-12) + 1) / (PERM + 1)
                bs = np.array([spearman(xv[i], yv[i]) for i in boot_idx])
                lo, hi = np.nanpercentile(bs, [2.5, 97.5])
            rows.append({"scope": scope, "n_samples": n, "x": x, "y": y,
                         "spearman_rho": round(rho, 4), "rho_CIl": round(float(lo), 4),
                         "rho_CIh": round(float(hi), 4), "perm_p": round(float(p), 4)})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ checks
def check_against_published(contrast_df, per_sample, low, high, tools):
    require((len(low), len(high)) == (11, 19), f"arm sizes {len(low)}/{len(high)} != 11/19")
    s31 = pd.read_csv(S31, sep="\t").set_index("tool")
    c = contrast_df.set_index("tool")
    for t in tools:
        for mine, theirs in (("pooled_MCC_low", "MCC_low"), ("pooled_MCC_high", "MCC_high"),
                             ("pooled_delta", "delta_high_minus_low"),
                             ("pooled_CIl", "delta_CIl"), ("pooled_CIh", "delta_CIh")):
            require(abs(c.loc[t, mine] - s31.loc[t, theirs]) < 6e-5,
                    f"{t} {mine}={c.loc[t, mine]} does not reproduce S31 {theirs}={s31.loc[t, theirs]}")
        require(bool(c.loc[t, "pooled_excludes_0"]) == (str(s31.loc[t, "delta_excludes_0"]) == "True"),
                f"{t} excludes-zero flag differs from S31")
    s30 = pd.read_csv(S30, sep="\t").set_index("sample")
    for t in tools:
        for s in low + high:
            require(abs(round(per_sample[(t, s)], 4) - s30.loc[s, t]) < 6e-5,
                    f"per-sample MCC {t}/{s} does not reproduce S30")
    # prevalence standardisation at the arm's own prevalence must return pooled MCC
    toy = (40, 60, 880, 20)
    require(abs(mcc_at_prevalence(*toy, 60 / 1000) - mcc_s31(*toy)) < 1e-12,
            "prevalence standardisation does not reduce to MCC at native prevalence")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    arm, low, high = load_arms()
    conf, tools = load_confusion()
    require(len(tools) == 14, f"expected 14 tools, got {len(tools)}")

    contrast_df, per_sample = contrast(conf, tools, low, high)
    check_against_published(contrast_df, per_sample, low, high, tools)
    samples_df = per_sample_table(conf, tools, arm, per_sample)
    require(int(samples_df.n_viral.sum()) == 3582 and int(samples_df.n_contigs.sum()) == 75562,
            "label totals differ from the frozen MiTCH benchmark (3,582 viral / 75,562 contigs)")
    corr_df = correlations(samples_df, tools)

    fmt = {"float_format": "%.4f"}
    contrast_df.to_csv(OUT / "contrast_estimands.tsv", sep="\t", index=False)
    samples_df.to_csv(OUT / "per_sample_diversity.tsv", sep="\t", index=False, **fmt)
    corr_df.to_csv(OUT / "diversity_correlations.tsv", sep="\t", index=False)
    for src, dest in EXPORT.items():
        shutil.copyfile(OUT / src, TABLES / dest)

    meta = {
        "script": "scripts/10_vmgc_crosscheck/mitch_prevalence_diversity.py",
        "bootstrap": {"B": B, "seed": SEED, "unit": "sample block, within arm",
                      "draw_order": "replays scaleup_diversity_contrast.py (Table S31)"},
        "permutation": {"n": PERM, "seed": PERM_SEED},
        "richness_min_relative_abundance": RICHNESS_MIN_RELABUND,
        "arms": {"low_CST-I/III/V": low, "high_CST-IV": high},
        "checks": ["pooled estimand reproduces Table S31 point estimates and CIs",
                   "per-sample MCC reproduces Table S30",
                   "prevalence standardisation reduces to MCC at native prevalence",
                   "label totals match frozen benchmark"],
        "inputs_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in (PRED, TRUTH, CST, TAXA, S30, S31)},
        "outputs_sha256": {src: sha256(OUT / src) for src in EXPORT},
    }
    (OUT / "design_metadata.json").write_text(json.dumps(meta, indent=2) + "\n")

    show = ["tool", "pooled_delta", "pooled_CIl", "pooled_CIh",
            "mean_per_sample_delta", "mean_per_sample_CIl", "mean_per_sample_CIh",
            "prev_std_low_delta", "prev_std_low_CIl", "prev_std_low_CIh",
            "prev_std_cohort_delta", "prev_std_cohort_CIl", "prev_std_cohort_CIh"]
    print(contrast_df[show].to_string(index=False))
    print(samples_df.groupby("arm")[["n_viral", "n_contigs"]].sum().assign(
        frac=lambda d: d.n_viral / d.n_contigs).to_string())
    print(corr_df[corr_df.y.isin(["n_viral", "viral_fraction", "n_contigs", "MCC_geNomad", "MCC_VirSorter2",
                                   "MCC_VIBRANT", "MCC_ViraLM", "MCC_Jaeger"])].to_string(index=False))
    print("All checks passed. Exported:", ", ".join(EXPORT.values()))


if __name__ == "__main__":
    main()
