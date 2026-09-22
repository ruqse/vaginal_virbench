#!/usr/bin/env python3
"""ViraLM real-metagenome investigation — analysis engine (Parts 2-4).

Adversarial re-evaluation of the claim "ViraLM is no longer top-tier on the
real-metagenome multi-evidence benchmark". Pure re-evaluation of EXISTING
per-sample tool outputs; no tool is re-executed. See
docs/viralm_top_tier_investigation_plan.md.

Inputs (all already on disk):
  * results/test_real/secondary_benchmark_pooled13/pooled_tool_predictions.tsv.gz
        long format: sample, contig_id, orig_contig_id, tool, label, score,
        pvalue, threshold, y_pred, emitted_call    (min-length 1500 master)
  * results/test_real/secondary_benchmark_pooled13/pooled_benchmark_ready_ground_truth.tsv
        sample, contig_id, orig_contig_id, tier, category, label, length, evidence_lines
  * results/test_real/ground_truth/<sample>/ground_truth_with_kraken2.tsv
        ... notes  (kraken2 host/bacterial or checkv:0_viral_genes)  -> FP strata
  * results/test_real/viralm_investigation/length_<L>/overall_metrics.tsv & per_sample_summary.tsv
        length sweep produced by 01_length_sweep.sbatch    -> Part 2a

Outputs -> results/test_real/viralm_investigation/analysis/*.tsv  + SUMMARY.txt + summary.json

Metric definitions copied verbatim from pooled_multievidence_benchmark.py /
evaluate_metagenome.py so numbers are identical to the locked benchmark.
"""
from __future__ import annotations
import csv, gzip, json, math, sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
POOL = ROOT / "results/test_real/secondary_benchmark_pooled13"
GTROOT = ROOT / "results/test_real/ground_truth"
INV = ROOT / "results/test_real/viralm_investigation"
OUT = INV / "analysis"
OUT.mkdir(parents=True, exist_ok=True)

SAMPLES = ["UC028_V2", "UC055_V1", "UC055_V2", "UC062_V2", "UC065_V2", "UC074_V2",
           "UC084_V2", "UC093_V2", "UC093_V3", "UC096_V2", "UC115_V2", "UC139_V2",
           "UC164_V2"]
LENGTHS = [0, 500, 1000, 1500, 3000, 5000]
# tools that emit a score for EVERY contig they can process (vs marker tools that
# only report positive calls). Used by the Part-2a coverage audit: for these,
# n_emitted < n_total means the tool's length window dropped contigs (-> conditional).
SCORE_ALL_TOOLS = {"DeepVirFinder", "VirFinder", "PPR-Meta", "MetaPhinder",
                   "HVSeeker", "Seeker", "TransGINmer", "ViraLM", "Jaeger"}
MARKER_TOOLS = {"geNomad", "VirSorter2", "VirSorter", "VIBRANT"}
RNG = np.random.default_rng(20260611)
B_BOOT = 2000
MARGIN = 0.05   # PRE-REGISTERED ΔMCC non-inferiority margin (Part 3b), set before bootstrap

SUMMARY: dict = {}
_lines: list[str] = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s); _lines.append(s)

# ---------------------------------------------------------------- metrics
def _mcc(tp, fp, tn, fn):
    denom = math.sqrt(max(tp + fp, 1) * max(tp + fn, 1) * max(tn + fp, 1) * max(tn + fn, 1))
    return (tp * tn - fp * fn) / denom if denom > 0 else 0.0

def _auprc(y_true, y_scores):
    """Area under the precision-recall curve (average precision; no sklearn).

    Tied scores are collapsed into a single operating point, so the result is
    invariant to input row order. The previous implementation ranked one row at
    a time via argsort, which broke ties by array position; for tools emitting
    few distinct scores (pre-filtered or binary output) that inflated AUPRC
    toward 1.0 whenever the input happened to list positives first.
    """
    y_true = np.asarray(y_true).astype(int)
    y_scores = np.asarray(y_scores, dtype=float)
    n_pos = int(y_true.sum())
    if n_pos == 0:
        return 0.0
    order = np.argsort(-y_scores, kind="mergesort")
    yt = y_true[order]
    ys = y_scores[order]
    tp = np.cumsum(yt)
    fp = np.cumsum(1 - yt)
    # Keep only the last index of each run of tied scores.
    last = np.r_[np.where(np.diff(ys) != 0)[0], len(ys) - 1]
    tp = tp[last]
    fp = fp[last]
    precision = tp / np.maximum(tp + fp, 1)
    recall = np.r_[0.0, tp / n_pos]
    return float(np.sum(np.diff(recall) * precision))

def confusion(yt, yp):
    yt = yt.astype(bool); yp = yp.astype(bool)
    tp = int(np.sum(yt & yp)); fp = int(np.sum(~yt & yp))
    tn = int(np.sum(~yt & ~yp)); fn = int(np.sum(yt & ~yp))
    return tp, fp, tn, fn

def metrics(yt, yp, sc):
    tp, fp, tn, fn = confusion(yt, yp)
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    au = _auprc(yt.astype(int), sc) if len(np.unique(sc)) > 1 else float("nan")
    return dict(TP=tp, FP=fp, TN=tn, FN=fn, precision=prec, recall=rec, F1=f1,
                MCC=_mcc(tp, fp, tn, fn), AUPRC=au)

def ranks_by(d, key="MCC"):
    order = sorted(d, key=lambda t: -d[t][key])
    return {t: i + 1 for i, t in enumerate(order)}

# ---------------------------------------------------------------- load master
def load():
    pr("[load] reading pooled predictions (min-length 1500 master) ...")
    df = pd.read_csv(POOL / "pooled_tool_predictions.tsv.gz", sep="\t",
                     compression="gzip",
                     usecols=["sample", "contig_id", "orig_contig_id", "tool",
                              "label", "score", "y_pred"])
    tools = sorted(df["tool"].unique())
    # one row per contig, fixed order
    base = (df[["contig_id", "sample", "orig_contig_id", "label"]]
            .drop_duplicates("contig_id").reset_index(drop=True))
    base["y_true"] = (base["label"] == "viral").astype(int)
    idx = {c: i for i, c in enumerate(base["contig_id"])}
    n = len(base)
    S = np.zeros((n, len(tools)), float)   # scores
    P = np.zeros((n, len(tools)), int)     # y_pred at default thresholds
    for j, t in enumerate(tools):
        sub = df[df["tool"] == t]
        rows = sub["contig_id"].map(idx).to_numpy()
        S[rows, j] = sub["score"].to_numpy()
        P[rows, j] = sub["y_pred"].to_numpy()
    # GT attributes (tier, category, evidence_lines, length)
    gt = pd.read_csv(POOL / "pooled_benchmark_ready_ground_truth.tsv", sep="\t")
    gt = gt.set_index("contig_id").reindex(base["contig_id"])
    base["tier"] = gt["tier"].to_numpy()
    base["category"] = gt["category"].to_numpy()
    base["evidence_lines"] = gt["evidence_lines"].fillna("-").to_numpy()
    base["length"] = gt["length"].to_numpy()
    # per-sample notes -> FP strata, joined by (sample, orig_contig_id)
    notes = {}
    for s in SAMPLES:
        f = GTROOT / s / "ground_truth_with_kraken2.tsv"
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                notes[(s, r["contig_id"])] = r.get("notes", "") or ""
    def stratum(note):
        if "Homo sapiens" in note or "9606" in note:
            return "kraken2_human"
        if note.startswith("kraken2:"):
            return "kraken2_bacterial"
        if note.startswith("checkv"):
            return "checkv_0_viral_genes"
        return "other_or_blank"
    base["note"] = [notes.get((s, o), "") for s, o in
                    zip(base["sample"], base["orig_contig_id"])]
    base["neg_stratum"] = base["note"].map(stratum)
    pr(f"[load] n_contigs={n}  tools={len(tools)}  "
       f"viral={int(base['y_true'].sum())}  neg={int((1-base['y_true']).sum())}")
    return base, tools, S, P, {t: j for j, t in enumerate(tools)}

# ---------------------------------------------------------------- Part 2a
def part2a_length_coverage(tools):
    pr("\n========== PART 2a: length-policy sweep + coverage audit ==========")
    rows = []; cov_rows = []; viralm_track = []
    for L in LENGTHS:
        d = INV / f"length_{L}"
        om = d / "overall_metrics.tsv"
        ps = d / "per_sample_summary.tsv"
        if not om.exists():
            pr(f"[2a] WARN missing {om} (length sweep not finished?)"); continue
        m = pd.read_csv(om, sep="\t")
        mcc = dict(zip(m["tool"], m["MCC"]))
        rk = {t: i + 1 for i, t in enumerate(sorted(mcc, key=lambda t: -mcc[t]))}
        n_total = int(m["n_total"].iloc[0])
        emit = dict(zip(m["tool"], m["n_emitted"]))
        # coverage from per_sample_summary matched_<tool> (sum over samples)
        psdf = pd.read_csv(ps, sep="\t")
        for t in tools:
            col = f"matched_{t}"
            matched = int(psdf[col].sum()) if col in psdf else emit.get(t, 0)
            frac = matched / n_total if n_total else 0.0
            cov_rows.append(dict(min_length=L, tool=t, n_total=n_total,
                                 n_matched=matched, emitted_fraction=round(frac, 4),
                                 score_all=t in SCORE_ALL_TOOLS))
        # coverage_ok for this length: every SCORE_ALL tool >= 0.99 emitted
        sa_fracs = {t: (int(psdf[f"matched_{t}"].sum()) / n_total)
                    for t in SCORE_ALL_TOOLS if f"matched_{t}" in psdf}
        viralm_frac = sa_fracs.get("ViraLM", float("nan"))
        jaeger_frac = sa_fracs.get("Jaeger", float("nan"))
        cov_ok = all(v >= 0.99 for v in sa_fracs.values())
        viralm_track.append(dict(min_length=L, n_total=n_total,
                                 viralm_MCC=round(mcc.get("ViraLM", float("nan")), 4),
                                 viralm_rank=rk.get("ViraLM"),
                                 viralm_AUPRC=round(float(m.loc[m.tool == "ViraLM", "AUPRC"].iloc[0]), 4),
                                 viralm_emitted_frac=round(viralm_frac, 4),
                                 jaeger_emitted_frac=round(jaeger_frac, 4),
                                 coverage_ok=cov_ok,
                                 top4=", ".join(f"{t}:{mcc[t]:.3f}" for t in
                                                sorted(mcc, key=lambda x: -mcc[x])[:4])))
        rows.append((L, mcc, rk))
    pd.DataFrame(cov_rows).to_csv(OUT / "part2a_coverage_audit.tsv", sep="\t", index=False)
    vt = pd.DataFrame(viralm_track)
    vt.to_csv(OUT / "part2a_viralm_rank_by_length.tsv", sep="\t", index=False)
    pr("ViraLM rank / MCC by length cutoff (coverage_ok flags interpretable cutoffs):")
    for r in viralm_track:
        flag = "" if r["coverage_ok"] else "  <-- CONDITIONAL (uneven coverage)"
        pr(f"  L>={r['min_length']:<5} rank {r['viralm_rank']:>2}  MCC {r['viralm_MCC']:.3f}"
           f"  AUPRC {r['viralm_AUPRC']:.3f}  ViraLM_cov {r['viralm_emitted_frac']:.3f}"
           f"  Jaeger_cov {r['jaeger_emitted_frac']:.3f}{flag}")
        pr(f"          top4: {r['top4']}")
    SUMMARY["part2a"] = viralm_track
    return vt

# ---------------------------------------------------------------- Part 2b
def part2b_persample(base, tools, S, P, tj):
    pr("\n========== PART 2b: per-sample ranking consistency ==========")
    yt = base["y_true"].to_numpy()
    rows = []; viralm_ranks = []
    for s in SAMPLES:
        mask = (base["sample"] == s).to_numpy()
        if mask.sum() == 0:
            continue
        mm = {}
        for t in tools:
            j = tj[t]
            mm[t] = metrics(yt[mask], P[mask, j], S[mask, j])
        rk = ranks_by(mm, "MCC")
        viralm_ranks.append(rk["ViraLM"])
        n_pos = int(yt[mask].sum())
        rows.append(dict(sample=s, n=int(mask.sum()), viral=n_pos,
                         viralm_MCC=round(mm["ViraLM"]["MCC"], 4),
                         viralm_rank=rk["ViraLM"],
                         best_tool=min(rk, key=rk.get),
                         genomad_rank=rk["geNomad"], vs2_rank=rk["VirSorter2"]))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "part2b_viralm_persample_rank.tsv", sep="\t", index=False)
    vr = np.array(viralm_ranks)
    pr(f"ViraLM per-sample rank (n=13): median {np.median(vr):.0f}  "
       f"mean {vr.mean():.1f}  min {vr.min()}  max {vr.max()}")
    pr(f"  distribution: " + ", ".join(f"{s}:r{r}" for s, r in zip(SAMPLES, viralm_ranks)))
    pr(f"  samples where ViraLM top-3: {int((vr<=3).sum())}/13 ; "
       f"top-tier(rank1): {int((vr==1).sum())}/13")
    uc028_rank = rows[0]["viralm_rank"]
    pr(f"  UC028_V2 (original single-sample) ViraLM rank @>=1500: {uc028_rank} "
       f"(median across samples {np.median(vr):.0f}) -> "
       f"{'UC028 NOT atypically favourable' if uc028_rank >= np.median(vr) else 'UC028 somewhat favourable to ViraLM'}")
    SUMMARY["part2b"] = dict(ranks=viralm_ranks, median=float(np.median(vr)),
                             mean=float(vr.mean()), min=int(vr.min()), max=int(vr.max()),
                             uc028_rank=int(uc028_rank))
    return df

# ---------------------------------------------------------------- Part 2c
def part2c_classratio(base, tools, S, P, tj):
    pr("\n========== PART 2c: class-ratio robustness ==========")
    yt = base["y_true"].to_numpy()
    pos = np.where(yt == 1)[0]; neg = np.where(yt == 0)[0]
    n_pos = len(pos)
    pr(f"positives={n_pos}  negatives={len(neg)}  (full ratio 1:{len(neg)/n_pos:.1f})")
    rows = []
    for ratio in [1, 2, 4, 8]:
        n_neg_target = min(len(neg), ratio * n_pos)
        B = 200
        accum = {t: [] for t in tools}
        for _ in range(B):
            negs = RNG.choice(neg, size=n_neg_target, replace=False)
            sel = np.concatenate([pos, negs])
            ytb = yt[sel]
            for t in tools:
                j = tj[t]
                tp, fp, tn, fn = confusion(ytb, P[sel, j])
                accum[t].append(_mcc(tp, fp, tn, fn))
        meanmcc = {t: float(np.mean(accum[t])) for t in tools}
        rk = {t: i + 1 for i, t in enumerate(sorted(meanmcc, key=lambda x: -meanmcc[x]))}
        rows.append(dict(ratio=f"1:{ratio}", n_neg=n_neg_target,
                         viralm_MCC=round(meanmcc["ViraLM"], 4), viralm_rank=rk["ViraLM"],
                         genomad_rank=rk["geNomad"], vs2_rank=rk["VirSorter2"],
                         vibrant_rank=rk["VIBRANT"], jaeger_rank=rk["Jaeger"]))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "part2c_class_ratio.tsv", sep="\t", index=False)
    for r in rows:
        pr(f"  ratio {r['ratio']:>4}: ViraLM rank {r['viralm_rank']} (MCC {r['viralm_MCC']:.3f}) "
           f"| geNomad r{r['genomad_rank']} VS2 r{r['vs2_rank']} VIBRANT r{r['vibrant_rank']} Jaeger r{r['jaeger_rank']}")
    # AUPRC (threshold-free) ranking for reference
    au = {}
    for t in tools:
        j = tj[t]
        au[t] = _auprc(yt, S[:, j]) if len(np.unique(S[:, j])) > 1 else float("nan")
    aurk = {t: i + 1 for i, t in enumerate(sorted(au, key=lambda x: -(au[x] if au[x] == au[x] else -1)))}
    pr(f"  AUPRC (threshold-free, full set): ViraLM rank {aurk['ViraLM']} (AUPRC {au['ViraLM']:.3f}); "
       f"Jaeger rank {aurk['Jaeger']} (AUPRC {au['Jaeger']:.3f})")
    SUMMARY["part2c"] = dict(rows=rows, viralm_auprc_rank=int(aurk["ViraLM"]))
    return df

# ---------------------------------------------------------------- Part 2d
def part2d_threshold(base, tools, S, P, tj):
    pr("\n========== PART 2d: ViraLM threshold sweep (+ best-MCC for all) ==========")
    yt = base["y_true"].to_numpy()
    # ViraLM detailed sweep
    j = tj["ViraLM"]; sc = S[:, j]
    grid = np.round(np.linspace(0.05, 0.95, 19), 2)
    rows = []
    for thr in grid:
        yp = (sc >= thr).astype(int)
        tp, fp, tn, fn = confusion(yt, yp)
        rows.append(dict(threshold=float(thr), TP=tp, FP=fp,
                         precision=round(tp/(tp+fp), 4) if (tp+fp) else 0.0,
                         recall=round(tp/(tp+fn), 4) if (tp+fn) else 0.0,
                         MCC=round(_mcc(tp, fp, tn, fn), 4)))
    pd.DataFrame(rows).to_csv(OUT / "part2d_viralm_threshold_sweep.tsv", sep="\t", index=False)
    best = max(rows, key=lambda r: r["MCC"])
    default = [r for r in rows if abs(r["threshold"] - 0.5) < 1e-9][0]
    pr(f"  ViraLM default(0.5) MCC {default['MCC']:.3f} ; best-MCC {best['MCC']:.3f} @ thr {best['threshold']}")
    # best-achievable MCC for every tool (fine grid over that tool's own scores)
    comp = []
    for t in tools:
        jj = tj[t]; s = S[:, jj]
        uq = np.unique(s)
        if len(uq) <= 1:
            comp.append(dict(tool=t, default_MCC=round(metrics(yt, P[:, jj], s)["MCC"], 4),
                             best_MCC=float("nan"), best_thr=float("nan"))); continue
        cand = np.unique(np.quantile(uq, np.linspace(0, 1, min(len(uq), 200))))
        bm, bt = -2, None
        for thr in cand:
            tp, fp, tn, fn = confusion(yt, (s >= thr).astype(int))
            mm = _mcc(tp, fp, tn, fn)
            if mm > bm:
                bm, bt = mm, thr
        comp.append(dict(tool=t, default_MCC=round(metrics(yt, P[:, jj], s)["MCC"], 4),
                         best_MCC=round(bm, 4), best_thr=round(float(bt), 4)))
    cdf = pd.DataFrame(comp).sort_values("best_MCC", ascending=False)
    cdf.to_csv(OUT / "part2d_best_mcc_all_tools.tsv", sep="\t", index=False)
    brk = {r["tool"]: i+1 for i, r in enumerate(cdf.to_dict("records"))}
    pr(f"  ViraLM best-MCC rank among all tools' best-MCC: {brk['ViraLM']}")
    pr(f"  (geNomad best {cdf[cdf.tool=='geNomad']['best_MCC'].iloc[0]}, "
       f"VS2 best {cdf[cdf.tool=='VirSorter2']['best_MCC'].iloc[0]}, "
       f"VIBRANT best {cdf[cdf.tool=='VIBRANT']['best_MCC'].iloc[0]})")
    SUMMARY["part2d"] = dict(viralm_default_mcc=default["MCC"], viralm_best_mcc=best["MCC"],
                             viralm_best_thr=best["threshold"], viralm_best_rank=int(brk["ViraLM"]))
    return cdf

# ---------------------------------------------------------------- Part 3
def part3_bootstrap(base, tools, S, P, tj):
    pr("\n========== PART 3: nested sample-block bootstrap + paired deltas ==========")
    yt = base["y_true"].to_numpy()
    smp = base["sample"].to_numpy()
    # precompute per-sample class index arrays
    blocks = {}
    for s in SAMPLES:
        m = np.where(smp == s)[0]
        blocks[s] = (m[yt[m] == 1], m[yt[m] == 0])
    nt = len(tools)
    point = {t: metrics(yt, P[:, tj[t]], S[:, tj[t]]) for t in tools}

    def run_bootstrap(mode):
        mcc_b = np.full((B_BOOT, nt), np.nan)
        au_b = np.full((B_BOOT, nt), np.nan)
        n_all = len(yt)
        for b in range(B_BOOT):
            if mode == "nested":
                chosen = RNG.choice(SAMPLES, size=len(SAMPLES), replace=True)
                parts = []
                for s in chosen:
                    pos, neg = blocks[s]
                    if len(pos): parts.append(RNG.choice(pos, len(pos), replace=True))
                    if len(neg): parts.append(RNG.choice(neg, len(neg), replace=True))
                sel = np.concatenate(parts)
            else:  # naive contig bootstrap (sensitivity only)
                sel = RNG.choice(n_all, n_all, replace=True)
            ytb = yt[sel]
            for j, t in enumerate(tools):
                tp, fp, tn, fn = confusion(ytb, P[sel, tj[t]])
                mcc_b[b, j] = _mcc(tp, fp, tn, fn)
                scb = S[sel, tj[t]]
                au_b[b, j] = _auprc(ytb.astype(int), scb) if len(np.unique(scb)) > 1 else np.nan
        return mcc_b, au_b

    rows = []; boot = {}
    for mode in ["nested", "naive"]:
        mcc_b, au_b = run_bootstrap(mode)
        boot[mode] = (mcc_b, au_b)
        for j, t in enumerate(tools):
            rows.append(dict(
                mode=mode, tool=t,
                MCC=round(point[t]["MCC"], 4),
                MCC_lo=round(float(np.nanpercentile(mcc_b[:, j], 2.5)), 4),
                MCC_hi=round(float(np.nanpercentile(mcc_b[:, j], 97.5)), 4),
                AUPRC=round(point[t]["AUPRC"], 4) if point[t]["AUPRC"] == point[t]["AUPRC"] else float("nan"),
                AUPRC_lo=round(float(np.nanpercentile(au_b[:, j], 2.5)), 4),
                AUPRC_hi=round(float(np.nanpercentile(au_b[:, j], 97.5)), 4)))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "part3a_bootstrap_CIs.tsv", sep="\t", index=False)
    pr("Nested sample-block bootstrap 95% CIs (MCC):")
    for r in sorted([x for x in rows if x["mode"] == "nested"], key=lambda x: -x["MCC"]):
        pr(f"  {r['tool']:14s} MCC {r['MCC']:.3f} [{r['MCC_lo']:.3f}, {r['MCC_hi']:.3f}]"
           f"   AUPRC {r['AUPRC'] if r['AUPRC']==r['AUPRC'] else 'NA'} "
           f"[{r['AUPRC_lo']:.3f}, {r['AUPRC_hi']:.3f}]")

    # ---- Part 3b paired deltas (nested replicates) ----
    pr(f"\n--- Part 3b: paired ΔMCC / ΔAUPRC (PRE-REGISTERED margin = {MARGIN} MCC units) ---")
    mcc_b, au_b = boot["nested"]
    vj = tools.index("ViraLM")
    drows = []
    for comp in ["VirSorter2", "VIBRANT", "geNomad", "Jaeger"]:
        cj = tools.index(comp)
        dm = mcc_b[:, vj] - mcc_b[:, cj]
        da = au_b[:, vj] - au_b[:, cj]
        dm_pt = point["ViraLM"]["MCC"] - point[comp]["MCC"]
        da_pt = point["ViraLM"]["AUPRC"] - point[comp]["AUPRC"]
        p_viralm_ge = float(np.mean(dm >= 0))
        p_below_margin = float(np.mean(dm <= -MARGIN))
        # verdict on MCC
        hi = np.nanpercentile(dm, 97.5)
        if hi < -MARGIN:
            verdict = "ViraLM MEANINGFULLY below (95% CI upper < -margin)"
        elif np.nanpercentile(dm, 2.5) > MARGIN:
            verdict = "ViraLM meaningfully ABOVE"
        else:
            verdict = "within-tier / not meaningfully separated"
        drows.append(dict(comparison=f"ViraLM-{comp}",
                          dMCC=round(dm_pt, 4),
                          dMCC_lo=round(float(np.nanpercentile(dm, 2.5)), 4),
                          dMCC_hi=round(float(np.nanpercentile(dm, 97.5)), 4),
                          P_ViraLM_ge=round(p_viralm_ge, 4),
                          P_below_margin=round(p_below_margin, 4),
                          dAUPRC=round(da_pt, 4),
                          dAUPRC_lo=round(float(np.nanpercentile(da, 2.5)), 4),
                          dAUPRC_hi=round(float(np.nanpercentile(da, 97.5)), 4),
                          MCC_verdict=verdict))
    ddf = pd.DataFrame(drows)
    ddf.to_csv(OUT / "part3b_paired_deltas.tsv", sep="\t", index=False)
    for r in drows:
        pr(f"  {r['comparison']:18s} ΔMCC {r['dMCC']:+.3f} [{r['dMCC_lo']:+.3f},{r['dMCC_hi']:+.3f}]"
           f"  P(ViraLM≥)={r['P_ViraLM_ge']:.3f}  ΔAUPRC {r['dAUPRC']:+.3f} "
           f"[{r['dAUPRC_lo']:+.3f},{r['dAUPRC_hi']:+.3f}]  -> {r['MCC_verdict']}")
    SUMMARY["part3"] = dict(margin=MARGIN, point_mcc={t: round(point[t]["MCC"], 4) for t in tools},
                            paired=drows)
    return df, ddf

# ---------------------------------------------------------------- Part 4a/4b
def part4ab_fp(base, tools, S, P, tj):
    pr("\n========== PART 4a/4b: precision decomposition + ViraLM FP characterization ==========")
    yt = base["y_true"].to_numpy()
    # 4a decomposition
    dec = []
    for t in ["geNomad", "VirSorter2", "VIBRANT", "Jaeger", "ViraLM"]:
        mm = metrics(yt, P[:, tj[t]], S[:, tj[t]])
        dec.append(dict(tool=t, TP=mm["TP"], FP=mm["FP"], precision=round(mm["precision"], 4),
                        recall=round(mm["recall"], 4), MCC=round(mm["MCC"], 4)))
    pd.DataFrame(dec).to_csv(OUT / "part4a_precision_decomposition.tsv", sep="\t", index=False)
    pr("4a precision/recall decomposition (top tools):")
    for d in dec:
        pr(f"  {d['tool']:12s} TP {d['TP']:4d} FP {d['FP']:5d}  prec {d['precision']:.3f}  rec {d['recall']:.3f}")
    # 4b FP strata (ViraLM)
    vj = tj["ViraLM"]
    fp_mask = (yt == 0) & (P[:, vj] == 1)
    strata = base.loc[fp_mask, "neg_stratum"].value_counts()
    total_fp = int(fp_mask.sum())
    # the full negative pool composition for reference
    neg_comp = base.loc[yt == 0, "neg_stratum"].value_counts()
    rows = []
    for st in ["kraken2_human", "kraken2_bacterial", "checkv_0_viral_genes", "other_or_blank"]:
        fpc = int(strata.get(st, 0))
        negc = int(neg_comp.get(st, 0))
        rows.append(dict(stratum=st, viralm_FP=fpc,
                         FP_pct=round(100*fpc/total_fp, 1) if total_fp else 0.0,
                         neg_pool=negc,
                         FP_rate_in_stratum=round(fpc/negc, 4) if negc else 0.0))
    pd.DataFrame(rows).to_csv(OUT / "part4b_viralm_fp_strata.tsv", sep="\t", index=False)
    confirmed = sum(r["viralm_FP"] for r in rows if r["stratum"] in
                    ("kraken2_human", "kraken2_bacterial"))
    checkv = next(r["viralm_FP"] for r in rows if r["stratum"] == "checkv_0_viral_genes")
    pr(f"ViraLM FP = {total_fp}. By stratum:")
    for r in rows:
        pr(f"  {r['stratum']:22s} {r['viralm_FP']:5d} ({r['FP_pct']:4.1f}%)  "
           f"[FP-rate within stratum {r['FP_rate_in_stratum']:.3f}]")
    pr(f"  -> confirmed host/bacterial (genuine specificity failure): "
       f"{confirmed}/{total_fp} = {100*confirmed/total_fp:.1f}%")
    pr(f"  -> CheckV-0-genes only (weaker negative; GT-FN upper bound): "
       f"{checkv}/{total_fp} = {100*checkv/total_fp:.1f}%")
    SUMMARY["part4b"] = dict(total_fp=total_fp, strata=rows,
                             confirmed_frac=round(confirmed/total_fp, 4),
                             checkv_frac=round(checkv/total_fp, 4))
    return rows

# ---------------------------------------------------------------- Part 4c
def _has(ev, lines):
    toks = set(ev.split(",")) if ev and ev != "-" else set()
    return any(l in toks for l in lines)

def part4c_markerfree(base, tools, S, P, tj):
    pr("\n========== PART 4c: marker-bias / GT fairness re-ranking ==========")
    yt = base["y_true"].to_numpy()
    ev = base["evidence_lines"].to_numpy()
    is_pos = yt == 1
    is_neg = yt == 0
    MARKER = ["E1a", "E1b", "E2"]; COMPFREE = ["E3", "E4", "E5"]
    # positive subsets
    has_compfree = np.array([_has(e, COMPFREE) for e in ev])
    marker_only_pos = is_pos & ~has_compfree            # EXCLUDE these in marker-free variant
    compfree_pos = is_pos & has_compfree
    pr(f"positives total {int(is_pos.sum())}: composition-independent (E3/E4/E5 present) "
       f"{int(compfree_pos.sum())}; marker-only (excluded in marker-free variant) {int(marker_only_pos.sum())}")

    def rerank(keep_mask, label):
        mm = {}
        sub = keep_mask
        for t in tools:
            j = tj[t]
            mm[t] = metrics(yt[sub], P[sub, j], S[sub, j])
        rk_mcc = ranks_by(mm, "MCC")
        # AUPRC rank (NA-safe)
        au = {t: (mm[t]["AUPRC"] if mm[t]["AUPRC"] == mm[t]["AUPRC"] else -1) for t in tools}
        rk_au = {t: i+1 for i, t in enumerate(sorted(au, key=lambda x: -au[x]))}
        return mm, rk_mcc, rk_au

    variants = {
        "primary (all Tier1+2)": is_pos | is_neg,
        "marker-free (E3/E4/E5 pos only)": compfree_pos | is_neg,
        "RCA-only (E5 pos only)": np.array([_has(e, ["E5"]) for e in ev]) & is_pos | is_neg,
        "prophage/CRISPR-only (E3/E4 pos)": np.array([_has(e, ["E3", "E4"]) for e in ev]) & is_pos | is_neg,
    }
    rows = []
    for label, mask in variants.items():
        npos = int((mask & is_pos).sum())
        mm, rk_mcc, rk_au = rerank(mask, label)
        rows.append(dict(variant=label, n_pos=npos, n_neg=int((mask & is_neg).sum()),
                         viralm_MCC=round(mm["ViraLM"]["MCC"], 4), viralm_rank_MCC=rk_mcc["ViraLM"],
                         viralm_AUPRC=round(mm["ViraLM"]["AUPRC"], 4) if mm["ViraLM"]["AUPRC"]==mm["ViraLM"]["AUPRC"] else float("nan"),
                         viralm_rank_AUPRC=rk_au["ViraLM"],
                         genomad_rank=rk_mcc["geNomad"], vs2_rank=rk_mcc["VirSorter2"],
                         vibrant_rank=rk_mcc["VIBRANT"], jaeger_rank_AUPRC=rk_au["Jaeger"],
                         jaeger_AUPRC=round(mm["Jaeger"]["AUPRC"], 4) if mm["Jaeger"]["AUPRC"]==mm["Jaeger"]["AUPRC"] else float("nan")))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "part4c_markerfree_rerank.tsv", sep="\t", index=False)
    pr("Re-ranking under marker-defined vs composition-independent positive sets:")
    for r in rows:
        pr(f"  {r['variant']:34s} (pos={r['n_pos']:4d}) ViraLM MCC-rank {r['viralm_rank_MCC']:>2} "
           f"(MCC {r['viralm_MCC']:.3f}) AUPRC-rank {r['viralm_rank_AUPRC']:>2} | "
           f"geNomad r{r['genomad_rank']} VS2 r{r['vs2_rank']} | Jaeger AUPRC r{r['jaeger_rank_AUPRC']} ({r['jaeger_AUPRC']})")

    # recall by Tier and by evidence stratum (does ViraLM track weak-marker contigs better?)
    pr("\nPer-tool recall by positive stratum (Tier1/Tier2; marker-only vs E3/E4/E5):")
    tier = base["tier"].to_numpy()
    strata = {
        "Tier1": is_pos & (tier == 1),
        "Tier2": is_pos & (tier == 2),
        "marker-only pos": marker_only_pos,
        "comp-indep pos (E3/E4/E5)": compfree_pos,
        "RCA(E5) pos": np.array([_has(e, ["E5"]) for e in ev]) & is_pos,
    }
    rrows = []
    for t in ["geNomad", "VirSorter2", "VIBRANT", "Jaeger", "ViraLM"]:
        j = tj[t]; row = {"tool": t}
        for sname, sm in strata.items():
            n = int(sm.sum())
            rec = float(P[sm, j].sum())/n if n else float("nan")
            row[f"recall_{sname}"] = round(rec, 3)
            row[f"n_{sname}"] = n
        rrows.append(row)
    pd.DataFrame(rrows).to_csv(OUT / "part4c_recall_by_stratum.tsv", sep="\t", index=False)
    for r in rrows:
        pr(f"  {r['tool']:12s} " + "  ".join(
            f"{k.replace('recall_','')}:{r[k]:.3f}" for k in r if k.startswith("recall_")))
    SUMMARY["part4c"] = dict(variants=rows, recall_by_stratum=rrows)
    return df

# ---------------------------------------------------------------- Part 4d
def part4d_tracka_contrast():
    pr("\n========== PART 4d: Track A (clean panel) vs real-metagenome contrast ==========")
    ta = pd.read_csv(ROOT / "results/tables/table2_track_a_overall_metrics_L1500.tsv", sep="\t")
    ta = ta[ta["tool"].isin(["geNomad", "ViraLM", "VirSorter2"])].copy()
    me = pd.read_csv(POOL / "overall_metrics.tsv", sep="\t")
    me = me[me["tool"].isin(["geNomad", "ViraLM", "VirSorter2"])].set_index("tool")
    rows = []
    for _, r in ta.iterrows():
        t = r["tool"]
        rows.append(dict(tool=t,
                         trackA_precision=round(float(r["precision"]), 4),
                         trackA_recall=round(float(r["recall"]), 4),
                         trackA_MCC=round(float(r["MCC"]), 4),
                         realmeta_precision=round(float(me.loc[t, "precision"]), 4),
                         realmeta_recall=round(float(me.loc[t, "recall"]), 4),
                         realmeta_MCC=round(float(me.loc[t, "MCC"]), 4)))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "part4d_tracka_vs_realmeta.tsv", sep="\t", index=False)
    for r in rows:
        pr(f"  {r['tool']:12s} Track A: prec {r['trackA_precision']:.3f} rec {r['trackA_recall']:.3f} "
           f"MCC {r['trackA_MCC']:.3f}  |  real-meta: prec {r['realmeta_precision']:.3f} "
           f"rec {r['realmeta_recall']:.3f} MCC {r['realmeta_MCC']:.3f}")
    SUMMARY["part4d"] = rows
    return df

# ---------------------------------------------------------------- main
def main():
    base, tools, S, P, tj = load()
    part2a_length_coverage(tools)
    part2b_persample(base, tools, S, P, tj)
    part2c_classratio(base, tools, S, P, tj)
    part2d_threshold(base, tools, S, P, tj)
    part3_bootstrap(base, tools, S, P, tj)
    part4ab_fp(base, tools, S, P, tj)
    part4c_markerfree(base, tools, S, P, tj)
    part4d_tracka_contrast()
    (OUT / "SUMMARY.txt").write_text("\n".join(_lines) + "\n")
    (OUT / "summary.json").write_text(json.dumps(SUMMARY, indent=2, default=str))
    pr(f"\nAll analysis outputs in {OUT}")

if __name__ == "__main__":
    main()
