#!/usr/bin/env python3
"""Stratify the multi-evidence (real-metagenome) benchmark by CST and per sample.

Extension to the ViraLM investigation. Tests whether ViraLM's mid-tier rank and the
manuscript's CST-IV precision-collapse finding (established on Track B spike-ins) also
appear in the real-metagenome multi-evidence benchmark, and gives the full per-sample
per-tool table the user asked for.

CST assignments (VERIFIED, not fabricated): from VIRGO2 + VISTA/mgCST
results/mgcst/vista/mgCSTs_5Mar2026.csv, collapsed to traditional CST via the mgCST
groupings in results/tables/table_s2_track_c_by_cst.tsv
(CST-I = mgCST 1,3; CST-III = mgCST 9,10,11; CST-IV = mgCST 16,21,22). Tally
5 CST-I / 4 CST-III / 4 CST-IV matches Methods (manuscript L125); pins UC115_V2=CST-I
and UC093_V3=CST-IV match Track B.

Pure re-evaluation of the locked min-length-1500 pooled predictions. No tool re-runs.
Outputs -> results/test_real/viralm_investigation/analysis/cst/*.tsv + SUMMARY_CST.txt
"""
from __future__ import annotations
import csv, gzip, math
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
POOL = ROOT / "results/test_real/secondary_benchmark_pooled13"
OUT = ROOT / "results/test_real/viralm_investigation/analysis/cst"
OUT.mkdir(parents=True, exist_ok=True)

CST_MAP = {
    "UC084_V2": "CST-I", "UC093_V2": "CST-I", "UC115_V2": "CST-I",
    "UC139_V2": "CST-I", "UC164_V2": "CST-I",
    "UC028_V2": "CST-III", "UC055_V1": "CST-III", "UC065_V2": "CST-III",
    "UC096_V2": "CST-III",
    "UC055_V2": "CST-IV", "UC062_V2": "CST-IV", "UC074_V2": "CST-IV",
    "UC093_V3": "CST-IV",
}
MGCST = {  # for provenance in the output
    "UC028_V2": 11, "UC055_V1": 9, "UC055_V2": 22, "UC062_V2": 22, "UC065_V2": 10,
    "UC074_V2": 16, "UC084_V2": 1, "UC093_V2": 3, "UC093_V3": 21, "UC096_V2": 10,
    "UC115_V2": 1, "UC139_V2": 3, "UC164_V2": 1,
}

_lines: list[str] = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s); _lines.append(s)

def _mcc(tp, fp, tn, fn):
    d = math.sqrt(max(tp+fp, 1)*max(tp+fn, 1)*max(tn+fp, 1)*max(tn+fn, 1))
    return (tp*tn - fp*fn)/d if d > 0 else 0.0

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
        return float("nan")
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

def metrics(yt, yp, sc):
    yt = yt.astype(bool); yp = yp.astype(bool)
    tp = int((yt & yp).sum()); fp = int((~yt & yp).sum())
    tn = int((~yt & ~yp).sum()); fn = int((yt & ~yp).sum())
    prec = tp/(tp+fp) if (tp+fp) else 0.0
    rec = tp/(tp+fn) if (tp+fn) else 0.0
    f1 = 2*prec*rec/(prec+rec) if (prec+rec) else 0.0
    au = _auprc(yt.astype(int), sc) if len(np.unique(sc)) > 1 else float("nan")
    return dict(TP=tp, FP=fp, TN=tn, FN=fn, precision=prec, recall=rec, F1=f1,
                MCC=_mcc(tp, fp, tn, fn), AUPRC=au)

# ---- load master predictions into per-tool aligned arrays ----
df = pd.read_csv(POOL / "pooled_tool_predictions.tsv.gz", sep="\t", compression="gzip",
                 usecols=["sample", "contig_id", "tool", "label", "score", "y_pred"])
tools = sorted(df["tool"].unique())
base = df[["contig_id", "sample", "label"]].drop_duplicates("contig_id").reset_index(drop=True)
base["y_true"] = (base["label"] == "viral").astype(int)
base["cst"] = base["sample"].map(CST_MAP)
idx = {c: i for i, c in enumerate(base["contig_id"])}
n = len(base)
S = np.zeros((n, len(tools))); P = np.zeros((n, len(tools)), int)
for j, t in enumerate(tools):
    sub = df[df["tool"] == t]
    rows = sub["contig_id"].map(idx).to_numpy()
    S[rows, j] = sub["score"].to_numpy(); P[rows, j] = sub["y_pred"].to_numpy()
tj = {t: j for j, t in enumerate(tools)}
yt = base["y_true"].to_numpy()
pr(f"[load] n={n} tools={len(tools)} | CST sizes: " +
   ", ".join(f"{c}={int((base['cst']==c).sum())}" for c in ["CST-I", "CST-III", "CST-IV"]))

# ---- per-CST per-tool metrics ----
cst_rows = []
for cst in ["CST-I", "CST-III", "CST-IV"]:
    m = (base["cst"] == cst).to_numpy()
    nv = int(yt[m].sum()); nn = int((1-yt[m]).sum())
    mm = {t: metrics(yt[m], P[m, tj[t]], S[m, tj[t]]) for t in tools}
    rk = {t: i+1 for i, t in enumerate(sorted(mm, key=lambda x: -mm[x]["MCC"]))}
    for t in tools:
        d = mm[t]
        cst_rows.append(dict(cst=cst, n_total=int(m.sum()), n_viral=nv, n_neg=nn,
                             ratio=f"1:{nn/max(nv,1):.1f}", tool=t, rank_in_cst=rk[t],
                             TP=d["TP"], FP=d["FP"], precision=round(d["precision"], 4),
                             recall=round(d["recall"], 4), F1=round(d["F1"], 4),
                             MCC=round(d["MCC"], 4),
                             AUPRC=round(d["AUPRC"], 4) if d["AUPRC"] == d["AUPRC"] else float("nan")))
pd.DataFrame(cst_rows).to_csv(OUT / "per_cst_tool_metrics.tsv", sep="\t", index=False)

# ---- full per-sample per-tool MCC (+AUPRC) ----
samp_rows = []
for s in sorted(CST_MAP, key=lambda x: (CST_MAP[x], x)):
    m = (base["sample"] == s).to_numpy()
    nv = int(yt[m].sum())
    mm = {t: metrics(yt[m], P[m, tj[t]], S[m, tj[t]]) for t in tools}
    rk = {t: i+1 for i, t in enumerate(sorted(mm, key=lambda x: -mm[x]["MCC"]))}
    for t in tools:
        d = mm[t]
        samp_rows.append(dict(sample=s, cst=CST_MAP[s], mgCST=MGCST[s], n=int(m.sum()),
                              n_viral=nv, tool=t, rank=rk[t], MCC=round(d["MCC"], 4),
                              precision=round(d["precision"], 4), recall=round(d["recall"], 4),
                              AUPRC=round(d["AUPRC"], 4) if d["AUPRC"] == d["AUPRC"] else float("nan")))
sdf = pd.DataFrame(samp_rows)
sdf.to_csv(OUT / "per_sample_tool_metrics_long.tsv", sep="\t", index=False)
# wide MCC matrix (tool x sample) for readability
wide = sdf.pivot(index="tool", columns="sample", values="MCC")
order = ["geNomad", "VirSorter2", "VIBRANT", "Jaeger", "VirSorter", "PPR-Meta", "ViraLM",
         "DeepVirFinder", "MetaPhinder", "VirFinder", "HVSeeker", "Seeker", "Sourmash", "TransGINmer"]
wide = wide.reindex([t for t in order if t in wide.index])
cols = [s for s in sorted(CST_MAP, key=lambda x: (CST_MAP[x], x))]
wide = wide[cols]
wide.to_csv(OUT / "per_sample_MCC_matrix.tsv", sep="\t")

# ---- focused summary ----
pr("\n===== Per-CST ranking (top tools by MCC) =====")
topset = ["geNomad", "VirSorter2", "VIBRANT", "Jaeger", "ViraLM", "PPR-Meta"]
for cst in ["CST-I", "CST-III", "CST-IV"]:
    sub = [r for r in cst_rows if r["cst"] == cst]
    nv = sub[0]["n_viral"]; nn = sub[0]["n_neg"]
    pr(f"\n{cst}  (n={sub[0]['n_total']}, viral={nv}, neg={nn}, {sub[0]['ratio']})")
    for r in sorted(sub, key=lambda x: x["rank_in_cst"])[:7]:
        pr(f"  {r['rank_in_cst']:>2}. {r['tool']:14s} MCC {r['MCC']:.3f}  "
           f"prec {r['precision']:.3f}  rec {r['recall']:.3f}  AUPRC {r['AUPRC']}  (FP {r['FP']})")
    vr = next(r for r in sub if r["tool"] == "ViraLM")
    pr(f"   -> ViraLM rank {vr['rank_in_cst']} (MCC {vr['MCC']:.3f}, prec {vr['precision']:.3f})")

pr("\n===== Precision-collapse check on the REAL-metagenome benchmark (top tools by CST) =====")
pr(f"{'tool':14s} | " + "  ".join(f"{c:>26s}" for c in ["CST-I (prec/MCC)", "CST-III (prec/MCC)", "CST-IV (prec/MCC)"]))
for t in topset:
    cells = []
    for cst in ["CST-I", "CST-III", "CST-IV"]:
        r = next(x for x in cst_rows if x["cst"] == cst and x["tool"] == t)
        cells.append(f"{r['precision']:.3f}/{r['MCC']:.3f}")
    pr(f"{t:14s} | " + "  ".join(f"{c:>26s}" for c in cells))
# mean precision/MCC of the 4 manuscript "top tools" by CST
mtop = ["geNomad", "VirSorter2", "VIBRANT", "ViraLM"]
pr("\nmean over geNomad/VS2/VIBRANT/ViraLM:")
for cst in ["CST-I", "CST-III", "CST-IV"]:
    ps = [next(x for x in cst_rows if x["cst"] == cst and x["tool"] == t)["precision"] for t in mtop]
    ms = [next(x for x in cst_rows if x["cst"] == cst and x["tool"] == t)["MCC"] for t in mtop]
    pr(f"  {cst}: mean precision {np.mean(ps):.3f}  mean MCC {np.mean(ms):.3f}")

pr("\n===== ViraLM rank across CSTs and samples =====")
for cst in ["CST-I", "CST-III", "CST-IV"]:
    ranks = sdf[(sdf.cst == cst) & (sdf.tool == "ViraLM")].sort_values("sample")
    pr(f"  {cst}: " + ", ".join(f"{r['sample']}:r{r['rank']}(MCC {r['MCC']:.2f})"
                                for _, r in ranks.iterrows()))

(OUT / "SUMMARY_CST.txt").write_text("\n".join(_lines) + "\n")
pr(f"\nOutputs in {OUT}")
