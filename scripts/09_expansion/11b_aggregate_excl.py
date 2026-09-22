#!/usr/bin/env python3
"""Re-score + re-aggregate the Track B multi-background diversity contrast after
applying the exact UC093_V3-style native-phage exclusions (ground_truth_excl.tsv).

Outputs, all _excl:
  results/expansion/track_b/trackB_expanded_per_background_excl.tsv
  results/tables/table_s9b_track_b_multibackground_excl.tsv
The curated manifest contains four CST-I and four CST-IV-B backgrounds.
Existing per-contig labels and predictions are reused; no tool rerun is needed.
Backgrounds are independently resampled within each arm (B=2000, seed 12345),
preserving the original primary bootstrap and parser/tool iteration order.

The per-background table also carries the exact confusion counts (TP, FP, TN, FN).
Table S9b additionally reports the precision and false-positive-rate (FPR)
contrasts with the same background bootstrap, each metric drawing from its own
generator (seed 12345) so the primary MCC intervals are unchanged, and the range
of the precision contrast when each of the eight backgrounds is omitted in turn.
"""
from __future__ import annotations
import csv, sys
from pathlib import Path
from collections import defaultdict
import numpy as np
from track_b_manifest import load_backgrounds, restore_metric_precision, require_successful_executions

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"scripts/07_evaluation"))
from benchmark_track_b import load_ground_truth, evaluate_assembly, TOOL_PARSERS, TOOL_THRESHOLDS

BG_TSV=ROOT/"data/expansion_cohorts/track_b_backgrounds.tsv"
ASM=ROOT/"results/expansion/track_b/assemblies"
BENCH=ROOT/"results/expansion/track_b/benchmark"
COV="10"
OUT_PB=ROOT/"results/expansion/track_b/trackB_expanded_per_background_excl.tsv"
OUT_S9B=ROOT/"results/tables/table_s9b_track_b_multibackground_excl.tsv"

backgrounds = load_backgrounds(BG_TSV)
require_successful_executions(backgrounds, ROOT)  # Fail before replacing any metrics.
backs=[(r["run"],r["stratum"]) for r in backgrounds]

per=[]; mcc_by=defaultdict(lambda:defaultdict(list))
for srr,arm in backs:
    gt_path=ASM/srr/f"assembly_cov{COV}"/"ground_truth_excl.tsv"
    tool_dir=BENCH/f"trackB_{srr}_cov{COV}"
    if not gt_path.exists() or not tool_dir.exists():
        raise FileNotFoundError(f"{srr}: missing excluded-label GT or tool dir")
    gt=load_ground_truth(str(gt_path))
    nv=sum(1 for v in gt.values() if v['label']=='viral')
    nn=sum(1 for v in gt.values() if v['label']=='negative')
    nx=sum(1 for v in gt.values() if v['label']=='excluded')
    res=evaluate_assembly(gt,str(tool_dir),TOOL_PARSERS,TOOL_THRESHOLDS)
    for r in res:
        r = restore_metric_precision(r)
        tp, fp, tn, fn = (int(r[k]) for k in ("TP", "FP", "TN", "FN"))
        if tp + fn != nv or fp + tn != nn:
            raise ValueError(f"{srr} {r['tool']}: confusion counts do not match labels "
                             f"(TP+FN={tp+fn} vs {nv} viral; FP+TN={fp+tn} vs {nn} negative)")
        per.append((srr,arm,r['tool'],r['MCC'],r['precision'],r['recall'],nv,nn,nx,tp,fp,tn,fn))
        mcc_by[r['tool']][arm].append(r['MCC'])
    print(f"  {srr:12s} {arm:9s}: {nv} viral, {nn} negative, {nx} excluded(native), {len(res)} tools")

with open(OUT_PB,'w',newline='') as fh:
    w=csv.writer(fh,delimiter='\t',lineterminator='\n')
    w.writerow(["background","stratum","tool","MCC","precision","recall","n_viral","n_negative","n_excluded_native",
                "TP","FP","TN","FN"])
    for row in per: w.writerow(row)  # Retain full precision for reproducible CIs.

def contrast(mby, B=2000, seed=12345):
    rng=np.random.default_rng(seed); out={}
    for t,armd in mby.items():
        vi=[m for m in armd.get('CST-I',[])]
        viv=armd.get('CST-IV-B',[])
        if len(vi) != 4 or len(viv) != 4:
            raise ValueError(f"{t}: incomplete four-versus-four metrics")
        vi=np.array(vi,float); viv=np.array(viv,float)
        bi=vi[rng.integers(0,len(vi),(B,len(vi)))].mean(1); biv=viv[rng.integers(0,len(viv),(B,len(viv)))].mean(1)
        bd=biv-bi
        out[t]=(vi.mean(),viv.mean(),viv.mean()-vi.mean(),float(np.percentile(bd,2.5)),float(np.percentile(bd,97.5)),len(vi),len(viv))
    return out

full=contrast(mcc_by)
order=sorted(full,key=lambda t:-full[t][0])

# per-arm precision/recall/FPR means over the backgrounds in each arm (backgrounds =
# unit of replication), aggregated from `per`; precision_rel_drop_pct is the
# relative reduction 100*(CST-I - CST-IV-B)/CST-I, backing the manuscript figures.
# FPR = FP/(FP+TN) from the exact confusion counts.
pr_by=defaultdict(lambda:defaultdict(list)); rc_by=defaultdict(lambda:defaultdict(list))
fpr_by=defaultdict(lambda:defaultdict(list)); pr_bg=defaultdict(lambda:defaultdict(list))
for srr,arm,tool,mcc,prec,rec,nv,nn,nx,tp,fp,tn,fn in per:
    pr_by[tool][arm].append(prec); rc_by[tool][arm].append(rec)
    fpr_by[tool][arm].append(fp/(fp+tn) if fp+tn else 0.0)
    pr_bg[tool][arm].append(srr)
def amean(d,t,arm):
    v=d[t].get(arm,[]); return float(np.mean(v)) if v else float('nan')

# Same background bootstrap as MCC; contrast() creates a fresh generator per call,
# so these draws do not alter the primary MCC intervals computed above.
prec_ct=contrast(pr_by)
fpr_ct=contrast(fpr_by)

def loo_range(vals, ids):
    """Precision contrast (CST-IV-B minus CST-I) with each background omitted in turn.
    Returns the smallest and largest of the eight estimates and the background whose
    omission gives the largest (least negative) estimate; blank when all are equal."""
    vi=np.array(vals['CST-I'],float); viv=np.array(vals['CST-IV-B'],float)
    est=[(viv.mean()-np.delete(vi,k).mean(), ids['CST-I'][k]) for k in range(len(vi))]
    est+=[(np.delete(viv,k).mean()-vi.mean(), ids['CST-IV-B'][k]) for k in range(len(viv))]
    lo=min(e for e,_ in est); hi=max(e for e,_ in est)
    return lo, hi, (max(est,key=lambda x:x[0])[1] if hi>lo else "")

with open(OUT_S9B,'w',newline='') as fh:
    w=csv.writer(fh,delimiter='\t',lineterminator='\n')
    w.writerow(["tool","MCC_CST_I_mean","MCC_CST_IVB_mean","delta_IVB_minus_I","delta_CIl","delta_CIh","delta_excludes_0","n_CST_I","n_CST_IVB",
                "precision_CST_I_mean","precision_CST_IVB_mean","precision_rel_drop_pct","recall_CST_I_mean","recall_CST_IVB_mean",
                "precision_delta_IVB_minus_I","precision_delta_CIl","precision_delta_CIh","precision_delta_excludes_0",
                "precision_delta_loo_min","precision_delta_loo_max","precision_delta_loo_max_omitted",
                "FPR_CST_I_mean","FPR_CST_IVB_mean","FPR_delta_IVB_minus_I","FPR_delta_CIl","FPR_delta_CIh","FPR_delta_excludes_0"])
    for t in order:
        d=full[t]
        pi=amean(pr_by,t,'CST-I'); pv=amean(pr_by,t,'CST-IV-B')
        ri=amean(rc_by,t,'CST-I'); rv=amean(rc_by,t,'CST-IV-B')
        reldrop=round(100*(pi-pv)/pi,2) if pi else float('nan')
        p=prec_ct[t]; f=fpr_ct[t]; llo,lhi,lom=loo_range(pr_by[t],pr_bg[t])
        w.writerow([t,round(d[0],4),round(d[1],4),round(d[2],4),round(d[3],4),round(d[4],4),(d[3]>0 or d[4]<0),d[5],d[6],
                    round(pi,4),round(pv,4),reldrop,round(ri,4),round(rv,4),
                    round(p[2],4),round(p[3],4),round(p[4],4),(p[3]>0 or p[4]<0),round(llo,4),round(lhi,4),lom,
                    round(f[0],4),round(f[1],4),round(f[2],4),round(f[3],4),round(f[4],4),(f[3]>0 or f[4]<0)])

print("\n=== Track B: 4 CST-I + 4 CST-IV-B, native-homolog-excluded labels ===")
print("Delta = CST-IV-B minus CST-I; bootstrap B=2000, seed 12345.")
print(f"{'tool':14s} | {'CST-I':>8s} | {'CST-IV-B':>8s} | {'delta [95% CI]':>26s}")
for t in order:
    d=full[t]
    sig='*' if (d[3]>0 or d[4]<0) else ' '
    print(f"{t:14s} | {d[0]:8.3f} | {d[1]:8.3f} | {d[2]:+6.3f} [{d[3]:+6.3f},{d[4]:+6.3f}] {sig}")
print("\n* delta 95% CI excludes 0.")
print("\nPrecision and FPR contrasts (same bootstrap; LOO = eight single-background omissions):")
print(f"{'tool':14s} | {'d precision [95% CI]':>26s} | {'LOO range':>17s} | {'d FPR [95% CI]':>26s}")
for t in order:
    p=prec_ct[t]; f=fpr_ct[t]; llo,lhi,_=loo_range(pr_by[t],pr_bg[t])
    ps='*' if (p[3]>0 or p[4]<0) else ' '; fs='*' if (f[3]>0 or f[4]<0) else ' '
    print(f"{t:14s} | {p[2]:+6.3f} [{p[3]:+6.3f},{p[4]:+6.3f}] {ps} | {llo:+6.3f} to {lhi:+6.3f} | {f[2]:+6.3f} [{f[3]:+6.3f},{f[4]:+6.3f}] {fs}")
print(f"\nPer-background (excl): {OUT_PB}")
print(f"Table S9b (excl):      {OUT_S9B}")
