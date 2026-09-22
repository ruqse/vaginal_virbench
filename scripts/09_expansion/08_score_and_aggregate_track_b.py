#!/usr/bin/env python3
"""Track B expansion: score the eight curated backgrounds and compute the
multi-background CST-I vs CST-IV-B diversity contrast with CIs.

Reuses the published Track B parsers and ground-truth loader (no scoring-logic
change); only the iteration set changes from the original 2 backgrounds to the
eight expansion backgrounds at 10x. Backgrounds are the unit of replication:
per-arm MCC is the mean over the four backgrounds in that arm, and the delta
(CST-IV-B minus CST-I) carries a bootstrap CI resampling backgrounds within arm
(B=2000, seed 12345), mirroring scaleup_diversity_contrast.py.

Outputs:
  results/expansion/track_b/trackB_expanded_per_background.tsv
  results/secondary_tables/table_s9b_track_b_multibackground.tsv

Note: this native-phage-inclusive aggregate is NOT the published Table S9.
The manuscript ships the native-phage-excluded variant written by
11b_aggregate_excl.py. This one is retained for provenance only, which is why
it writes into results/secondary_tables/ rather than alongside the display tables.
"""
from __future__ import annotations
import csv, sys
from pathlib import Path
import numpy as np
from track_b_manifest import load_backgrounds, restore_metric_precision, require_successful_executions

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"scripts/07_evaluation"))
from benchmark_track_b import load_ground_truth, evaluate_assembly, TOOL_PARSERS, TOOL_THRESHOLDS

BG_TSV = ROOT/"data/expansion_cohorts/track_b_backgrounds.tsv"
ASM    = ROOT/"results/expansion/track_b/assemblies"
BENCH  = ROOT/"results/expansion/track_b/benchmark"
COV    = "10"
OUT_PB = ROOT/"results/expansion/track_b/trackB_expanded_per_background.tsv"
OUT_S9B= ROOT/"results/secondary_tables/table_s9b_track_b_multibackground.tsv"

# --- read backgrounds + arm (stratum) ---
backgrounds = load_backgrounds(BG_TSV)
require_successful_executions(backgrounds, ROOT)  # Fail before replacing any metrics.
backs = [(r["run"], r["stratum"]) for r in backgrounds]

# --- score each background (reuse published parsers) ---
per=[]  # rows: background, stratum, tool, MCC, precision, recall, n_viral, n_negative
mcc_by={}  # (tool,arm) -> list of per-background MCC
for srr, arm in backs:
    gt_path = ASM/srr/f"assembly_cov{COV}"/"ground_truth.tsv"
    tool_dir= BENCH/f"trackB_{srr}_cov{COV}"
    if not gt_path.exists() or not tool_dir.exists():
        raise FileNotFoundError(f"{srr}: missing GT or tool dir")
    gt=load_ground_truth(str(gt_path))
    nv=sum(1 for v in gt.values() if v['label']=='viral')
    nn=sum(1 for v in gt.values() if v['label']=='negative')
    res=evaluate_assembly(gt, str(tool_dir), TOOL_PARSERS, TOOL_THRESHOLDS)
    for r in res:
        r = restore_metric_precision(r)
        per.append((srr, arm, r['tool'], r['MCC'], r['precision'], r['recall'], nv, nn))
        mcc_by.setdefault((r['tool'],arm),[]).append(r['MCC'])
    print(f"  {srr:12s} {arm:9s}: {nv} viral, {nn} negative, {len(res)} tools scored")

tools=sorted({t for (t,a) in mcc_by})
arms=("CST-I","CST-IV-B")

# --- write per-background metrics ---
OUT_PB.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_PB,'w',newline='') as fh:
    w=csv.writer(fh,delimiter='\t',lineterminator='\n')
    w.writerow(["background","stratum","tool","MCC","precision","recall","n_viral","n_negative"])
    for row in per:
        w.writerow(row)  # Preserve full precision for reproducible downstream CIs.

# --- multi-background stats: mean per arm + delta, bootstrap over backgrounds ---
B=2000; rng=np.random.default_rng(12345)
def boot_mean(vals):
    vals=np.asarray(vals,float)
    idx=rng.integers(0,len(vals),(B,len(vals)))
    return vals[idx].mean(axis=1)   # B bootstrap means

rows=[]
for t in tools:
    vi = mcc_by.get((t,"CST-I"),[]); viv= mcc_by.get((t,"CST-IV-B"),[])
    mi=float(np.mean(vi)) if vi else float('nan'); miv=float(np.mean(viv)) if viv else float('nan')
    bi=boot_mean(vi); biv=boot_mean(viv); bd=biv-bi
    ci=lambda a:(float(np.percentile(a,2.5)),float(np.percentile(a,97.5)))
    ilo,ihi=ci(bi); ivlo,ivhi=ci(biv); dlo,dhi=ci(bd)
    delta=miv-mi; sig=(dlo>0) or (dhi<0)
    rows.append((t,len(vi),len(viv),mi,ilo,ihi,miv,ivlo,ivhi,delta,dlo,dhi,sig))
rows.sort(key=lambda r:-r[3])  # by CST-I mean (matches original Table S9 ordering)

# per-arm precision/recall means over the backgrounds in each arm (backgrounds =
# unit of replication); precision_rel_drop_pct is the relative reduction
# 100*(CST-I - CST-IV-B)/CST-I.
from collections import defaultdict as _dd
pr_by=_dd(lambda:_dd(list)); rc_by=_dd(lambda:_dd(list))
for srr,arm,tool,mcc,prec,rec,nv,nn in per:
    pr_by[tool][arm].append(prec); rc_by[tool][arm].append(rec)
def amean(d,t,arm):
    v=d[t].get(arm,[]); return float(np.mean(v)) if v else float('nan')

OUT_S9B.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_S9B,'w',newline='') as fh:
    w=csv.writer(fh,delimiter='\t',lineterminator='\n')
    w.writerow(["tool","n_CST_I","n_CST_IVB","MCC_CST_I_mean","CST_I_CIl","CST_I_CIh",
                "MCC_CST_IVB_mean","CST_IVB_CIl","CST_IVB_CIh",
                "delta_IVB_minus_I","delta_CIl","delta_CIh","delta_excludes_0",
                "precision_CST_I_mean","precision_CST_IVB_mean","precision_rel_drop_pct",
                "recall_CST_I_mean","recall_CST_IVB_mean"])
    for r in rows:
        t=r[0]
        pi=amean(pr_by,t,'CST-I'); pv=amean(pr_by,t,'CST-IV-B')
        ri=amean(rc_by,t,'CST-I'); rv=amean(rc_by,t,'CST-IV-B')
        reldrop=round(100*(pi-pv)/pi,2) if pi else float('nan')
        w.writerow([r[0],r[1],r[2]]+[round(x,4) if isinstance(x,float) else x for x in r[3:]]
                   +[round(pi,4),round(pv,4),reldrop,round(ri,4),round(rv,4)])

print(f"\n=== Track B multi-background diversity contrast (4 CST-I + 4 CST-IV-B, 10x) ===")
print(f"Backgrounds = unit of replication; bootstrap B={B}, seed 12345.\n")
print(f"{'tool':14s} {'CST-I mean [95% CI]':>24s} {'CST-IV-B mean [95% CI]':>26s} {'delta [95% CI]':>24s} sig")
for r in rows:
    t,ni,niv,mi,ilo,ihi,miv,ivlo,ivhi,d,dlo,dhi,sig=r
    print(f"{t:14s} {mi:6.3f} [{ilo:6.3f},{ihi:6.3f}] {miv:6.3f} [{ivlo:6.3f},{ivhi:6.3f}] "
          f"{d:+6.3f} [{dlo:+6.3f},{dhi:+6.3f}] {'*' if sig else ''}")
print("\n* delta 95% CI excludes 0.")
print(f"\nPer-background: {OUT_PB}\nTable S9b:      {OUT_S9B}")
