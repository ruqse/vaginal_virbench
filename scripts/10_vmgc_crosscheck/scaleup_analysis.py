#!/usr/bin/env python3
"""Post-scoring analysis of the 30-sample MITCH external validation.

Reads pooled_tool_predictions.tsv.gz (per-sample long format) and produces:
  * per-sample x per-tool MCC matrix
  * pooled MCC with 95% bootstrap CIs (resampling the 30 sample-blocks)
  * Spearman vs the discovery-cohort no_e5 ranking, with bootstrap CI
  * per-mgCST stratified ranking (25 scale-up samples carry mgCST labels)
Deterministic bootstrap (fixed integer seeds; no Math.random equivalent issues).
"""
from __future__ import annotations
import csv, gzip, math
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
OUT  = ROOT/"results/test_real/viralm_investigation/mitch_scaleup30"
PRED = OUT/"pooled_tool_predictions.tsv.gz"
DISC = ROOT/"results/test_real/secondary_benchmark_pooled13_no_e5/overall_metrics.tsv"
# Real VISTA mgCST + Valencia CST (keyed by MITCHxx), independently reproduced.
# Per-community-type analysis is skipped entirely until this file exists
# (never falls back to any other label source).
MGCST_REAL = ROOT/"results/mgcst_mitch30/mgCSTs_with_CSTs_mitch30.csv"

def mcc(tp,fp,tn,fn):
    d=math.sqrt(max(tp+fp,1)*max(tp+fn,1)*max(tn+fp,1)*max(tn+fn,1))
    return (tp*tn-fp*fn)/d if d>0 else 0.0

def spearman(a,b):
    return float(spearmanr(a, b).statistic)

# --- load per-sample confusion: conf[(tool,sample)] = [tp,fp,tn,fn] ---
conf={}; tools=set(); samples=set()
with gzip.open(PRED,'rt') as fh:
    for r in csv.DictReader(fh,delimiter='\t'):
        t=r['tool']; s=r['sample']; tools.add(t); samples.add(s)
        k=(t,s); v=conf.setdefault(k,[0,0,0,0])
        viral=(r['label']=='viral'); pred=(r['y_pred']=='1')
        if viral: v[0]+=pred; v[3]+=(not pred)
        else:     v[1]+=pred; v[2]+=(not pred)
tools=sorted(tools); samples=sorted(samples)

def pooled_mcc(tool, samp_list):
    tp=fp=tn=fn=0
    for s in samp_list:
        a=conf.get((tool,s));
        if a: tp+=a[0]; fp+=a[1]; tn+=a[2]; fn+=a[3]
    return mcc(tp,fp,tn,fn)

# --- pooled MCC point estimate + rank ---
point={t:pooled_mcc(t,samples) for t in tools}
order=sorted(tools,key=lambda t:-point[t])

# --- bootstrap over sample blocks (B resamples, fixed seeds) ---
B=2000; rng=np.random.default_rng(12345)
idx=np.arange(len(samples))
boot_mcc={t:[] for t in tools}
disc_rows={r['tool']:float(r['MCC']) for r in csv.DictReader(open(DISC),delimiter='\t')}
disc_order=sorted(disc_rows,key=lambda t:-disc_rows[t])
disc_rank={t:i for i,t in enumerate(disc_order)}
common=[t for t in order if t in disc_rank]
boot_spear=[]
for _ in range(B):
    samp=[samples[i] for i in rng.integers(0,len(samples),len(samples))]
    m={t:pooled_mcc(t,samp) for t in tools}
    for t in tools: boot_mcc[t].append(m[t])
    boot_spear.append(spearman([m[t] for t in common],[disc_rows[t] for t in common]))

def ci(v): a=np.percentile(v,[2.5,97.5]); return float(a[0]),float(a[1])

# --- write pooled ranking with CIs ---
with open(OUT/"scaleup30_ranking_ci.tsv","w",newline='') as fh:
    w=csv.writer(fh,delimiter='\t'); w.writerow(["rank","tool","MCC","CI_low","CI_high","disc_no_e5_MCC","disc_rank"])
    for i,t in enumerate(order,1):
        lo,hi=ci(boot_mcc[t])
        w.writerow([i,t,round(point[t],4),round(lo,4),round(hi,4),
                    round(disc_rows.get(t,float('nan')),4), disc_rank.get(t,'NA')+1 if t in disc_rank else 'NA'])
sp_point=spearman([point[t] for t in common],[disc_rows[t] for t in common])
sp_lo,sp_hi=ci(boot_spear)

print("=== 30-sample MITCH external validation: pooled MCC (95% block-bootstrap CI) ===")
print(f"{'rk':>2} {'tool':14s} {'MCC':>7} {'95% CI':>17} | {'disc':>6} disc_rk")
for i,t in enumerate(order,1):
    lo,hi=ci(boot_mcc[t])
    print(f"{i:2d} {t:14s} {point[t]:7.3f} [{lo:6.3f},{hi:6.3f}] | {disc_rows.get(t,float('nan')):6.3f} #{disc_rank.get(t,-1)+1 if t in disc_rank else 'NA'}")
print(f"\nSpearman vs discovery no_e5 = {sp_point:.3f}  (95% CI [{sp_lo:.3f}, {sp_hi:.3f}])")

# --- per-community-type stratified ranking (REAL VISTA mgCST + Valencia CST) ---
from collections import defaultdict
strat={'mgCST':{}, 'CST':{}}
if MGCST_REAL.exists():
    for r in csv.DictReader(open(MGCST_REAL)):
        for key in ('mgCST','CST'):
            if r.get(key): strat[key][r['sampleID']]=r[key]
if not strat['mgCST'] and not strat['CST']:
    print("\n=== per-community-type: SKIPPED ===")
    print(f"  Real mgCST/CST table not found at {MGCST_REAL}.")
    print("  Re-run scripts/10_vmgc_crosscheck/scaleup_analysis.py once it exists.")
else:
    for key,fname in (('mgCST','scaleup30_per_mgCST.tsv'),('CST','scaleup30_per_CST.tsv')):
        lab=strat[key]
        if not lab: continue
        bycst=defaultdict(list)
        for s in samples:
            if s in lab: bycst[lab[s]].append(s)
        with open(OUT/fname,"w",newline='') as fh:
            w=csv.writer(fh,delimiter='\t'); w.writerow([key,"n_samples","top1","top2","top3"])
            print(f"\n=== per-{key} top tools (MCC) [real] ===")
            for c in sorted(bycst,key=lambda c:-len(bycst[c])):
                ss=bycst[c]; m={t:pooled_mcc(t,ss) for t in tools}
                top=sorted(tools,key=lambda t:-m[t])[:3]
                w.writerow([c,len(ss)]+[f"{t}:{m[t]:.3f}" for t in top])
                print(f"  {str(c):12s} (n={len(ss)}): " + ", ".join(f"{t} {m[t]:.3f}" for t in top))

# --- per-sample MCC matrix ---
with open(OUT/"scaleup30_per_sample_MCC.tsv","w",newline='') as fh:
    w=csv.writer(fh,delimiter='\t'); w.writerow(["sample"]+order)
    for s in samples:
        w.writerow([s]+[round(mcc(*conf.get((t,s),[0,0,0,0])),4) for t in order])
print(f"\nOutputs in {OUT}")
