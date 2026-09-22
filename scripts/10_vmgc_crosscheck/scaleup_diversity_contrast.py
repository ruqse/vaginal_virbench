#!/usr/bin/env python3
"""MITCH external validation: low- vs high-diversity tool performance with CIs.

Collapses the verified community state types into low-diversity (Lactobacillus-
dominated: CST I/II/III/V) vs high-diversity (BV-associated: CST IV*) and reports
per-tool pooled MCC in each stratum plus the high-minus-low delta, each with a
block bootstrap 95% CI that resamples sample-blocks WITHIN the stratum
(B=2000, seed 12345; same scheme as scaleup_analysis.py). The delta CI is the
inferential quantity: it tells whether a tool's performance differs by diversity.

This is the observational generalisation of the controlled Track B diversity
contrast (CST-I vs CST-IV backgrounds, n=1 per arm). Here the axis is tested
across 30 real communities (low n, high n printed at runtime).
"""
from __future__ import annotations
import csv, gzip, math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT  = ROOT/"results/test_real/viralm_investigation/mitch_scaleup30"
PRED = OUT/"pooled_tool_predictions.tsv.gz"
CST  = ROOT/"results/mgcst_mitch30/mgCSTs_with_CSTs_mitch30.csv"

def mcc(tp,fp,tn,fn):
    d=math.sqrt(max(tp+fp,1)*max(tp+fn,1)*max(tn+fp,1)*max(tn+fn,1))
    return (tp*tn-fp*fn)/d if d>0 else 0.0

# diversity class from Ravel CST: IV* = high diversity (BV), else Lactobacillus low
def divclass(c):
    c=(c or "").strip()
    return "high" if c.startswith("IV") else ("low" if c else None)

samp_div={}
for r in csv.DictReader(open(CST)):
    d=divclass(r.get("CST"))
    if d: samp_div[r["sampleID"]]=d

conf={}; tools=set()
with gzip.open(PRED,'rt') as fh:
    for r in csv.DictReader(fh,delimiter='\t'):
        t=r['tool']; s=r['sample']; tools.add(t)
        v=conf.setdefault((t,s),[0,0,0,0])
        viral=(r['label']=='viral'); pred=(r['y_pred']=='1')
        if viral: v[0]+=pred; v[3]+=(not pred)
        else:     v[1]+=pred; v[2]+=(not pred)
tools=sorted(tools)
low =[s for s,d in samp_div.items() if d=='low']
high=[s for s,d in samp_div.items() if d=='high']

def pooled(tool, samp_list):
    tp=fp=tn=fn=0
    for s in samp_list:
        a=conf.get((tool,s))
        if a: tp+=a[0];fp+=a[1];tn+=a[2];fn+=a[3]
    return mcc(tp,fp,tn,fn)

# point estimates
pt_low ={t:pooled(t,low)  for t in tools}
pt_high={t:pooled(t,high) for t in tools}

# block bootstrap within each stratum (independent resamples; disjoint samples)
B=2000; rng=np.random.default_rng(12345)
bl={t:[] for t in tools}; bh={t:[] for t in tools}; bd={t:[] for t in tools}
lo_arr=np.array(low); hi_arr=np.array(high)
for _ in range(B):
    rs_low =[low[i]  for i in rng.integers(0,len(low),len(low))]
    rs_high=[high[i] for i in rng.integers(0,len(high),len(high))]
    for t in tools:
        ml=pooled(t,rs_low); mh=pooled(t,rs_high)
        bl[t].append(ml); bh[t].append(mh); bd[t].append(mh-ml)

def ci(v): a=np.percentile(v,[2.5,97.5]); return float(a[0]),float(a[1])

rows=[]
for t in tools:
    llo,lhi=ci(bl[t]); hlo,hhi=ci(bh[t]); dlo,dhi=ci(bd[t])
    sig = (dlo>0) or (dhi<0)   # delta CI excludes 0
    rows.append((t, pt_low[t], llo, lhi, pt_high[t], hlo, hhi,
                 pt_high[t]-pt_low[t], dlo, dhi, sig))
rows.sort(key=lambda r:-max(r[1],r[4]))

with open(OUT/"scaleup30_diversity_contrast.tsv","w",newline='') as fh:
    w=csv.writer(fh,delimiter='\t')
    w.writerow(["tool","MCC_low","low_CIl","low_CIh","MCC_high","high_CIl","high_CIh",
                "delta_high_minus_low","delta_CIl","delta_CIh","delta_excludes_0"])
    for r in rows:
        w.writerow([r[0]]+[round(x,4) if isinstance(x,float) else x for x in r[1:]])

print(f"MITCH low-diversity (Lactobacillus CST I/II/III/V): n={len(low)} samples")
print(f"MITCH high-diversity (BV-associated CST IV*):       n={len(high)} samples")
print(f"Block bootstrap B={B}, seed 12345, resampling sample-blocks within each stratum.\n")
print(f"{'tool':14s} {'low MCC [95% CI]':>22s} {'high MCC [95% CI]':>23s} {'delta [95% CI]':>22s} sig")
for r in rows:
    t,ml,llo,lhi,mh,hlo,hhi,dl,dlo,dhi,sig=r
    print(f"{t:14s} {ml:6.3f} [{llo:6.3f},{lhi:6.3f}] {mh:6.3f} [{hlo:6.3f},{hhi:6.3f}] "
          f"{dl:+6.3f} [{dlo:+6.3f},{dhi:+6.3f}] {'*' if sig else ''}")
print("\n* delta 95% CI excludes 0 (diversity effect resolved in sign for this tool).")
print(f"\nWritten: {OUT/'scaleup30_diversity_contrast.tsv'}")
