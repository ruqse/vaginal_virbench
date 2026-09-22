#!/usr/bin/env bash
# After evidence completes: (1) filter E4 to same-sample spacer matches,
# (2) build the merged 4-line GT, (3) demux GT per-sample by MITCH0X__ prefix,
# (4) extract the Tier 0/1/2 benchmark-subset fasta for the tool stage.
set -euo pipefail
cd ${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}
module load SciPy-bundle/2024.05-gfbf-2024a 2>/dev/null
S=MITCH_scaleup
GT=results/test_real/ground_truth/$S
MERGED=results/test_real/spades/${S}_contigs.fasta
SUBSET=results/test_real/spades/${S}_bench_contigs.fasta

# (0) CRISPR spacer -> contig matching (run_crisprcasfinder only extracts spacers)
if [[ -s "$GT/crispr_spacers.fasta" && ! -f "$GT/evidence_4_crispr.tsv" ]]; then
  module load BLAST+/2.17.0-gompi-2024a 2>/dev/null
  python scripts/04_ground_truth/match_crispr_spacers.py \
    --spacers "$GT/crispr_spacers.fasta" --contigs "$MERGED" \
    --blast-output "$GT/blast_spacers_raw.tsv" \
    --output "$GT/evidence_4_crispr.tsv" --threads 16 2>&1 | tail -6
fi

# (1) same-sample E4 filter (prefix before '__' must match)
if [[ -f "$GT/evidence_4_crispr.tsv" ]]; then
  python3 - "$GT/evidence_4_crispr.tsv" <<'PY'
import csv,sys
p=sys.argv[1]; rows=list(csv.DictReader(open(p),delimiter='\t')); flds=rows[0].keys() if rows else []
def pre(x): return x.split('__')[0]
kept=[r for r in rows if pre(r.get('target_contig_id',''))==pre(r.get('source_array_contig',''))]
import os; os.replace(p, p+'.cross') # keep original as .cross
with open(p,'w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(flds),delimiter='\t'); w.writeheader(); w.writerows(kept)
print(f"  [E4 filter] kept {len(kept)}/{len(rows)} same-sample matches (cross-sample saved to {p}.cross)")
PY
fi

# (2) build merged 4-line GT (E5 N/A for shotgun-only)
python scripts/04_ground_truth/build_ground_truth.py \
  --evidence-dir "$GT" \
  --contigs "$MERGED" \
  --kraken2-taxonomy "$GT/kraken2_output.tsv" \
  --output "$GT/ground_truth_4line_kraken2.tsv" 2>&1 | tail -20

# (3) demux GT per-sample + (4) benchmark subset
python3 - <<'PY'
import csv, os
from pathlib import Path
ROOT=Path(__import__("os").environ["VBENCH_ROOT"])
S="MITCH_scaleup"
gtp=ROOT/f"results/test_real/ground_truth/{S}/ground_truth_4line_kraken2.tsv"
rows=list(csv.DictReader(open(gtp),delimiter='\t')); hdr=rows[0].keys()
by=dict()
for r in rows:
    s=r['contig_id'].split('__')[0]
    by.setdefault(s,[]).append(r)
# per-sample GT dirs
for s,rs in sorted(by.items()):
    d=ROOT/f"results/test_real/ground_truth/{s}"; d.mkdir(parents=True,exist_ok=True)
    with open(d/"ground_truth_4line_kraken2.tsv","w",newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=list(hdr),delimiter='\t'); w.writeheader(); w.writerows(rs)
print(f"[demux] wrote per-sample GT for {len(by)} samples: {sorted(by)}")
# benchmark subset ids = tier in {0,1,2}
bench=set(r['contig_id'] for r in rows if r['tier'] in ('0','1','2'))
print(f"[subset] benchmark contigs (Tier0/1/2): {len(bench)} / {len(rows)} total")
# write subset fasta from merged
src=ROOT/f"results/test_real/spades/{S}_contigs.fasta"
out=ROOT/f"results/test_real/spades/{S}_bench_contigs.fasta"
keep=False; n=0
with open(src) as i, open(out,'w') as o:
    for line in i:
        if line.startswith('>'):
            keep=line[1:].split()[0].strip() in bench
            if keep: n+=1
        if keep: o.write(line)
print(f"[subset] wrote {n} contigs -> {out}")
PY
echo "GT_SUBSET_DONE"
