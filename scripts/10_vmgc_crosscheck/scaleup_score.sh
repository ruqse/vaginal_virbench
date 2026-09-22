#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -J scaleup_score
#SBATCH -o logs/scaleup_score_%j.out
#SBATCH -e logs/scaleup_score_%j.err
# Final scoring of the 30-sample external validation (pilot MITCH01-05 + scaleup MITCH06-30).
set -uo pipefail
cd ${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}
module load SciPy-bundle/2024.05-gfbf-2024a 2>/dev/null
FR=results/test_real/full_run
GTR=results/test_real/ground_truth
GTNAME=ground_truth_4line_kraken2.tsv

# --- (a) demux pilot 4-line GT into MITCH01..05, symlink pilot tool dir ---
python3 - <<'PY'
import csv
from pathlib import Path
ROOT=Path(__import__("os").environ["VBENCH_ROOT"])
gtp=ROOT/"results/test_real/ground_truth/MITCH_pilot/ground_truth_4line_kraken2.tsv"
rows=list(csv.DictReader(open(gtp),delimiter='\t')); hdr=list(rows[0].keys())
by={}
for r in rows: by.setdefault(r['contig_id'].split('__')[0],[]).append(r)
for s,rs in by.items():
    d=ROOT/f"results/test_real/ground_truth/{s}"; d.mkdir(parents=True,exist_ok=True)
    with open(d/"ground_truth_4line_kraken2.tsv","w",newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=hdr,delimiter='\t'); w.writeheader(); w.writerows(rs)
print(f"[pilot demux] {sorted(by)}")
PY
for s in MITCH01 MITCH02 MITCH03 MITCH04 MITCH05; do
  [[ -e "$FR/$s" ]] || ln -s MITCH_pilot "$FR/$s"
done

# --- (b) symlink scaleup per-sample tool dirs to their chunk dir ---
while IFS=$'\t' read -r s chunk; do
  [[ "$s" == "sample" ]] && continue
  tgt="MITCH_scaleup${chunk}"
  [[ -e "$FR/$s" ]] || ln -s "$tgt" "$FR/$s"
done < .mitch_scaleup_chunkmap.tsv

# --- (c) score all 30 ---
SAMP=$(printf "MITCH%02d " $(seq 1 30))
python scripts/07_evaluation/pooled_multievidence_benchmark.py \
    --samples $SAMP \
    --results-dir $FR --gt-root $GTR --gt-name $GTNAME \
    --min-length 1500 --min-tier 2 \
    --output-dir results/test_real/viralm_investigation/mitch_scaleup30 2>&1 | tail -10

# --- (d) per-block CIs + per-mgCST + Spearman vs discovery ---
python scripts/10_vmgc_crosscheck/scaleup_analysis.py 2>&1 | tail -40
echo "SCORE_DONE"
