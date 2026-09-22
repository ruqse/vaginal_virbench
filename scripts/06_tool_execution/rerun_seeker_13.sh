#!/usr/bin/env bash
# Re-run Seeker (predict-metagenome) on all 13 real-metagenome assemblies.
# Fixes vs the original failed run:
#   (1) --bind the real project path (container only auto-mounts $HOME/tmp/CWD;
#       the project path may be a symlink the container can't resolve).
#   (2) Seeker asserts every input sequence >= 200 bp, so filter the assembly to
#       contigs >= 1500 bp first (this is exactly the Gate-0 benchmark-ready set,
#       and also makes the run fast).
set -uo pipefail
PROJ="$(cd "$(dirname "$0")/../.." && pwd -P)"
IMG="${PROJ}/What_the_Phage/singularity_images/seeker_0.1.img"
MINLEN="${MINLEN:-1500}"
SAMPLES=(UC028_V2 UC055_V1 UC055_V2 UC062_V2 UC065_V2 UC074_V2 UC084_V2 \
         UC093_V2 UC093_V3 UC096_V2 UC115_V2 UC139_V2 UC164_V2)
echo "[$(date)] Seeker re-run START (${#SAMPLES[@]} samples, >=${MINLEN}bp), PROJ=$PROJ"
for s in "${SAMPLES[@]}"; do
  FASTA="${PROJ}/results/test_real/spades/${s}_contigs.fasta"
  OUT="${PROJ}/results/test_real/full_run/${s}/seeker"
  [[ -f "$FASTA" ]] || { echo "[$(date)] $s SKIP (no fasta)"; continue; }
  mkdir -p "$OUT"
  # filtered input must live under the bound PROJ dir (the container cannot see node-local scratch)
  FILT="${OUT}/_input_ge${MINLEN}.fa"
  # filter by ACTUAL sequence length >= MINLEN
  awk -v ML="$MINLEN" '
    /^>/{ if(h!=""&&length(seq)>=ML){print h; print seq} h=$0; seq=""; next }
    { seq=seq $0 }
    END{ if(h!=""&&length(seq)>=ML){print h; print seq} }' "$FASTA" > "$FILT"
  n=$(grep -c '^>' "$FILT")
  apptainer exec --bind "$PROJ" "$IMG" predict-metagenome "$FILT" \
      > "${OUT}/seeker_results.tsv" 2> "${OUT}/seeker_rerun.log"
  rc=$?
  rows=$(($(wc -l < "${OUT}/seeker_results.tsv") - 1))
  rm -f "$FILT"
  echo "[$(date)] $s DONE rc=$rc input=$n rows=$rows"
done
echo "[$(date)] Seeker re-run COMPLETE"
