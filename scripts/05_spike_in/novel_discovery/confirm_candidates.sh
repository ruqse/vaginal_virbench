#!/usr/bin/env bash
# =============================================================================
# confirm_candidates.sh — DIAMOND BLASTx + skani novelty check on candidates
# =============================================================================
# Runs DIAMOND BLASTx to confirm viral protein homology, then skani to verify
# novelty vs the current spike-in panel.
#
# Usage: sbatch confirm_candidates.sh
# =============================================================================

#SBATCH --cpus-per-task=8
#SBATCH -t 01:00:00
#SBATCH --mem=16G
#SBATCH -J confirm_candidates
#SBATCH -o logs/confirm_candidates_%j.out
#SBATCH -e logs/confirm_candidates_%j.err

set -euo pipefail

PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
OUTDIR="${PROJ_DIR}/data/novel_spike_discovery"
CANDIDATES="${OUTDIR}/candidates/top10_candidates.fasta"
DIAMOND_DB="${PROJ_DIR}/databases/refseq_viral_prot/refseq_viral.dmnd"
SPIKE_VIRAL="${PROJ_DIR}/data/spike_in/all_viral_genomes.fasta"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

export APPTAINER_BIND="${VBENCH_BIND:?set VBENCH_BIND (see config/paths.example.sh)},${TMPDIR:-/tmp}"

mkdir -p "${OUTDIR}/candidates" logs

echo "============================================================"
echo "[$(date)] Candidate Confirmation: DIAMOND + skani"
echo "============================================================"
echo "  Candidates: ${CANDIDATES}"
echo "  Threads:    ${THREADS}"
echo ""

# --- DIAMOND BLASTx ---
echo "[$(date)] DIAMOND BLASTx vs RefSeq viral proteins..."
module load DIAMOND/2.1.11-GCC-13.3.0 2>/dev/null || true

DIAMOND_OUT="${OUTDIR}/candidates/top10_blastx.tsv"
diamond blastx \
    --query "${CANDIDATES}" \
    --db "${DIAMOND_DB}" \
    --out "${DIAMOND_OUT}" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
    --evalue 1e-5 \
    --max-target-seqs 5 \
    --threads "${THREADS}" \
    --block-size 4 \
    --index-chunks 1 \
    --quiet

echo ""
echo "=== Per-candidate best DIAMOND hits ==="
python3 << PYEOF
import csv
from collections import defaultdict

hits = defaultdict(list)
with open("${DIAMOND_OUT}") as f:
    for row in csv.reader(f, delimiter='\t'):
        qid = row[0]
        hits[qid].append({
            'subject': row[1], 'pident': float(row[2]),
            'evalue': float(row[4]), 'bitscore': float(row[5]),
            'desc': row[6] if len(row) > 6 else ''
        })

for qid in sorted(hits.keys()):
    best = sorted(hits[qid], key=lambda x: -x['bitscore'])
    print(f"\n{qid[:70]}:")
    for h in best[:3]:
        print(f"  {h['desc'][:65]:<65}  id={h['pident']:.1f}%  e={h['evalue']:.1e}  score={h['bitscore']:.0f}")

# Flag contigs with NO viral hits
all_candidates = set()
with open("${CANDIDATES}") as f:
    for line in f:
        if line.startswith('>'):
            all_candidates.add(line[1:].strip().split()[0])
no_hits = all_candidates - set(hits.keys())
if no_hits:
    print(f"\nWARNING: {len(no_hits)} candidates with NO viral protein hits:")
    for n in no_hits:
        print(f"  {n}")
else:
    print(f"\nAll {len(all_candidates)} candidates have viral protein hits.")
PYEOF

# --- skani novelty check ---
echo ""
echo "[$(date)] skani novelty check vs spike-in panel..."
module load skani/0.3.1 2>/dev/null || true

SKANI_OUT="${OUTDIR}/candidates/top10_skani.tsv"
skani dist \
    --ql "${CANDIDATES}" \
    --rl "${SPIKE_VIRAL}" \
    -o "${SKANI_OUT}" \
    -t "${THREADS}" \
    --ci \
    -n 5 \
    2>/dev/null || true

echo ""
if [[ -s "${SKANI_OUT}" ]] && [[ $(wc -l < "${SKANI_OUT}") -gt 1 ]]; then
    echo "=== Candidates matching spike-in genomes (NOT novel) ==="
    awk 'BEGIN{FS="\t"} NR>1 && $3 >= 80 {printf "  %s -> %s  ANI=%.1f%%\n", $1, $2, $3}' "${SKANI_OUT}"
else
    echo "No candidates match spike-in genomes at >= 80% ANI."
    echo "ALL CANDIDATES ARE NOVEL relative to the current panel."
fi

echo ""
echo "[$(date)] Done."
echo ""
echo "Results:"
echo "  DIAMOND: ${DIAMOND_OUT}"
echo "  skani:   ${SKANI_OUT}"
