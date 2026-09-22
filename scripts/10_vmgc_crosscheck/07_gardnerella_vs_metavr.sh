#!/usr/bin/env bash
# =============================================================================
# 07_gardnerella_vs_metavr.sh — direct sequence-level check of the Track A
#   Gardnerella phage (MW387018.1, vB_Gva_AB1) against MetaVR v5.
# =============================================================================
# Complements the host-label finding (0 of 24.4M MetaVR UViGs carry a
# Gardnerella host prediction): does the Gardnerella phage SEQUENCE have any
# MetaVR v5 match? Reports the best-subject coverage with the corrected
# interval-merged parser (skani is skipped: it OOMs on the 272 GB MetaVR FASTA).
#
# Usage: sbatch scripts/10_vmgc_crosscheck/07_gardnerella_vs_metavr.sh
# =============================================================================

#SBATCH --cpus-per-task=8
#SBATCH -t 02:00:00
#SBATCH --mem=32G
#SBATCH -J gard_metavr
#SBATCH -o logs/gard_metavr_%j.out
#SBATCH -e logs/gard_metavr_%j.err

set -euo pipefail
PROJ="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
cd "${PROJ}"
OUT="results/vmgc_crosscheck/gardnerella_check"
mkdir -p "${OUT}" logs
THREADS="${SLURM_CPUS_PER_TASK:-8}"

module load BLAST+/2.17.0-gompi-2024a 2>/dev/null || module load BLAST+/2.16.0-gompi-2024a 2>/dev/null || true

Q="data/spike_in/gardnerella_phage/MW387018.1.fasta"
RAW="${OUT}/MW387018_vs_metavr.tsv"

echo "[$(date)] BLASTn MW387018.1 (Gardnerella phage vB_Gva_AB1) vs MetaVR v5..."
blastn -query "${Q}" -db databases/metavr_v5/metavr_v5 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
    -evalue 1e-5 -max_target_seqs 25 -num_threads "${THREADS}" \
    -out "${RAW}"
echo "  raw hit lines: $(wc -l < "${RAW}")"

echo "[$(date)] Best-subject coverage (corrected interval-merged parser):"
python3 - "${RAW}" << 'PY'
import sys, os
sys.path.insert(0, "scripts/10_vmgc_crosscheck")
from vmgc_coverage import best_subject_per_query
f = sys.argv[1]
if not os.path.exists(f) or os.path.getsize(f) == 0:
    print("  RESULT: NO MetaVR v5 hit at e <= 1e-5 (Gardnerella phage absent from MetaVR by sequence)")
    sys.exit()
best = best_subject_per_query(f)
for q, m in best.items():
    species = "YES" if (m["weighted_identity"] >= 95 and m["af_shorter"] >= 85) else "no"
    print(f"  {q}")
    print(f"    best MetaVR subject : {m['sseqid']}")
    print(f"    weighted identity   : {m['weighted_identity']:.2f}%")
    print(f"    query coverage      : {m['q_af']:.2f}%   subject coverage: {m['s_af']:.2f}%")
    print(f"    AF of shorter seq   : {m['af_shorter']:.2f}%   (qlen {m['qlen']}, slen {m['slen']})")
    print(f"    species-level match : {species}  (>=95% ANI over >=85% AF of shorter)")
PY
echo "[$(date)] done. -> ${RAW}"
