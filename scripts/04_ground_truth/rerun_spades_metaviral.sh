#!/usr/bin/env bash
# =============================================================================
# rerun_spades_metaviral.sh — Re-assemble RCA reads with metaSPAdes
# =============================================================================
# Step 0: Re-assemble RCA-enriched virome reads using --meta -k 33,55,77
# (matching the reference paper approach). --metaviral was tried first but
# produced no circular contigs for this dataset.
#
# Usage: sbatch rerun_spades_metaviral.sh
# =============================================================================

#SBATCH --cpus-per-task=12
#SBATCH -t 4:00:00
#SBATCH --mem=64G
#SBATCH -J metaviral_spades
#SBATCH -o logs/spades_metaviral_%j.out
#SBATCH -e logs/spades_metaviral_%j.err

set -euo pipefail

# --- Configuration -----------------------------------------------------------
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
SAMPLE="UC_028_V2_RCA"
R1="${PROJ_DIR}/results/test_real/samtools/${SAMPLE}_host_removed_R1.fastq.gz"
R2="${PROJ_DIR}/results/test_real/samtools/${SAMPLE}_host_removed_R2.fastq.gz"
OUTDIR="${PROJ_DIR}/results/test_real/spades_metaviral/${SAMPLE}"
MIN_LENGTH=500  # RCA threshold per strategy §2.2 (lower than shotgun 1500 bp)
THREADS="${SLURM_CPUS_PER_TASK:-12}"
MEMORY=60  # GB — must stay below SBATCH --mem (64G)

mkdir -p logs

# --- Validate inputs ---------------------------------------------------------
for f in "${R1}" "${R2}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: Input file not found: ${f}" >&2
        exit 1
    fi
done

# --- Load modules ------------------------------------------------------------
module load SPAdes/4.2.0-GCC-13.3.0 2>/dev/null || true
module load Biopython/1.84-gfbf-2024a 2>/dev/null || true

# =============================================================================
# Run metaSPAdes with k = 33,55,77 (per reference paper)
# =============================================================================
echo "[$(date)] Running metaSPAdes on ${SAMPLE}..."
echo "  R1: ${R1}"
echo "  R2: ${R2}"
echo "  Output: ${OUTDIR}"
echo "  K-mers: 33,55,77"
echo "  Threads: ${THREADS}, Memory: ${MEMORY} GB"

# Clean output dir for fresh assembly (SPAdes refuses non-empty output)
rm -rf "${OUTDIR}"

spades.py \
    --meta \
    -k 33,55,77 \
    -1 "${R1}" \
    -2 "${R2}" \
    -o "${OUTDIR}" \
    -t "${THREADS}" \
    -m "${MEMORY}"

# Check SPAdes exit
if [[ ! -f "${OUTDIR}/contigs.fasta" ]]; then
    echo "ERROR: SPAdes did not produce contigs.fasta" >&2
    exit 1
fi

# =============================================================================
# Filter contigs by minimum length
# =============================================================================
echo "[$(date)] Filtering contigs >= ${MIN_LENGTH} bp..."

FILTERED="${OUTDIR}/contigs_filtered_${MIN_LENGTH}bp.fasta"

python3 -c "
from Bio import SeqIO
import sys

input_fa = '${OUTDIR}/contigs.fasta'
output_fa = '${FILTERED}'
min_len = ${MIN_LENGTH}

total = 0
kept = 0
with open(output_fa, 'w') as out:
    for rec in SeqIO.parse(input_fa, 'fasta'):
        total += 1
        if len(rec.seq) >= min_len:
            kept += 1
            SeqIO.write(rec, out, 'fasta')

print(f'  Total contigs: {total}')
print(f'  Kept (>= {min_len} bp): {kept}')
print(f'  Removed: {total - kept}')
"

# =============================================================================
# Assembly statistics comparison
# =============================================================================
echo ""
echo "[$(date)] Assembly statistics for NEW metaSPAdes assembly:"

python3 -c "
from Bio import SeqIO
import numpy as np

def assembly_stats(fasta_path, label):
    lengths = sorted([len(r.seq) for r in SeqIO.parse(fasta_path, 'fasta')], reverse=True)
    if not lengths:
        print(f'  {label}: No contigs')
        return
    total = sum(lengths)
    cum = 0
    n50 = 0
    for l in lengths:
        cum += l
        if cum >= total / 2:
            n50 = l
            break
    print(f'  {label}:')
    print(f'    Contigs:     {len(lengths)}')
    print(f'    Total bp:    {total:,}')
    print(f'    N50:         {n50:,}')
    print(f'    Longest:     {max(lengths):,}')
    print(f'    Shortest:    {min(lengths):,}')
    print(f'    Mean length: {np.mean(lengths):,.0f}')

assembly_stats('${OUTDIR}/contigs.fasta', 'metaSPAdes (raw)')
print()
assembly_stats('${FILTERED}', 'metaSPAdes (filtered >= ${MIN_LENGTH} bp)')

# Compare with old --meta assembly if it exists
import os
old_assembly = '${PROJ_DIR}/results/test_real/spades/${SAMPLE}_contigs.fasta'
if os.path.exists(old_assembly):
    print()
    assembly_stats(old_assembly, 'Old metaSPAdes assembly (for comparison)')
"

echo ""
echo "[$(date)] metaSPAdes re-assembly complete."
echo "  Raw contigs: ${OUTDIR}/contigs.fasta"
echo "  Filtered:    ${FILTERED}"
