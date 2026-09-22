#!/usr/bin/env bash
# =============================================================================
# run_standalone_assembly.sh — Standalone metaSPAdes Assembly for Background Samples
# =============================================================================
# Assembles host-removed reads for a single sample using metaSPAdes.
# Produces standalone assemblies needed for:
#   1. Spike-in homology check (check_spikein_homology.sh)
#   2. Secondary benchmark (ground truth construction + 14-tool evaluation)
#
# Uses scratch ($SNIC_TMP) for SPAdes working files, copies key outputs
# to persistent storage.
#
# Usage:
#   sbatch run_standalone_assembly.sh UC115_V2   # CST-I background
#   sbatch run_standalone_assembly.sh UC093_V3   # CST-IV background
#   sbatch run_standalone_assembly.sh UC028_V2   # Any sample with host-removed reads
# =============================================================================

#SBATCH --cpus-per-task=12
#SBATCH -t 1-00:00:00
#SBATCH --mem=72G
#SBATCH -J standalone_asm
#SBATCH -o logs/standalone_asm_%j.out
#SBATCH -e logs/standalone_asm_%j.err

set -euo pipefail

# --- Sample ID from argument ------------------------------------------------
SAMPLE="${1:?Usage: sbatch $0 SAMPLE_ID (e.g. UC115_V2, UC093_V3)}"

# --- Configuration -----------------------------------------------------------
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
R1="${PROJ_DIR}/results/test_real/samtools/${SAMPLE}_host_removed_R1.fastq.gz"
R2="${PROJ_DIR}/results/test_real/samtools/${SAMPLE}_host_removed_R2.fastq.gz"
OUTDIR="${PROJ_DIR}/results/test_real/spades"
ASM_WORK="${SNIC_TMP:-${TMPDIR:-/tmp}}/spades_${SAMPLE}_$$"
THREADS="${SLURM_CPUS_PER_TASK:-12}"
MEMORY=60  # GB for SPAdes (leave headroom from 72G SLURM allocation)

mkdir -p "${OUTDIR}" "${ASM_WORK}" logs

# --- Validate inputs ---------------------------------------------------------
if [[ ! -f "${R1}" ]]; then
    echo "ERROR: R1 not found: ${R1}" >&2
    echo "  Run run_host_removal.sh ${SAMPLE} first." >&2
    exit 1
fi
if [[ ! -f "${R2}" ]]; then
    echo "ERROR: R2 not found: ${R2}" >&2
    exit 1
fi

# Skip if assembly already exists
if [[ -f "${OUTDIR}/${SAMPLE}_contigs.fasta" ]]; then
    N_EXISTING=$(grep -c '^>' "${OUTDIR}/${SAMPLE}_contigs.fasta")
    echo "Assembly already exists: ${OUTDIR}/${SAMPLE}_contigs.fasta (${N_EXISTING} contigs)"
    echo "Delete it first to re-assemble."
    exit 0
fi

# --- Load SPAdes module ------------------------------------------------------
module load SPAdes/4.2.0-GCC-13.3.0

# --- Assembly info -----------------------------------------------------------
R1_SIZE=$(du -sh "${R1}" | cut -f1)
R2_SIZE=$(du -sh "${R2}" | cut -f1)

echo "================================================================"
echo "Standalone metaSPAdes assembly: ${SAMPLE}"
echo "  R1: ${R1} (${R1_SIZE})"
echo "  R2: ${R2} (${R2_SIZE})"
echo "  Work dir: ${ASM_WORK}"
echo "  Output: ${OUTDIR}/${SAMPLE}_contigs.fasta"
echo "  Threads: ${THREADS}, Memory: ${MEMORY}G"
echo "  SPAdes: $(spades.py --version 2>&1 || echo 'loaded')"
echo "  Started: $(date)"
echo "================================================================"
echo ""

# --- Run metaSPAdes ----------------------------------------------------------
echo "[$(date)] Running metaSPAdes..."
spades.py --meta -k 33,55,77 \
    -1 "${R1}" \
    -2 "${R2}" \
    -o "${ASM_WORK}" \
    -t "${THREADS}" \
    -m "${MEMORY}" \
    2>&1 | tail -20

echo ""
echo "[$(date)] Assembly complete."

# --- Validate output ---------------------------------------------------------
if [[ ! -f "${ASM_WORK}/contigs.fasta" ]]; then
    echo "ERROR: SPAdes did not produce contigs.fasta" >&2
    echo "  Check ${ASM_WORK}/spades.log for errors." >&2
    # Copy log even on failure
    cp "${ASM_WORK}/spades.log" "${OUTDIR}/${SAMPLE}_spades.log" 2>/dev/null || true
    exit 1
fi

# --- Copy key outputs to persistent storage ----------------------------------
echo "Copying outputs to ${OUTDIR}/"
cp "${ASM_WORK}/contigs.fasta"  "${OUTDIR}/${SAMPLE}_contigs.fasta"
cp "${ASM_WORK}/scaffolds.fasta" "${OUTDIR}/${SAMPLE}_scaffolds.fasta"
cp "${ASM_WORK}/assembly_graph_with_scaffolds.gfa" "${OUTDIR}/${SAMPLE}_graph.gfa"
cp "${ASM_WORK}/spades.log" "${OUTDIR}/${SAMPLE}_spades.log"

# --- Summary -----------------------------------------------------------------
N_CONTIGS=$(grep -c '^>' "${OUTDIR}/${SAMPLE}_contigs.fasta")
N_LONG=$(python3 -c "
count = 0
with open('${OUTDIR}/${SAMPLE}_contigs.fasta') as f:
    seq_len = 0
    for line in f:
        if line.startswith('>'):
            if seq_len >= 1500: count += 1
            seq_len = 0
        else:
            seq_len += len(line.strip())
    if seq_len >= 1500: count += 1
print(count)
")

echo ""
echo "================================================================"
echo "Assembly complete: ${SAMPLE}"
echo "  Total contigs: ${N_CONTIGS}"
echo "  Contigs >= 1500 bp: ${N_LONG}"
echo "  Output: ${OUTDIR}/${SAMPLE}_contigs.fasta"
echo "  Finished: $(date)"
echo ""
echo "Next steps:"
echo "  1. Homology check: sbatch scripts/05_spike_in/check_spikein_homology.sh ${SAMPLE}"
echo "  2. Ground truth:   bash scripts/04_ground_truth/run_all_samples.sh"
echo "  3. Tool execution: sbatch scripts/06_tool_execution/run_secondary_cpu_tools.sh ${SAMPLE}"
echo "================================================================"

# --- Cleanup scratch ---------------------------------------------------------
rm -rf "${ASM_WORK}"
