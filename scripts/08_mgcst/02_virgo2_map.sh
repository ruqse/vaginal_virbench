#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH -t 4:00:00
#SBATCH --mem=128G
#SBATCH -J mgcst_map
#SBATCH -o logs/mgcst_map_%A_%a.out
#SBATCH -e logs/mgcst_map_%A_%a.err

# SLURM array job: map each sample to VIRGO2 gene catalog.
# Array index selects sample from 01_array_samples.txt.
# Reuses VIRGO2 installed at the VISTA project location.

set -euo pipefail

PROJ="${PROJ:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
VIRGO2="${VIRGO2:?set VIRGO2 (see config/paths.example.sh)}"
MGCST_RESULTS="${MGCST_RESULTS:-${PROJ}/results/mgcst}"
CONCATDIR="${MGCST_RESULTS}/concat_reads"
OUTDIR="${MGCST_RESULTS}/virgo2_maps"
ARRAY_LIST="${ARRAY_LIST:-${PROJ}/scripts/08_mgcst/01_array_samples.txt}"

module load SciPy-bundle/2024.05-gfbf-2024a
module load Bowtie2/2.5.4-GCC-13.3.0
module load SAMtools/1.22-GCC-13.3.0

mkdir -p "${OUTDIR}"

# Get sample name from array index (1-based line number)
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ARRAY_LIST}")

if [[ -z "${SAMPLE}" ]]; then
    echo "ERROR: No sample found at array index ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

READS="${CONCATDIR}/${SAMPLE}.fq.gz"

if [[ ! -f "${READS}" ]]; then
    echo "ERROR: Read file not found: ${READS}"
    exit 1
fi

# Skip if already done
if [[ -s "${OUTDIR}/${SAMPLE}.out" ]]; then
    echo "[$(date)] ${SAMPLE} SKIPPED -- output already exists"
    exit 0
fi

echo "=== VIRGO2 mapping ==="
echo "Sample:     ${SAMPLE}"
echo "Reads:      ${READS} ($(du -h ${READS} | cut -f1))"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Threads:    ${SLURM_CPUS_PER_TASK}"
echo "Start:      $(date)"

python3 "${VIRGO2}" map \
    -r "${READS}" \
    -p ${SLURM_CPUS_PER_TASK} \
    -o "${OUTDIR}/${SAMPLE}"

echo "=== Mapping complete ==="
if [[ -f "${OUTDIR}/${SAMPLE}.out" ]]; then
    echo "Output: ${OUTDIR}/${SAMPLE}.out"
    echo "Genes detected: $(wc -l < ${OUTDIR}/${SAMPLE}.out)"
else
    echo "ERROR: Output file not produced"
    exit 1
fi
echo "End: $(date)"
