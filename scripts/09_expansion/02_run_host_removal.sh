#!/usr/bin/env bash
# =============================================================================
# 02_run_host_removal.sh — Per-sample fastp + Bowtie2 hg19 host removal
# =============================================================================
# SLURM array wrapper around the existing 05_spike_in/run_host_removal.sh.
# Reads sample IDs from data/expansion_cohorts/array_samples.txt; processes
# one sample per array task with all heavy intermediates (SAM, sorted BAM)
# on $SNIC_TMP so /proj never sees a SAM file.
#
# Output: results/expansion/clean_reads/<run>_host_removed_R{1,2}.fastq.gz
#
# Usage:
#   sbatch --array=1-138 scripts/09_expansion/02_run_host_removal.sh
#
# The existing run_host_removal.sh handles fastp + FastQC + Bowtie2 + samtools.
# We override READS_DIR, RESULTS_DIR, READS_SUFFIX, SCRATCH_DIR via env vars.
# =============================================================================

#SBATCH --cpus-per-task=12
#SBATCH -t 12:00:00
#SBATCH --mem=36G
#SBATCH -J expand_hostrm
#SBATCH -o logs/expand_hostrm_%A_%a.out
#SBATCH -e logs/expand_hostrm_%A_%a.err

set -euo pipefail

PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
EXPANSION_DIR="${PROJ_DIR}/data/expansion_cohorts"
ARRAY_TXT="${EXPANSION_DIR}/array_samples.txt"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: must be run as a SLURM array job (sbatch --array=...)" >&2
    exit 1
fi

if [[ ! -f "${ARRAY_TXT}" ]]; then
    echo "ERROR: array list missing: ${ARRAY_TXT}" >&2
    echo "Run 01_symlink_reads.sh first." >&2
    exit 1
fi

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ARRAY_TXT}")
if [[ -z "${SAMPLE}" ]]; then
    echo "ERROR: no sample at array index ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

# Skip if already done — idempotent re-run support
OUT_R1="${PROJ_DIR}/results/expansion/clean_reads/${SAMPLE}_host_removed_R1.fastq.gz"
OUT_R2="${PROJ_DIR}/results/expansion/clean_reads/${SAMPLE}_host_removed_R2.fastq.gz"
if [[ -s "${OUT_R1}" ]] && [[ -s "${OUT_R2}" ]]; then
    echo "[$(date)] SKIP ${SAMPLE} (already host-removed)"
    exit 0
fi

# Per-task scratch (auto-cleaned by SLURM at job end)
SCRATCH_BASE="${SNIC_TMP:-${TMPDIR:-/tmp}}"
mkdir -p "${SCRATCH_BASE}"

# Stage the input symlinks into a per-task working directory under scratch
# so the existing run_host_removal.sh (which writes alongside its READS_DIR
# when running fastp logs) cannot accidentally try to write back to /proj.
TASK_READS="${SCRATCH_BASE}/reads_${SAMPLE}"
TASK_RESULTS="${SCRATCH_BASE}/results_${SAMPLE}"
mkdir -p "${TASK_READS}" "${TASK_RESULTS}"

# Resolve symlinks → real source files (Bowtie2 can read through symlinks but
# we copy the symlink for resilience). Symlinks point at ENA fastq.gz.
ln -sf "${EXPANSION_DIR}/reads/symlinks/${SAMPLE}_1.fastq.gz" "${TASK_READS}/${SAMPLE}_1.fastq.gz"
ln -sf "${EXPANSION_DIR}/reads/symlinks/${SAMPLE}_2.fastq.gz" "${TASK_READS}/${SAMPLE}_2.fastq.gz"

echo "[$(date)] Host removal start: sample=${SAMPLE} array_task=${SLURM_ARRAY_TASK_ID}"
echo "  reads dir : ${TASK_READS}"
echo "  results   : ${TASK_RESULTS}"
echo "  scratch   : ${SCRATCH_BASE}"

# Invoke the existing script. We do NOT use sbatch from inside this array task —
# we run the script's body inline. The existing script auto-detects positional
# args as sample names, so we just call it with our SAMPLE.
export PROJ_DIR
export READS_DIR="${TASK_READS}"
export RESULTS_DIR="${TASK_RESULTS}"
export READS_SUFFIX_R1="_1.fastq.gz"
export READS_SUFFIX_R2="_2.fastq.gz"
export SCRATCH_DIR="${SCRATCH_BASE}/sam"

bash "${PROJ_DIR}/scripts/05_spike_in/run_host_removal.sh" "${SAMPLE}"

# Move only the cleaned FASTQs back to /proj
FINAL_DIR="${PROJ_DIR}/results/expansion/clean_reads"
mkdir -p "${FINAL_DIR}"

TASK_OUT_R1="${TASK_RESULTS}/samtools/${SAMPLE}_host_removed_R1.fastq.gz"
TASK_OUT_R2="${TASK_RESULTS}/samtools/${SAMPLE}_host_removed_R2.fastq.gz"

if [[ -s "${TASK_OUT_R1}" ]] && [[ -s "${TASK_OUT_R2}" ]]; then
    cp -p "${TASK_OUT_R1}" "${OUT_R1}"
    cp -p "${TASK_OUT_R2}" "${OUT_R2}"
    echo ""
    echo "[$(date)] Copied to /proj:"
    echo "  ${OUT_R1} ($(du -h "${OUT_R1}" | cut -f1))"
    echo "  ${OUT_R2} ($(du -h "${OUT_R2}" | cut -f1))"
else
    echo "ERROR: host-removed outputs missing or empty" >&2
    ls -la "${TASK_RESULTS}/samtools/" >&2 || true
    exit 1
fi

# Verify no SAM survived to /proj
if find "${PROJ_DIR}/results/expansion" -name '*.sam' 2>/dev/null | grep -q .; then
    echo "ERROR: SAM file leaked to /proj. Check SCRATCH_DIR routing." >&2
    exit 1
fi

# SLURM auto-cleans $SNIC_TMP, but be explicit
rm -rf "${TASK_READS}" "${TASK_RESULTS}"

echo "[$(date)] Done ${SAMPLE}"
