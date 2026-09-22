#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH -t 4:00:00
#SBATCH --mem=8G
#SBATCH -J mgcst_prepare
#SBATCH -o logs/mgcst_prepare_%j.out
#SBATCH -e logs/mgcst_prepare_%j.err

# Concatenate paired-end R1+R2 reads into a single file per sample for VIRGO2.
# VIRGO2 only supports single-end input (-U flag in bowtie2).
# Concatenation is standard for gene catalog mapping since each read maps independently.

set -euo pipefail

PROJ="${PROJ:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
READDIR="${READDIR:-${PROJ}/UC-fq-merged}"
# READS_SUFFIX_R1/R2 default to UChoose convention "_R{1,2}.fastq.gz"; override
# to "_host_removed_R{1,2}.fastq.gz" for expansion-cohort cleaned reads.
READS_SUFFIX_R1="${READS_SUFFIX_R1:-_R1.fastq.gz}"
READS_SUFFIX_R2="${READS_SUFFIX_R2:-_R2.fastq.gz}"
MGCST_RESULTS="${MGCST_RESULTS:-${PROJ}/results/mgcst}"
OUTDIR="${MGCST_RESULTS}/concat_reads"
SAMPLE_LIST="${SAMPLE_LIST:-${PROJ}/scripts/08_mgcst/00_sample_list.txt}"
ARRAY_LIST="${ARRAY_LIST:-${PROJ}/scripts/08_mgcst/01_array_samples.txt}"

mkdir -p "${OUTDIR}"

echo "=== Concatenating paired-end reads for VIRGO2 ==="
echo "Start: $(date)"

PROCESSED=0
SKIPPED=0
FAILED=0

while IFS= read -r sample; do
    R1="${READDIR}/${sample}${READS_SUFFIX_R1}"
    R2="${READDIR}/${sample}${READS_SUFFIX_R2}"
    OUTFILE="${OUTDIR}/${sample}.fq.gz"

    # Skip if already done
    if [[ -s "${OUTFILE}" ]]; then
        echo "[$(date)] ${sample} SKIPPED -- already exists ($(du -h ${OUTFILE} | cut -f1))"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    if [[ -f "${R1}" ]] && [[ -f "${R2}" ]]; then
        cat "${R1}" "${R2}" > "${OUTFILE}"
        echo "[$(date)] ${sample} OK ($(du -h ${OUTFILE} | cut -f1))"
        PROCESSED=$((PROCESSED + 1))
    elif [[ -f "${R1}" ]]; then
        cp "${R1}" "${OUTFILE}"
        echo "[$(date)] ${sample} WARN: R1 only"
        PROCESSED=$((PROCESSED + 1))
    else
        echo "[$(date)] ${sample} ERROR: No reads found"
        FAILED=$((FAILED + 1))
    fi
done < "${SAMPLE_LIST}"

# Write array sample list for the mapping step (only samples with reads).
# ARRAY_LIST is now configurable via env var (set in the configuration block above).
ls -1 "${OUTDIR}"/*.fq.gz | xargs -I{} basename {} .fq.gz > "${ARRAY_LIST}"

echo ""
echo "=== Summary ==="
echo "Processed: ${PROCESSED}"
echo "Skipped:   ${SKIPPED}"
echo "Failed:    ${FAILED}"
echo "Array list: ${ARRAY_LIST} ($(wc -l < ${ARRAY_LIST}) samples)"
echo "End: $(date)"
