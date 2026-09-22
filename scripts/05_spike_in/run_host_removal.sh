#!/usr/bin/env bash
# =============================================================================
# run_host_removal.sh — QC + Host Removal for Track B Background Samples
# =============================================================================
# Bash equivalent of the Nextflow pipeline: fastp -> FastQC -> host removal.
# Mirrors modules/local/{fastp,fastqc,host_removal,multiqc}/main.nf exactly.
#
# Host removal uses on-disk SAM (same as the Nextflow BOWTIE2_HOST_ALIGN module)
# to avoid SIGPIPE issues with piped workflows under set -euo pipefail.
#
# Steps (per sample):
#   1. fastp — QC trim + adapter removal (PE, --very-sensitive-local)
#   2. FastQC — quality assessment on trimmed reads
#   3. Bowtie2 — align to hg19 → SAM
#   4. samtools view -f 13 -F 256 — keep both-unmapped pairs
#   5. samtools sort -n — name-sort for FASTQ extraction
#   6. samtools fastq — BAM to host-removed FASTQ pairs
#   7. MultiQC — aggregate all reports
#
# Output:
#   results/test_real/fastp/{SAMPLE}.fastp.{json,html}
#   results/test_real/fastqc/{SAMPLE}_{1,2}_fastqc.{html,zip}
#   results/test_real/samtools/{SAMPLE}_host_removed_R{1,2}.fastq.gz
#   results/test_real/multiqc/multiqc_report.html
#
# Usage:
#   sbatch run_host_removal.sh                    # Both UC115_V2 + UC093_V3
#   sbatch run_host_removal.sh UC115_V2           # Single sample
# =============================================================================

#SBATCH --cpus-per-task=12
#SBATCH -t 08:00:00
#SBATCH --mem=36G
#SBATCH -J host_removal
#SBATCH -o logs/host_removal_%j.out
#SBATCH -e logs/host_removal_%j.err

set -euo pipefail

# --- Configuration -----------------------------------------------------------
# All paths overridable via env vars so the script can be reused for the
# expansion-cohort pipeline (scripts/09_expansion/) without breaking the
# original UC115_V2/UC093_V3 invocation.
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
READS_DIR="${READS_DIR:-${PROJ_DIR}/UC-fq-merged}"
RESULTS_DIR="${RESULTS_DIR:-${PROJ_DIR}/results/test_real}"
# READS_SUFFIX_R1/R2 default to "_R{1,2}.fastq.gz" (UChoose convention).
# Expansion cohorts use "_{1,2}.fastq.gz" — wrapper overrides these.
READS_SUFFIX_R1="${READS_SUFFIX_R1:-_R1.fastq.gz}"
READS_SUFFIX_R2="${READS_SUFFIX_R2:-_R2.fastq.gz}"
BOWTIE2_INDEX="${BOWTIE2_INDEX:-/sw/data/reference/Homo_sapiens/hg19/program_files/bowtie2/concat}"
THREADS="${SLURM_CPUS_PER_TASK:-12}"

# Samples to process (positional argument or default to both)
if [[ $# -ge 1 ]]; then
    SAMPLES=("$@")
else
    SAMPLES=("UC115_V2" "UC093_V3")
fi

# Output directories (match Nextflow publishDir layout)
FASTP_DIR="${RESULTS_DIR}/fastp"
FASTQC_DIR="${RESULTS_DIR}/fastqc"
SAMTOOLS_DIR="${RESULTS_DIR}/samtools"
MULTIQC_DIR="${RESULTS_DIR}/multiqc"

# SAM/BAM intermediates: by default land in SAMTOOLS_DIR (UChoose pipeline).
# Expansion pipeline overrides SCRATCH_DIR=$SNIC_TMP to keep heavy
# intermediates off /proj. SAM is auto-deleted after BAM is written either way.
SCRATCH_DIR="${SCRATCH_DIR:-${SAMTOOLS_DIR}}"

mkdir -p "${FASTP_DIR}" "${FASTQC_DIR}" "${SAMTOOLS_DIR}" "${MULTIQC_DIR}" "${SCRATCH_DIR}" logs

# --- Load modules (the HPC cluster: no bioinfo-tools needed) ---------------------------
module load fastp/1.0.1
module load FastQC/0.12.1
module load Bowtie2/2.5.4
module load SAMtools/1.21

echo "[$(date)] QC + Host removal pipeline"
echo "  fastp:   $(fastp --version 2>&1 || echo loaded)"
echo "  FastQC:  $(fastqc --version 2>&1 | head -1 || echo loaded)"
echo "  Bowtie2: $(bowtie2 --version 2>&1 | head -1 || echo loaded)"
echo "  samtools: $(samtools --version 2>&1 | head -1 || echo loaded)"
echo "  Bowtie2 index: ${BOWTIE2_INDEX}"
echo "  Threads:       ${THREADS}"
echo "  Samples:       ${SAMPLES[*]}"
echo ""

# =============================================================================
# Process each sample
# =============================================================================
for SAMPLE in "${SAMPLES[@]}"; do
    echo "================================================================"
    echo "[$(date)] Processing: ${SAMPLE}"
    echo "================================================================"

    # Input reads (suffix configurable via READS_SUFFIX_R1/R2 env vars)
    R1="${READS_DIR}/${SAMPLE}${READS_SUFFIX_R1}"
    R2="${READS_DIR}/${SAMPLE}${READS_SUFFIX_R2}"

    # Validate inputs
    if [[ ! -f "${R1}" ]] || [[ ! -f "${R2}" ]]; then
        echo "  ERROR: Input reads not found:" >&2
        echo "    ${R1}" >&2
        echo "    ${R2}" >&2
        continue
    fi

    # Final output files
    OUT_R1="${SAMTOOLS_DIR}/${SAMPLE}_host_removed_R1.fastq.gz"
    OUT_R2="${SAMTOOLS_DIR}/${SAMPLE}_host_removed_R2.fastq.gz"

    # Skip if already done
    if [[ -f "${OUT_R1}" ]] && [[ -s "${OUT_R1}" ]] && \
       [[ -f "${OUT_R2}" ]] && [[ -s "${OUT_R2}" ]]; then
        echo "  SKIP: Host-removed reads already exist"
        continue
    fi

    # ── Step 1: fastp — QC trimming ──────────────────────────────────
    # Mirrors: modules/local/fastp/main.nf + conf/modules.config args
    TRIM_R1="${FASTP_DIR}/${SAMPLE}_1.fastp.fastq.gz"
    TRIM_R2="${FASTP_DIR}/${SAMPLE}_2.fastp.fastq.gz"

    if [[ -f "${TRIM_R1}" ]] && [[ -s "${TRIM_R1}" ]]; then
        echo "  [1/6] fastp: SKIP (trimmed reads exist)"
    else
        echo "  [1/6] fastp: QC trimming..."
        fastp \
            --in1 "${R1}" \
            --in2 "${R2}" \
            --out1 "${TRIM_R1}" \
            --out2 "${TRIM_R2}" \
            --json "${FASTP_DIR}/${SAMPLE}.fastp.json" \
            --html "${FASTP_DIR}/${SAMPLE}.fastp.html" \
            --thread "${THREADS}" \
            --qualified_quality_phred 20 \
            --length_required 50 \
            --cut_front \
            --cut_tail \
            --cut_mean_quality 20 \
            --detect_adapter_for_pe \
            2> "${FASTP_DIR}/${SAMPLE}.fastp.log"
        echo "    Done: $(grep 'reads passed filter' "${FASTP_DIR}/${SAMPLE}.fastp.log" | tail -1 || echo 'check log')"
    fi

    # ── Step 2: FastQC — quality assessment on trimmed reads ─────────
    # Mirrors: modules/local/fastqc/main.nf
    if [[ -f "${FASTQC_DIR}/${SAMPLE}_1.fastp_fastqc.html" ]]; then
        echo "  [2/6] FastQC: SKIP (reports exist)"
    else
        echo "  [2/6] FastQC: quality assessment..."
        fastqc \
            --threads "${THREADS}" \
            --outdir "${FASTQC_DIR}" \
            "${TRIM_R1}" "${TRIM_R2}"
    fi

    # ── Step 3: Bowtie2 — align to hg19 (SAM to disk) ───────────────
    # Mirrors: BOWTIE2_HOST_ALIGN with ext.args = '--very-sensitive-local'
    # SAM goes to SCRATCH_DIR (== SAMTOOLS_DIR by default; == $SNIC_TMP under expansion).
    SAM_FILE="${SCRATCH_DIR}/${SAMPLE}.sam"

    if [[ -f "${SAM_FILE}" ]] && [[ -s "${SAM_FILE}" ]]; then
        echo "  [3/6] Bowtie2: SKIP (SAM exists)"
    else
        echo "  [3/6] Bowtie2: aligning to hg19..."
        bowtie2 \
            --very-sensitive-local \
            -p "${THREADS}" \
            -x "${BOWTIE2_INDEX}" \
            -1 "${TRIM_R1}" \
            -2 "${TRIM_R2}" \
            -S "${SAM_FILE}" \
            2> "${SAMTOOLS_DIR}/${SAMPLE}_bowtie2.log"

        ALIGN_RATE=$(grep "overall alignment rate" "${SAMTOOLS_DIR}/${SAMPLE}_bowtie2.log" | head -1)
        echo "    ${ALIGN_RATE}"
    fi

    # ── Step 4+5: samtools filter unmapped pairs + name-sort ─────────
    # Mirrors: SAMTOOLS_FILTER_UNMAPPED
    # -f 13: paired (1) + both unmapped (4+8=12) = 13
    # -F 256: exclude secondary alignments
    BAM_UNMAPPED="${SCRATCH_DIR}/${SAMPLE}_unmapped_sorted.bam"

    if [[ -f "${BAM_UNMAPPED}" ]] && [[ -s "${BAM_UNMAPPED}" ]]; then
        echo "  [4/6] samtools filter: SKIP (BAM exists)"
    else
        echo "  [4/6] samtools filter: extracting unmapped pairs (-f 13 -F 256)..."
        samtools view \
            -b -f 13 -F 256 \
            -@ "${THREADS}" \
            "${SAM_FILE}" \
        | samtools sort \
            -n \
            -@ "${THREADS}" \
            -o "${BAM_UNMAPPED}"

        echo "    Unmapped BAM: $(du -h "${BAM_UNMAPPED}" | cut -f1)"
    fi

    # Clean up SAM (large file, ~50-100 GB)
    if [[ -f "${BAM_UNMAPPED}" ]] && [[ -s "${BAM_UNMAPPED}" ]]; then
        rm -f "${SAM_FILE}"
        echo "  [5/6] Cleaned up SAM file"
    fi

    # ── Step 6: samtools fastq — BAM to FASTQ ───────────────────────
    # Mirrors: SAMTOOLS_BAM_TO_FASTQ
    echo "  [6/6] samtools fastq: converting to FASTQ..."
    samtools fastq \
        -@ "${THREADS}" \
        -1 "${OUT_R1}" \
        -2 "${OUT_R2}" \
        -s /dev/null \
        -n \
        "${BAM_UNMAPPED}"

    # Verify
    if [[ -f "${OUT_R1}" ]] && [[ -s "${OUT_R1}" ]]; then
        # Safe read counting without SIGPIPE
        OUT_READS=$(samtools view -c "${BAM_UNMAPPED}" 2>/dev/null | awk '{print int($1/2)}')
        echo ""
        echo "  Result: ~${OUT_READS} host-removed read pairs"
        echo "    R1: ${OUT_R1} ($(du -h "${OUT_R1}" | cut -f1))"
        echo "    R2: ${OUT_R2} ($(du -h "${OUT_R2}" | cut -f1))"
    else
        echo "  ERROR: Output files empty or missing!" >&2
    fi

    echo ""
done

# =============================================================================
# MultiQC — aggregate all reports
# =============================================================================
echo "[$(date)] Running MultiQC..."
module load MultiQC/1.28

multiqc \
    --force \
    --verbose \
    --outdir "${MULTIQC_DIR}" \
    "${FASTP_DIR}" "${FASTQC_DIR}" "${SAMTOOLS_DIR}" \
    2>&1 | tail -5

echo ""

# =============================================================================
# Summary
# =============================================================================
echo "[$(date)] Pipeline complete."
echo ""
echo "  Outputs:"
echo "    fastp:         ${FASTP_DIR}/"
echo "    FastQC:        ${FASTQC_DIR}/"
echo "    Host-removed:  ${SAMTOOLS_DIR}/"
echo "    MultiQC:       ${MULTIQC_DIR}/multiqc_report.html"
echo ""
for SAMPLE in "${SAMPLES[@]}"; do
    r1="${SAMTOOLS_DIR}/${SAMPLE}_host_removed_R1.fastq.gz"
    if [[ -f "${r1}" ]]; then
        echo "    ${SAMPLE}: $(du -h "${r1}" | cut -f1) (R1)"
    else
        echo "    ${SAMPLE}: MISSING"
    fi
done
echo ""
echo "  Next step: sbatch mix_and_assemble.sh all"
