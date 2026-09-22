#!/usr/bin/env bash
# =============================================================================
# run_all_rca_samples.sh — Submit RCA pipeline for all 13 samples
# =============================================================================
# RCA pipeline has two steps per sample:
#   Step 1: Host removal on RCA reads  (run_host_removal.sh with RCA sample IDs)
#   Step 2: RCA read mapping to shotgun contigs  (run_rca_read_mapping.sh)
#
# Step 2 requires BOTH host-removed RCA reads AND a shotgun assembly.
# When --submit-readmap is used, Step 2 is submitted with a SLURM dependency
# on the host-removal job (--dependency=afterok:JOBID).
#
# Naming convention:
#   Shotgun sample : UC028_V2
#   RCA sample     : UC_028_V2_RCA  (underscore after UC, _RCA suffix)
#   RCA raw reads  : UC-fq-merged/UC_028_V2_RCA_R{1,2}.fastq.gz
#   Host-removed   : results/test_real/samtools/UC_028_V2_RCA_host_removed_R{1,2}.fastq.gz
#   Shotgun contigs: results/test_real/spades/UC028_V2_contigs.fasta
#
# Usage:
#   bash scripts/04_ground_truth/run_all_rca_samples.sh                 # Host removal only (safe default)
#   bash scripts/04_ground_truth/run_all_rca_samples.sh --submit-readmap  # Host removal + read mapping (with SLURM dep)
#   bash scripts/04_ground_truth/run_all_rca_samples.sh --readmap-only    # Read mapping only (host removal must be done)
#   bash scripts/04_ground_truth/run_all_rca_samples.sh --status          # Print status only, submit nothing
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
READS_DIR="${PROJ_DIR}/UC-fq-merged"
SAMTOOLS_DIR="${PROJ_DIR}/results/test_real/samtools"
SPADES_DIR="${PROJ_DIR}/results/test_real/spades"
GT_DIR="${PROJ_DIR}/results/test_real/ground_truth"
HOST_REMOVAL_SCRIPT="${PROJ_DIR}/scripts/05_spike_in/run_host_removal.sh"
RCA_READMAP_SCRIPT="${SCRIPT_DIR}/run_rca_read_mapping.sh"

MODE="${1:---host-removal-only}"

# All 13 participant-timepoint pairs
# Format: SHOTGUN_SAMPLE:RCA_SAMPLE
SAMPLE_PAIRS=(
    "UC084_V2:UC_084_V2_RCA"
    "UC115_V2:UC_115_V2_RCA"
    "UC164_V2:UC_164_V2_RCA"
    "UC093_V2:UC_093_V2_RCA"
    "UC139_V2:UC_139_V2_RCA"
    "UC055_V1:UC_055_V1_RCA"
    "UC065_V2:UC_065_V2_RCA"
    "UC096_V2:UC_096_V2_RCA"
    "UC028_V2:UC_028_V2_RCA"
    "UC074_V2:UC_074_V2_RCA"
    "UC093_V3:UC_093_V3_RCA"
    "UC055_V2:UC_055_V2_RCA"
    "UC062_V2:UC_062_V2_RCA"
)

# =============================================================================
# Status check
# =============================================================================
check_status() {
    local shotgun="$1" rca="$2"

    local rca_raw_ok="NO"
    local rca_hr_ok="NO"
    local assembly_ok="NO"
    local e5_ok="NO"

    [[ -f "${READS_DIR}/${rca}_R1.fastq.gz" ]] && rca_raw_ok="YES"
    [[ -f "${SAMTOOLS_DIR}/${rca}_host_removed_R1.fastq.gz" ]] && \
        [[ -s "${SAMTOOLS_DIR}/${rca}_host_removed_R1.fastq.gz" ]] && rca_hr_ok="YES"
    [[ -f "${SPADES_DIR}/${shotgun}_contigs.fasta" ]] && assembly_ok="YES"
    [[ -f "${GT_DIR}/${shotgun}/evidence_5_rca_readlevel.tsv" ]] && e5_ok="YES"

    printf "  %-10s  %-14s  raw=%-3s  hr=%-3s  asm=%-3s  e5=%-3s\n" \
        "${shotgun}" "${rca}" "${rca_raw_ok}" "${rca_hr_ok}" "${assembly_ok}" "${e5_ok}"
}

echo "=========================================================================="
echo "RCA Pipeline — Batch Submission for 13 Samples"
echo "=========================================================================="
echo "  Mode: ${MODE}"
echo ""
echo "  Status (raw=RCA reads, hr=host-removed, asm=shotgun assembly, e5=RCA mapping):"
for pair in "${SAMPLE_PAIRS[@]}"; do
    IFS=':' read -r shotgun rca <<< "$pair"
    check_status "$shotgun" "$rca"
done
echo ""

if [[ "${MODE}" == "--status" ]]; then
    echo "Status check only — nothing submitted."
    exit 0
fi

# =============================================================================
# Step 1: Submit RCA host removal for samples that need it
# =============================================================================
if [[ "${MODE}" != "--readmap-only" ]]; then
    echo "=========================================================================="
    echo "Step 1: RCA Host Removal"
    echo "=========================================================================="

    n_submitted=0
    n_skipped=0
    declare -A HR_JOBIDS  # track SLURM job IDs for dependency chaining

    for pair in "${SAMPLE_PAIRS[@]}"; do
        IFS=':' read -r shotgun rca <<< "$pair"

        # Check if host-removed reads already exist
        hr_r1="${SAMTOOLS_DIR}/${rca}_host_removed_R1.fastq.gz"
        if [[ -f "${hr_r1}" ]] && [[ -s "${hr_r1}" ]]; then
            echo "  ${rca}: SKIP (host-removed reads exist)"
            HR_JOBIDS["${shotgun}"]="done"
            ((n_skipped++))
            continue
        fi

        # Check raw reads exist
        raw_r1="${READS_DIR}/${rca}_R1.fastq.gz"
        if [[ ! -f "${raw_r1}" ]]; then
            echo "  ${rca}: SKIP (raw reads not found: ${raw_r1})"
            ((n_skipped++))
            continue
        fi

        # Submit host removal
        jobid=$(sbatch --parsable "${HOST_REMOVAL_SCRIPT}" "${rca}")
        HR_JOBIDS["${shotgun}"]="${jobid}"
        echo "  ${rca}: SUBMITTED (job ${jobid})"
        ((n_submitted++))
    done

    echo ""
    echo "  Host removal: ${n_submitted} submitted, ${n_skipped} skipped"
    echo ""
fi

# =============================================================================
# Step 2: Submit RCA read mapping (requires host-removed RCA + shotgun assembly)
# =============================================================================
if [[ "${MODE}" == "--submit-readmap" ]] || [[ "${MODE}" == "--readmap-only" ]]; then
    echo "=========================================================================="
    echo "Step 2: RCA Read Mapping (E5)"
    echo "=========================================================================="

    n_submitted=0
    n_skipped=0
    n_blocked=0

    for pair in "${SAMPLE_PAIRS[@]}"; do
        IFS=':' read -r shotgun rca <<< "$pair"

        # Check E5 already exists
        e5_file="${GT_DIR}/${shotgun}/evidence_5_rca_readlevel.tsv"
        if [[ -f "${e5_file}" ]]; then
            echo "  ${shotgun}: SKIP (E5 already exists)"
            ((n_skipped++))
            continue
        fi

        # Check shotgun assembly exists
        contigs="${SPADES_DIR}/${shotgun}_contigs.fasta"
        if [[ ! -f "${contigs}" ]]; then
            echo "  ${shotgun}: BLOCKED (no shotgun assembly yet)"
            ((n_blocked++))
            continue
        fi

        # Determine dependency
        dep_flag=""
        if [[ "${MODE}" == "--submit-readmap" ]]; then
            hr_jobid="${HR_JOBIDS[${shotgun}]:-}"
            if [[ "${hr_jobid}" == "done" ]]; then
                dep_flag=""  # no dependency needed
            elif [[ -n "${hr_jobid}" ]]; then
                dep_flag="--dependency=afterok:${hr_jobid}"
            else
                # Check if host-removed reads exist already
                hr_r1="${SAMTOOLS_DIR}/${rca}_host_removed_R1.fastq.gz"
                if [[ ! -f "${hr_r1}" ]] || [[ ! -s "${hr_r1}" ]]; then
                    echo "  ${shotgun}: BLOCKED (no host-removed RCA reads, no pending job)"
                    ((n_blocked++))
                    continue
                fi
            fi
        else
            # --readmap-only: require host-removed reads to exist
            hr_r1="${SAMTOOLS_DIR}/${rca}_host_removed_R1.fastq.gz"
            if [[ ! -f "${hr_r1}" ]] || [[ ! -s "${hr_r1}" ]]; then
                echo "  ${shotgun}: BLOCKED (no host-removed RCA reads)"
                ((n_blocked++))
                continue
            fi
        fi

        # Submit read mapping
        if [[ -n "${dep_flag}" ]]; then
            jobid=$(sbatch --parsable ${dep_flag} "${RCA_READMAP_SCRIPT}" "${shotgun}" "${rca}")
            echo "  ${shotgun}: SUBMITTED (job ${jobid}, depends on ${hr_jobid})"
        else
            jobid=$(sbatch --parsable "${RCA_READMAP_SCRIPT}" "${shotgun}" "${rca}")
            echo "  ${shotgun}: SUBMITTED (job ${jobid})"
        fi
        ((n_submitted++))
    done

    echo ""
    echo "  Read mapping: ${n_submitted} submitted, ${n_skipped} skipped, ${n_blocked} blocked"
    echo ""
fi

# =============================================================================
# Summary
# =============================================================================
echo "=========================================================================="
echo "Done."
echo "=========================================================================="
echo ""
echo "  Monitor jobs:   squeue -u \$USER"
echo "  Check status:   bash $0 --status"
echo ""
echo "  After completion, rebuild ground truth:"
echo "    for s in UC084_V2 UC115_V2 UC164_V2 UC093_V2 UC139_V2 \\"
echo "             UC055_V1 UC065_V2 UC096_V2 UC028_V2 \\"
echo "             UC074_V2 UC093_V3 UC055_V2 UC062_V2; do"
echo "      python scripts/04_ground_truth/build_ground_truth.py \\"
echo "          --evidence-dir results/test_real/ground_truth/\$s/ \\"
echo "          --contigs results/test_real/spades/\${s}_contigs.fasta \\"
echo "          --output results/test_real/ground_truth/\$s/ground_truth.tsv \\"
echo "          --e5-mode read"
echo "    done"
echo ""
echo "  Then run RCA concordance (§7.5.2):"
echo "    python scripts/07_evaluation/rca_concordance.py \\"
echo "        --samples UC084_V2 UC115_V2 ... \\"
echo "        --evidence-root results/test_real/ground_truth/ \\"
echo "        --results-root results/test_real/full_run/ \\"
echo "        --e1a-root results/test_real/ground_truth/ \\"
echo "        --output-dir results/test_real/rca_concordance/"
