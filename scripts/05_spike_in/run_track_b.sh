#!/usr/bin/env bash
# =============================================================================
# run_track_b.sh — Track B Orchestrator: Simulate + Host-Remove + Mix + Assemble
# =============================================================================
# SLURM dependency chain for the full Track B (assembly-based) benchmark:
#
#   Job 1: simulate_reads.sh     — Simulate reads for all 14 viral genomes
#   Job 2: run_host_removal.sh   — Host removal for UC115_V2 + UC093_V3
#   Job 3: mix_and_assemble.sh   — Mix + co-assemble (depends on jobs 1+2)
#   Job 4: label_contigs.sh      — Label contigs with ground truth (depends on job 3)
#
# Jobs 1 and 2 run in parallel. Job 3 waits for both.
#
# Usage: bash run_track_b.sh              # Submit full chain
#        bash run_track_b.sh --dry-run    # Show commands without submitting
# =============================================================================

set -euo pipefail

PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
SCRIPT_DIR="${PROJ_DIR}/scripts/05_spike_in"
SPIKE_IN_DIR="${PROJ_DIR}/data/spike_in"

DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=true
    echo "=== DRY RUN MODE ==="
    echo ""
fi

cd "${PROJ_DIR}"
mkdir -p scripts/05_spike_in/logs

echo "================================================================"
echo "Track B: Assembly-Based Spike-In Benchmark"
echo "================================================================"
echo ""
echo "  Backgrounds: UC115_V2 (CST-I), UC093_V3 (CST-IV)"
echo "  Viral panel: 14 genomes (10 original + 4 MetaVR phages)"
echo "  Coverages:   0.1x, 0.5x, 1x, 5x, 10x, 50x"
echo "  Assembler:   metaSPAdes v4.2.0 (--meta -k 33,55,77)"
echo ""
echo "  Total assemblies: 2 backgrounds x 6 coverages + 6 controls = 18"
echo ""

# =============================================================================
# Job 1: Simulate reads for all 14 viral genomes
# =============================================================================
echo "--- Job 1: Read simulation (all 14 viral genomes) ---"

if [[ "${DRY_RUN}" == true ]]; then
    echo "  sbatch ${SCRIPT_DIR}/simulate_reads.sh"
    JOB1="DRY_JOB1"
else
    JOB1=$(sbatch --parsable "${SCRIPT_DIR}/simulate_reads.sh")
    echo "  Submitted: ${JOB1}"
fi

# =============================================================================
# Job 2: Host removal for UC115_V2 + UC093_V3
# =============================================================================
echo ""
echo "--- Job 2: Host removal (UC115_V2 + UC093_V3) ---"

if [[ "${DRY_RUN}" == true ]]; then
    echo "  sbatch ${SCRIPT_DIR}/run_host_removal.sh"
    JOB2="DRY_JOB2"
else
    JOB2=$(sbatch --parsable "${SCRIPT_DIR}/run_host_removal.sh")
    echo "  Submitted: ${JOB2}"
fi

# =============================================================================
# Job 3: Mix + assemble (depends on jobs 1 + 2)
# =============================================================================
echo ""
echo "--- Job 3: Mix & co-assemble (depends on jobs 1+2) ---"
echo "  Backgrounds: UC115_V2 + UC093_V3 + viral-only controls"

if [[ "${DRY_RUN}" == true ]]; then
    echo "  sbatch --dependency=afterok:JOB1:JOB2 ${SCRIPT_DIR}/mix_and_assemble.sh all"
    JOB3="DRY_JOB3"
else
    JOB3=$(sbatch --parsable \
        --dependency="afterok:${JOB1}:${JOB2}" \
        "${SCRIPT_DIR}/mix_and_assemble.sh" all)
    echo "  Submitted: ${JOB3} (depends on ${JOB1},${JOB2})"
fi

# =============================================================================
# Job 4: Label contigs with ground truth (depends on job 3)
# =============================================================================
echo ""
echo "--- Job 4: Ground truth labeling (depends on job 3) ---"

# Create an inline labeling script that processes all assemblies
LABEL_SCRIPT="${SCRIPT_DIR}/logs/label_all_contigs_$$.sh"

cat > "${LABEL_SCRIPT}" << 'LABEL_EOF'
#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH -t 02:00:00
#SBATCH --mem=16G
#SBATCH -J label_contigs
#SBATCH -o logs/label_contigs_%j.out
#SBATCH -e logs/label_contigs_%j.err

set -euo pipefail

PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
SPIKE_IN_DIR="${PROJ_DIR}/data/spike_in"
SCRIPT_DIR="${PROJ_DIR}/scripts/05_spike_in"
VIRAL="${SPIKE_IN_DIR}/all_viral_genomes.fasta"
NEGATIVES="${SPIKE_IN_DIR}/all_negative_controls.fasta"

module load minimap2/2.29

echo "[$(date)] Labeling all Track B assemblies..."

# Label co-assemblies for each background
for bg in UC115_V2 UC093_V3; do
    for cov in 0.1 0.5 1 5 10 50; do
        asm_dir="${SPIKE_IN_DIR}/assemblies/${bg}/assembly_cov${cov}"
        contigs="${asm_dir}/contigs_filtered.fasta"
        gt_out="${asm_dir}/ground_truth.tsv"

        if [[ ! -f "${contigs}" ]]; then
            echo "  SKIP: ${bg}/assembly_cov${cov} (no filtered contigs)"
            continue
        fi

        if [[ -f "${gt_out}" ]]; then
            echo "  SKIP: ${bg}/assembly_cov${cov} (ground truth exists)"
            continue
        fi

        echo "  Labeling: ${bg}/assembly_cov${cov}"
        python3 "${SCRIPT_DIR}/label_spike_in_contigs.py" \
            --contigs "${contigs}" \
            --viral-genomes "${VIRAL}" \
            --negative-genomes "${NEGATIVES}" \
            --output "${gt_out}"
    done
done

# Label viral-only controls
for cov in 0.1 0.5 1 5 10 50; do
    ctl_dir="${SPIKE_IN_DIR}/assemblies/control/assembly_control_cov${cov}"
    contigs="${ctl_dir}/contigs_filtered.fasta"
    gt_out="${ctl_dir}/ground_truth.tsv"

    if [[ ! -f "${contigs}" ]]; then
        echo "  SKIP: control/assembly_control_cov${cov} (no filtered contigs)"
        continue
    fi

    if [[ -f "${gt_out}" ]]; then
        echo "  SKIP: control/assembly_control_cov${cov} (ground truth exists)"
        continue
    fi

    echo "  Labeling: control/assembly_control_cov${cov}"
    python3 "${SCRIPT_DIR}/label_spike_in_contigs.py" \
        --contigs "${contigs}" \
        --viral-genomes "${VIRAL}" \
        --negative-genomes "${NEGATIVES}" \
        --output "${gt_out}"
done

echo ""
echo "[$(date)] All labeling complete."

# Summary
echo ""
echo "================================================================"
echo "Track B Labeling Summary"
echo "================================================================"
echo ""
printf "%-12s %-8s %8s %8s %8s %8s\n" "Background" "Coverage" "Total" "Viral" "Bacterial" "Background"
echo "------------------------------------------------------------------------"
for bg in UC115_V2 UC093_V3 control; do
    for cov in 0.1 0.5 1 5 10 50; do
        if [[ "${bg}" == "control" ]]; then
            gt="${SPIKE_IN_DIR}/assemblies/${bg}/assembly_control_cov${cov}/ground_truth.tsv"
        else
            gt="${SPIKE_IN_DIR}/assemblies/${bg}/assembly_cov${cov}/ground_truth.tsv"
        fi
        if [[ -f "${gt}" ]]; then
            total=$(tail -n+2 "${gt}" | wc -l)
            n_viral=$(awk -F'\t' '$2=="viral"' "${gt}" | wc -l)
            n_bact=$(awk -F'\t' '$2=="bacterial"' "${gt}" | wc -l)
            n_bg=$(awk -F'\t' '$2=="background"' "${gt}" | wc -l)
            printf "%-12s %-8s %8d %8d %8d %8d\n" "${bg}" "${cov}x" "${total}" "${n_viral}" "${n_bact}" "${n_bg}"
        fi
    done
done
LABEL_EOF

chmod +x "${LABEL_SCRIPT}"

if [[ "${DRY_RUN}" == true ]]; then
    echo "  sbatch --dependency=afterok:JOB3 ${LABEL_SCRIPT}"
    JOB4="DRY_JOB4"
else
    JOB4=$(sbatch --parsable \
        --dependency="afterok:${JOB3}" \
        "${LABEL_SCRIPT}")
    echo "  Submitted: ${JOB4} (depends on ${JOB3})"
fi

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "================================================================"
echo "Track B SLURM Dependency Chain"
echo "================================================================"
echo ""
echo "  Job 1 (simulate reads):   ${JOB1}"
echo "  Job 2 (host removal):     ${JOB2}      [parallel with job 1]"
echo "  Job 3 (mix + assemble):   ${JOB3}   [depends on 1+2]"
echo "  Job 4 (label contigs):    ${JOB4}   [depends on 3]"
echo ""
echo "  Estimated wall-clock: 8-12 hours"
echo ""
echo "  Monitor: squeue -u \$USER"
echo "  Cancel:  scancel ${JOB1:-} ${JOB2:-} ${JOB3:-} ${JOB4:-}"
echo ""
echo "  Verification after completion:"
echo "    1. Host-removed reads:  ls results/test_real/samtools/*_host_removed_*.fastq.gz"
echo "    2. Simulated reads:     ls data/spike_in/simulated_reads/IMGVR_*/  (4 new dirs)"
echo "    3. Pooled reads:        ls data/spike_in/simulated_reads/pooled_*.fastq.gz"
echo "    4. Assemblies:          ls data/spike_in/assemblies/{UC115_V2,UC093_V3,control}/"
echo "    5. Ground truth:        find data/spike_in/assemblies/ -name ground_truth.tsv | wc -l  (expect 18)"
