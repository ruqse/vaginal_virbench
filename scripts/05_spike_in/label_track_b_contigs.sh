#!/usr/bin/env bash
# =============================================================================
# label_track_b_contigs.sh — Label all Track B co-assemblies with ground truth
# =============================================================================
# Runs label_spike_in_contigs.py on all Track B co-assemblies and controls.
#
# Excludes genomes with HIGH homology risk per homology check results:
#   UC093_V3: vB_Gva_AB1 (MW387018.1) — 80.3% coverage at 99.5% identity
#             with native Gardnerella phage (NODE_361_length_19947_cov_3.884046)
#   UC115_V2: No exclusions (all genomes LOW risk)
#
# Usage:
#   sbatch --dependency=afterok:ASSEMBLY_JOB label_track_b_contigs.sh
# =============================================================================

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
echo ""

# --- UC115_V2: no exclusions (all genomes LOW risk) -------------------------
echo "=== UC115_V2 (no exclusions) ==="
for cov in 0.1 0.5 1 5 10 50; do
    asm_dir="${SPIKE_IN_DIR}/assemblies/UC115_V2/assembly_cov${cov}"
    contigs="${asm_dir}/contigs_filtered.fasta"
    gt_out="${asm_dir}/ground_truth.tsv"
    [[ ! -f "${contigs}" ]] && echo "  SKIP: UC115_V2/cov${cov} (no contigs)" && continue
    [[ -f "${gt_out}" ]] && echo "  SKIP: UC115_V2/cov${cov} (GT exists)" && continue
    echo "  Labeling: UC115_V2/cov${cov}"
    python3 "${SCRIPT_DIR}/label_spike_in_contigs.py" \
        --contigs "${contigs}" --viral-genomes "${VIRAL}" \
        --negative-genomes "${NEGATIVES}" --output "${gt_out}"
done
echo ""

# --- UC093_V3: exclude vB_Gva_AB1 (MW387018.1) — HIGH chimeric risk ---------
echo "=== UC093_V3 (excluding MW387018.1 / vB_Gva_AB1) ==="
for cov in 0.1 0.5 1 5 10 50; do
    asm_dir="${SPIKE_IN_DIR}/assemblies/UC093_V3/assembly_cov${cov}"
    contigs="${asm_dir}/contigs_filtered.fasta"
    gt_out="${asm_dir}/ground_truth.tsv"
    [[ ! -f "${contigs}" ]] && echo "  SKIP: UC093_V3/cov${cov} (no contigs)" && continue
    [[ -f "${gt_out}" ]] && echo "  SKIP: UC093_V3/cov${cov} (GT exists)" && continue
    echo "  Labeling: UC093_V3/cov${cov}"
    python3 "${SCRIPT_DIR}/label_spike_in_contigs.py" \
        --contigs "${contigs}" --viral-genomes "${VIRAL}" \
        --negative-genomes "${NEGATIVES}" \
        --exclude-genomes "MW387018.1" \
        --output "${gt_out}"
done
echo ""

# --- Controls: no exclusions (viral reads only, no background) ---------------
echo "=== Controls (no exclusions) ==="
for cov in 0.1 0.5 1 5 10 50; do
    ctl_dir="${SPIKE_IN_DIR}/assemblies/control/assembly_control_cov${cov}"
    contigs="${ctl_dir}/contigs_filtered.fasta"
    gt_out="${ctl_dir}/ground_truth.tsv"
    [[ ! -f "${contigs}" ]] && echo "  SKIP: control/cov${cov} (no contigs)" && continue
    [[ -f "${gt_out}" ]] && echo "  SKIP: control/cov${cov} (GT exists)" && continue
    echo "  Labeling: control/cov${cov}"
    python3 "${SCRIPT_DIR}/label_spike_in_contigs.py" \
        --contigs "${contigs}" --viral-genomes "${VIRAL}" \
        --negative-genomes "${NEGATIVES}" --output "${gt_out}"
done

echo ""
echo "[$(date)] All labeling complete."
echo ""
echo "Note: UC093_V3 contigs mapping to vB_Gva_AB1 (MW387018.1) labeled as"
echo "  'excluded_chimeric_risk' due to 80.3% homology at 99.5% identity"
echo "  with native Gardnerella phage in UC093_V3 standalone assembly."
echo "  See: results/spike_in/homology_check/UC093_V3/HOMOLOGY_REPORT.md"
