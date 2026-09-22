#!/usr/bin/env bash
# =============================================================================
# 07b_label_track_b.sh — Label expansion-cohort Track B contigs with ground truth
# =============================================================================
# Direct adaptation of scripts/05_spike_in/label_track_b_contigs.sh, scoped to
# the eight expansion-cohort backgrounds at single 10x coverage. Reuses the
# existing label_spike_in_contigs.py (minimap2 asm5) — no logic change.
#
# Backgrounds and exclusions are read from the TSV produced by
# 06_select_track_b_backgrounds.sh. By default no per-background exclusions
# are applied; the existing pre-screen (mash) already removed backgrounds
# whose native phages collide with the spike-in panel.
# =============================================================================

#SBATCH --cpus-per-task=4
#SBATCH -t 02:00:00
#SBATCH --mem=16G
#SBATCH -J expand_label
#SBATCH -o logs/expand_label_%j.out
#SBATCH -e logs/expand_label_%j.err

set -euo pipefail

PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
SPIKE_IN_DIR="${PROJ_DIR}/data/spike_in"
SCRIPT_DIR="${PROJ_DIR}/scripts/05_spike_in"
ASM_DIR="${ASM_DIR:-${PROJ_DIR}/results/expansion/track_b/assemblies}"
BACKGROUNDS_LIST="${BACKGROUNDS_LIST:-${PROJ_DIR}/data/expansion_cohorts/track_b_backgrounds.tsv}"
COVERAGES="${COVERAGES:-10}"
VIRAL="${SPIKE_IN_DIR}/all_viral_genomes.fasta"
NEGATIVES="${SPIKE_IN_DIR}/all_negative_controls.fasta"

module load minimap2/2.29

echo "[$(date)] Labeling expansion-cohort Track B assemblies..."
echo "  ASM_DIR    : ${ASM_DIR}"
echo "  Backgrounds: ${BACKGROUNDS_LIST}"
echo "  Coverages  : ${COVERAGES}"
echo ""

# Read the background list (skip header)
BACKGROUND_RUNS=$(python3 "${PROJ_DIR}/scripts/09_expansion/track_b_manifest.py" --manifest "${BACKGROUNDS_LIST}")
mapfile -t BACKGROUNDS <<< "${BACKGROUND_RUNS}"

for bg in "${BACKGROUNDS[@]}"; do
    for cov in ${COVERAGES}; do
        asm_dir="${ASM_DIR}/${bg}/assembly_cov${cov}"
        contigs="${asm_dir}/contigs_filtered.fasta"
        gt_out="${asm_dir}/ground_truth.tsv"
        if [[ ! -f "${contigs}" ]]; then
            echo "  SKIP: ${bg}/cov${cov} (no contigs)"
            continue
        fi
        if [[ -f "${gt_out}" ]]; then
            echo "  SKIP: ${bg}/cov${cov} (GT exists)"
            continue
        fi
        echo "  Labeling: ${bg}/cov${cov}"
        python3 "${SCRIPT_DIR}/label_spike_in_contigs.py" \
            --contigs "${contigs}" --viral-genomes "${VIRAL}" \
            --negative-genomes "${NEGATIVES}" --output "${gt_out}"
    done
done

echo ""
echo "[$(date)] All expansion labelling complete."
