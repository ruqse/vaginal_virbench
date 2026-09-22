#!/usr/bin/env bash
# =============================================================================
# 04_run_track_a_novel.sh — Mine ANI-novel phages from expansion assemblies
# =============================================================================
# Wraps the existing 05_spike_in/novel_discovery/ pipeline (CheckV + DIAMOND +
# BLASTn vs MetaVR v5) to ingest expansion-cohort assemblies. Reuses the
# existing scripts unmodified — they accept OUTDIR/INDEX/PROJ env-var overrides
# (added in this expansion phase, backward-compatible).
#
# Pipeline:
#   1. Build assembly index from data/expansion_cohorts/assemblies/*.contigs.fa.gz
#   2. Submit step2_checkv_array.sh with the new index
#   3. After CheckV completes, submit step3_merge_characterize.sh
#
# Output base: data/expansion_cohorts/track_a_panel/
#   ├── assembly_index.tsv
#   ├── filtered/<sample>.fasta
#   ├── checkv/<sample>/quality_summary.tsv
#   ├── candidates/
#   │   ├── candidate_summary.tsv          # per-contig ANI to MetaVR v5
#   │   ├── hq_viral_contigs.fasta
#   │   ├── diamond_blastx.tsv
#   │   └── blastn_metavr_v5.tsv
#   └── selected/                          # filled by 04b script (post-CheckV)
#       └── expansion_novel_genomes.fasta
#
# Usage:
#   bash scripts/09_expansion/04_run_track_a_novel.sh             # submit chain
#   bash scripts/09_expansion/04_run_track_a_novel.sh --dry-run
# =============================================================================

set -euo pipefail

PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
EXPANSION_DIR="${PROJ_DIR}/data/expansion_cohorts"
ASM_DIR="${EXPANSION_DIR}/assemblies"
TRACK_A_DIR="${EXPANSION_DIR}/track_a_panel"
SCRIPT_DIR="${PROJ_DIR}/scripts/05_spike_in/novel_discovery"

# SLURM accounts (per staging-README convention)
ACCOUNT="${ACCOUNT:-${SBATCH_ACCOUNT:?set SBATCH_ACCOUNT}}"

DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true

mkdir -p "${TRACK_A_DIR}"/{filtered,checkv,candidates,selected} \
         "${PROJ_DIR}/scripts/09_expansion/logs"

echo "============================================================"
echo "[$(date)] Track A novel-panel discovery (expansion cohorts)"
echo "  Assemblies in: ${ASM_DIR}"
echo "  Pipeline out : ${TRACK_A_DIR}"
echo "  Account      : ${ACCOUNT}"
echo "  Dry run      : ${DRY_RUN}"
echo "============================================================"
echo ""

# --- 1. Build assembly index ------------------------------------------------
INDEX="${TRACK_A_DIR}/assembly_index.tsv"
ls "${ASM_DIR}"/*.contigs.fa.gz 2>/dev/null | sort > "${INDEX}.tmp"

N_ASM=$(wc -l < "${INDEX}.tmp")
if [[ ${N_ASM} -eq 0 ]]; then
    echo "ERROR: no assemblies in ${ASM_DIR}/*.contigs.fa.gz" >&2
    echo "Run 03_run_metaspades.sh first." >&2
    rm -f "${INDEX}.tmp"
    exit 1
fi
mv "${INDEX}.tmp" "${INDEX}"
echo "[$(date)] Assembly index: ${INDEX} (${N_ASM} assemblies)"

# --- 2. Submit step2 CheckV array (with env-var overrides) ------------------
echo ""
echo "[$(date)] Submitting CheckV array (1-${N_ASM})..."

if ${DRY_RUN}; then
    echo "[DRY] sbatch -A ${ACCOUNT} --array=1-${N_ASM} \\"
    echo "      --export=ALL,PROJ=${PROJ_DIR},OUTDIR=${TRACK_A_DIR},INDEX=${INDEX} \\"
    echo "      ${SCRIPT_DIR}/step2_checkv_array.sh"
    JOB_CHECKV="DRY"
else
    JOB_CHECKV=$(sbatch --parsable \
        -A "${ACCOUNT}" \
        --array="1-${N_ASM}" \
        --export="ALL,PROJ=${PROJ_DIR},OUTDIR=${TRACK_A_DIR},INDEX=${INDEX}" \
        --chdir="${PROJ_DIR}" \
        "${SCRIPT_DIR}/step2_checkv_array.sh")
    echo "  Job ${JOB_CHECKV} submitted"
fi

# --- 3. Submit step3 merge+characterize (depends on step2) ------------------
echo ""
echo "[$(date)] Submitting merge+characterize (depends on ${JOB_CHECKV})..."

if ${DRY_RUN}; then
    echo "[DRY] sbatch -A ${ACCOUNT} --dependency=afterok:${JOB_CHECKV} \\"
    echo "      --export=ALL,PROJ=${PROJ_DIR},OUTDIR=${TRACK_A_DIR} \\"
    echo "      ${SCRIPT_DIR}/step3_merge_characterize.sh"
    JOB_MERGE="DRY"
else
    JOB_MERGE=$(sbatch --parsable \
        -A "${ACCOUNT}" \
        --dependency="afterok:${JOB_CHECKV}" \
        --export="ALL,PROJ=${PROJ_DIR},OUTDIR=${TRACK_A_DIR}" \
        --chdir="${PROJ_DIR}" \
        "${SCRIPT_DIR}/step3_merge_characterize.sh")
    echo "  Job ${JOB_MERGE} submitted"
fi

# --- 4. Submit panel selection (depends on step3) ---------------------------
SEL_SCRIPT="${PROJ_DIR}/scripts/09_expansion/04b_select_expansion_panel.sh"
echo ""
echo "[$(date)] Submitting expansion-panel selector (depends on ${JOB_MERGE})..."

if ${DRY_RUN}; then
    echo "[DRY] sbatch -A ${ACCOUNT} --dependency=afterok:${JOB_MERGE} ${SEL_SCRIPT}"
    JOB_SEL="DRY"
else
    JOB_SEL=$(sbatch --parsable \
        -A "${ACCOUNT}" \
        --dependency="afterok:${JOB_MERGE}" \
        --export="ALL,PROJ=${PROJ_DIR}" \
        --chdir="${PROJ_DIR}" \
        "${SEL_SCRIPT}")
    echo "  Job ${JOB_SEL} submitted"
fi

echo ""
echo "============================================================"
echo "[$(date)] Track A submission chain:"
echo "  CheckV array (${N_ASM} tasks): ${JOB_CHECKV}"
echo "  Merge+characterize           : ${JOB_MERGE}"
echo "  Expansion-panel selection    : ${JOB_SEL}"
echo ""
echo "  Monitor: squeue -u \$USER"
echo "  Outputs: ${TRACK_A_DIR}/candidates/candidate_summary.tsv"
echo "           ${TRACK_A_DIR}/selected/expansion_novel_genomes.fasta"
echo "============================================================"
