#!/usr/bin/env bash
# =============================================================================
# run_expansion_pipeline.sh — Top-level orchestrator for the expansion
# =============================================================================
# Submits the 8-stage expansion pipeline as a SLURM dependency chain:
#
#   01 symlink reads          (no SLURM; runs in driver shell)
#   02 host removal (array)
#   --- Track A path ---
#   03 metaSPAdes (PRJNA1170175 only by default; array=1-24)
#   04 Track A novel-panel chain (CheckV array + merge + selection)
#   --- CST + Track B path (parallel where possible) ---
#   05 CST pipeline (mgcst chain + Valencia)
#   06 select Track B backgrounds
#   07 run Track B expansion (eight backgrounds, single 10x coverage)
#   --- Stats ---
#   08 score inclusive metrics + 11 exact native-homolog labels/primary S9
#
# Usage:
#   bash scripts/09_expansion/run_expansion_pipeline.sh
#   bash scripts/09_expansion/run_expansion_pipeline.sh --dry-run
#
# Skip phases by exporting SKIP_TRACK_A=1, SKIP_TRACK_B=1, etc.
# =============================================================================

set -euo pipefail

PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
SCRIPT_DIR="${PROJ_DIR}/scripts/09_expansion"
EXPANSION_DIR="${PROJ_DIR}/data/expansion_cohorts"
ACCOUNT="${ACCOUNT:-${SBATCH_ACCOUNT:?set SBATCH_ACCOUNT}}"
DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true

cd "${PROJ_DIR}"
mkdir -p "${SCRIPT_DIR}/logs"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

echo "============================================================"
log "Expansion pipeline orchestrator"
echo "  Account : ${ACCOUNT}"
echo "  Dry run : ${DRY_RUN}"
echo "============================================================"

# -------- Phase 0: symlink reads (interactive, fast) -----------------------
if [[ ! -f "${EXPANSION_DIR}/array_samples.txt" ]]; then
    log "Phase 0: symlink reads"
    if ${DRY_RUN}; then
        echo "[DRY] bash ${SCRIPT_DIR}/01_symlink_reads.sh"
    else
        bash "${SCRIPT_DIR}/01_symlink_reads.sh"
    fi
else
    log "Phase 0: symlinks already exist (skip)"
fi

N_TOTAL=$(wc -l < "${EXPANSION_DIR}/array_samples.txt" 2>/dev/null || echo 138)

# -------- Phase 1: host removal (array, all 138 samples) -------------------
log "Phase 1: host removal array (${N_TOTAL} tasks)"
if ${DRY_RUN}; then
    echo "[DRY] sbatch -A ${ACCOUNT} --array=1-${N_TOTAL} ${SCRIPT_DIR}/02_run_host_removal.sh"
    JOB_HOST="DRY_HOST"
else
    JOB_HOST=$(sbatch --parsable -A "${ACCOUNT}" --array="1-${N_TOTAL}" \
        --chdir="${PROJ_DIR}" \
        "${SCRIPT_DIR}/02_run_host_removal.sh")
fi
log "  host-removal job: ${JOB_HOST}"

# -------- Phase 2: Track A — metaSPAdes (PRJNA1170175 by default) ----------
if [[ "${SKIP_TRACK_A:-0}" != "1" ]]; then
    # Build PRJNA1170175 sample list if absent (idempotent)
    awk -F'\t' 'NR>1 && $1=="PRJNA1170175" {print $2}' \
        "${EXPANSION_DIR}/samples.tsv" | sort -u > "${EXPANSION_DIR}/track_a_samples.txt"
    N_A=$(wc -l < "${EXPANSION_DIR}/track_a_samples.txt")

    log "Phase 2: Track A metaSPAdes (${N_A} samples)"
    if ${DRY_RUN}; then
        echo "[DRY] sbatch -A ${ACCOUNT} --array=1-${N_A} --dependency=afterok:${JOB_HOST} ${SCRIPT_DIR}/03_run_metaspades.sh"
        JOB_SPADES="DRY_SPADES"
    else
        JOB_SPADES=$(sbatch --parsable -A "${ACCOUNT}" --array="1-${N_A}" \
            --dependency="afterok:${JOB_HOST}" \
            --chdir="${PROJ_DIR}" \
            "${SCRIPT_DIR}/03_run_metaspades.sh")
    fi
    log "  metaSPAdes job: ${JOB_SPADES}"

    log "Phase 3: Track A novel-panel discovery"
    # 04_run_track_a_novel.sh submits step2 + step3 + selection internally;
    # it must run AFTER metaSPAdes completes.
    if ${DRY_RUN}; then
        echo "[DRY] bash ${SCRIPT_DIR}/04_run_track_a_novel.sh (after ${JOB_SPADES})"
        JOB_TRACK_A="DRY_TRACK_A"
    else
        # Wait for metaSPAdes via a tiny dependency-only stub
        JOB_TRACK_A=$(sbatch --parsable -A "${ACCOUNT}" \
            --dependency="afterok:${JOB_SPADES}" \
            --time=01:00:00 --mem=4G --cpus-per-task=1 -J expand_trackA_kick \
            --chdir="${PROJ_DIR}" \
            --wrap "bash ${SCRIPT_DIR}/04_run_track_a_novel.sh")
    fi
    log "  Track A discovery kick: ${JOB_TRACK_A}"
else
    log "Phase 2-3: SKIP_TRACK_A=1 (Track A skipped)"
    JOB_TRACK_A=""
fi

# -------- Phase 4: CST classification (parallel with Track A) --------------
log "Phase 4: CST classification on all 138 samples"
if ${DRY_RUN}; then
    echo "[DRY] bash ${SCRIPT_DIR}/05_run_cst_pipeline.sh (after ${JOB_HOST})"
    JOB_CST_KICK="DRY_CST"
else
    JOB_CST_KICK=$(sbatch --parsable -A "${ACCOUNT}" \
        --dependency="afterok:${JOB_HOST}" \
        --time=01:00:00 --mem=4G --cpus-per-task=1 -J expand_cst_kick \
        --chdir="${PROJ_DIR}" \
        --wrap "bash ${SCRIPT_DIR}/05_run_cst_pipeline.sh")
fi
log "  CST kick job: ${JOB_CST_KICK}"

# -------- Phase 5: Track B (depends on CST output) -------------------------
if [[ "${SKIP_TRACK_B:-0}" != "1" ]]; then
    log "Phase 5: Track B background selection + run"
    # Background selection runs after CST chain completes (no easy dependency
    # tracking on submitted-but-pending children; rely on a script that
    # polls for the Valencia output. The CST kick job submits the chain;
    # the Valencia step is the last in that chain — so we depend on the
    # kick completing AND the output existing).

    # We cannot trivially express "after the entire CST chain" via SLURM
    # without capturing the inner job IDs. The kick job submits an inner
    # chain whose final job ID is not known at this orchestrator level.
    # Practical solution: a small wait-for-file SBATCH that polls until
    # mgCSTs_with_CSTs.csv appears, then triggers selection.

    if ${DRY_RUN}; then
        echo "[DRY] sbatch -A ${ACCOUNT} ${SCRIPT_DIR}/06_select_track_b_backgrounds.sh"
        JOB_BG="DRY_BG"
        echo "[DRY] bash ${SCRIPT_DIR}/07_run_track_b_expanded.sh"
        JOB_TB="DRY_TB"
    else
        JOB_WAIT=$(sbatch --parsable -A "${ACCOUNT}" \
            --time=04:00:00 --mem=4G --cpus-per-task=1 \
            -J expand_wait_cst \
            --chdir="${PROJ_DIR}" \
            --wrap "while [[ ! -s ${PROJ_DIR}/results/mgcst_expansion/mgCSTs_with_CSTs.csv ]]; do sleep 60; done; echo CST output ready")
        log "  CST wait job  : ${JOB_WAIT}"

        JOB_BG=$(sbatch --parsable -A "${ACCOUNT}" \
            --dependency="afterok:${JOB_WAIT}" \
            --chdir="${PROJ_DIR}" \
            "${SCRIPT_DIR}/06_select_track_b_backgrounds.sh")
        log "  Background sel: ${JOB_BG}"

        JOB_TB=$(sbatch --parsable -A "${ACCOUNT}" \
            --dependency="afterok:${JOB_BG}" \
            --time=01:00:00 --mem=4G --cpus-per-task=1 -J expand_trackB_kick \
            --chdir="${PROJ_DIR}" \
            --wrap "bash ${SCRIPT_DIR}/07_run_track_b_expanded.sh")
        log "  Track B kick  : ${JOB_TB}"
    fi
else
    log "Phase 5: SKIP_TRACK_B=1 (Track B skipped)"
    JOB_TB=""
fi

# -------- Phase 6: primary Track B aggregation ------------------------------------------
if [[ -n "${JOB_TB}" ]] && [[ "${SKIP_STATS:-0}" != "1" ]]; then
    log "Phase 6: primary Track B score and aggregate (after all selected tool records exist)"
    if ${DRY_RUN}; then
        echo "[DRY] python3 ${SCRIPT_DIR}/08_score_and_aggregate_track_b.py; bash ${SCRIPT_DIR}/11_native_phage_exact.sh (after all selected tool records exist)"
    else
        # Scope completion checks to the curated manifest. Historical saved run
        # directories must never satisfy readiness or leak into aggregation.
        # The exact native-homolog step writes labels and primary Table S9.
        sbatch -A "${ACCOUNT}" \
            --dependency="afterok:${JOB_TB}" \
            --time=04:00:00 --mem=32G --cpus-per-task=8 -J expand_dmcc_kick \
            --chdir="${PROJ_DIR}" \
            --wrap "set -e; while ! python3 ${SCRIPT_DIR}/track_b_manifest.py --check-completion; do sleep 60; done; module load SciPy-bundle/2024.05-gfbf-2024a; python3 ${SCRIPT_DIR}/08_score_and_aggregate_track_b.py; bash ${SCRIPT_DIR}/11_native_phage_exact.sh"
        log "  dMCC bootstrap submitted"
    fi
fi

echo ""
echo "============================================================"
log "All phases submitted."
echo "  Monitor      : squeue -u \$USER"
echo "  Track A out  : ${EXPANSION_DIR}/track_a_panel/selected/expansion_novel_genomes.fasta"
echo "  CST out      : ${PROJ_DIR}/results/mgcst_expansion/mgCSTs_with_CSTs.csv"
echo "  Track B BG   : ${EXPANSION_DIR}/track_b_backgrounds.tsv"
echo "  Track B asm  : ${PROJ_DIR}/results/expansion/track_b/assemblies/"
echo "  Track B bench: ${PROJ_DIR}/results/expansion/track_b/benchmark/"
echo "  Primary S9   : ${PROJ_DIR}/results/tables/table_s9b_track_b_multibackground_excl.tsv"
echo "============================================================"
