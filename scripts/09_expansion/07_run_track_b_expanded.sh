#!/usr/bin/env bash
# =============================================================================
# 07_run_track_b_expanded.sh — Track B on the eight curated backgrounds
# =============================================================================
# Wraps the existing Track B chain to run at single 10x coverage across the eight
# selected expansion-cohort backgrounds (4 CST-I + 4 CST-IV-B).
#
# Pipeline:
#   1. simulate_reads.sh        — reuse, single 10x coverage only
#   2. mix_and_assemble.sh      — once per background (case "*" branch covers all SRR* IDs)
#   3. label_track_b_contigs.sh — same minimap2 asm5 logic as the original
#   4. run_all_track_b.sh       — submits CPU + GPU tool runs over the eight expansion BG
#
# All HOST_REMOVED_DIR / OUTDIR / RES_DIR are env-overridden so expansion
# outputs land under results/expansion/track_b/ and never touch the original
# spike_in/assemblies/ tree.
#
# Usage:
#   bash scripts/09_expansion/07_run_track_b_expanded.sh
#   bash scripts/09_expansion/07_run_track_b_expanded.sh --dry-run
# =============================================================================

set -euo pipefail

PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
EXPANSION_DIR="${PROJ_DIR}/data/expansion_cohorts"
TRACK_B_BG_TSV="${EXPANSION_DIR}/track_b_backgrounds.tsv"
SPIKE_IN_DIR="${PROJ_DIR}/data/spike_in"
ASM_OUT_DIR="${PROJ_DIR}/results/expansion/track_b/assemblies"
RES_OUT_DIR="${PROJ_DIR}/results/expansion/track_b/benchmark"
HOST_REMOVED_DIR="${PROJ_DIR}/results/expansion/clean_reads"
SCRIPT_DIR="${PROJ_DIR}/scripts/05_spike_in"
TOOL_DIR="${PROJ_DIR}/scripts/06_tool_execution"
ACCOUNT="${ACCOUNT:-${SBATCH_ACCOUNT:?set SBATCH_ACCOUNT}}"
DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true

mkdir -p "${ASM_OUT_DIR}" "${RES_OUT_DIR}" "${PROJ_DIR}/scripts/09_expansion/logs"

if [[ ! -s "${TRACK_B_BG_TSV}" ]]; then
    echo "ERROR: background list missing: ${TRACK_B_BG_TSV}" >&2
    echo "Run 06_select_track_b_backgrounds.sh first." >&2
    exit 1
fi

# Read selected backgrounds (skip header)
BACKGROUND_RUNS=$(python3 "${PROJ_DIR}/scripts/09_expansion/track_b_manifest.py" --manifest "${TRACK_B_BG_TSV}")
mapfile -t BACKGROUNDS <<< "${BACKGROUND_RUNS}"
N_BG=${#BACKGROUNDS[@]}

if [[ ${N_BG} -lt 1 ]]; then
    echo "ERROR: no backgrounds selected" >&2
    exit 1
fi

echo "============================================================"
echo "[$(date)] Track B expansion (${N_BG} backgrounds, 10x only)"
echo "  Background list: ${TRACK_B_BG_TSV}"
echo "  Backgrounds    : ${BACKGROUNDS[*]}"
echo "  Asm output     : ${ASM_OUT_DIR}"
echo "  Bench output   : ${RES_OUT_DIR}"
echo "  Host-rm reads  : ${HOST_REMOVED_DIR}"
echo "  Account        : ${ACCOUNT}"
echo "  Dry run        : ${DRY_RUN}"
echo "============================================================"

# Verify clean reads exist for each background
missing=()
for bg in "${BACKGROUNDS[@]}"; do
    if [[ ! -s "${HOST_REMOVED_DIR}/${bg}_host_removed_R1.fastq.gz" ]]; then
        missing+=("${bg}")
    fi
done
if [[ ${#missing[@]} -gt 0 ]]; then
    echo "ERROR: clean reads missing for ${#missing[@]} backgrounds:" >&2
    printf '  %s\n' "${missing[@]}" >&2
    echo "Run 02_run_host_removal.sh first." >&2
    exit 1
fi

submit() {
    local desc="$1"; shift
    if ${DRY_RUN}; then
        echo "[DRY] sbatch -A ${ACCOUNT} $* ($desc)"
        echo "DRY"
    else
        sbatch --parsable -A "${ACCOUNT}" "$@" --chdir="${PROJ_DIR}"
    fi
}

# --- Step 1: simulate reads at single 10x coverage --------------------------
# simulate_reads.sh respects COVERAGE_LEVELS env var; exposes pooled FASTQs
# under data/spike_in/simulated_reads/pooled_cov10_R{1,2}.fastq.gz. If the 10x
# pool already exists from the original Track B run, we skip resubmission.
POOLED_R1="${SPIKE_IN_DIR}/simulated_reads/pooled_cov10_R1.fastq.gz"
if [[ -s "${POOLED_R1}" ]]; then
    JOB_SIM="0"
    echo ""
    echo "[1/4] Simulate reads      : SKIP (pooled_cov10 already exists)"
else
    echo ""
    echo "[1/4] Simulate reads (10x only)..."
    JOB_SIM=$(submit "simulate" \
        --export="ALL,COVERAGE_LEVELS=10" \
        "${SCRIPT_DIR}/simulate_reads.sh")
    echo "      job: ${JOB_SIM}"
fi

# --- Step 2: mix_and_assemble per background (10x only) ---------------------
# Each background runs as its own SBATCH job (small footprint); they all
# depend on simulate having completed. mix_and_assemble.sh's case "*" branch
# now accepts any background ID with reads at HOST_REMOVED_DIR.
echo ""
echo "[2/4] Mix + co-assemble per background..."
ASM_DEPS=""
# Build a single dependency string only if simulation actually ran.
# Empty bash arrays expanded as "${arr[@]:-}" yield a literal "" that
# sbatch interprets as a filename argument, hence "Unable to open file".
if [[ "${JOB_SIM}" != "0" ]]; then
    SIM_DEP_ARG="--dependency=afterok:${JOB_SIM}"
else
    SIM_DEP_ARG=""
fi
for bg in "${BACKGROUNDS[@]}"; do
    if [[ -n "${SIM_DEP_ARG}" ]]; then
        JOB_BG=$(submit "mix+assemble ${bg}" \
            "${SIM_DEP_ARG}" \
            --export="ALL,COVERAGE_LEVELS=10,HOST_REMOVED_DIR=${HOST_REMOVED_DIR},OUTDIR=${ASM_OUT_DIR},PROJ_DIR=${PROJ_DIR},SPIKE_IN_DIR=${SPIKE_IN_DIR}" \
            "${SCRIPT_DIR}/mix_and_assemble.sh" "${bg}")
    else
        JOB_BG=$(submit "mix+assemble ${bg}" \
            --export="ALL,COVERAGE_LEVELS=10,HOST_REMOVED_DIR=${HOST_REMOVED_DIR},OUTDIR=${ASM_OUT_DIR},PROJ_DIR=${PROJ_DIR},SPIKE_IN_DIR=${SPIKE_IN_DIR}" \
            "${SCRIPT_DIR}/mix_and_assemble.sh" "${bg}")
    fi
    [[ -z "${ASM_DEPS}" ]] && ASM_DEPS="${JOB_BG}" || ASM_DEPS="${ASM_DEPS}:${JOB_BG}"
    echo "      ${bg}: ${JOB_BG}"
done

# --- Step 3: label contigs (after all assemblies complete) ------------------
echo ""
echo "[3/4] Label expansion-cohort contigs..."
JOB_LABEL=$(submit "label" \
    --dependency="afterok:${ASM_DEPS}" \
    --export="ALL,ASM_DIR=${ASM_OUT_DIR},BACKGROUNDS_LIST=${TRACK_B_BG_TSV},COVERAGES=10" \
    "${PROJ_DIR}/scripts/09_expansion/07b_label_track_b.sh")
echo "      job: ${JOB_LABEL}"

# --- Step 4: run all 14 tools (after labels exist) --------------------------
echo ""
echo "[4/4] Run 14 tools across all expansion backgrounds..."
BG_STR="${BACKGROUNDS[*]}"
JOB_TOOLS=$(submit "run_all_track_b" \
    --dependency="afterok:${JOB_LABEL}" \
    --export="ALL,BACKGROUNDS_OVERRIDE=${BG_STR},COVERAGES_OVERRIDE=10,ASM_DIR=${ASM_OUT_DIR},RES_DIR=${RES_OUT_DIR},RUN_CONTROLS=false,PROJ_DIR=${PROJ_DIR}" \
    --wrap "bash ${TOOL_DIR}/run_all_track_b.sh")
echo "      job: ${JOB_TOOLS}"

echo ""
echo "============================================================"
echo "Track B expansion submission chain:"
echo "  Simulate reads     : ${JOB_SIM}"
echo "  Mix+assemble (xN)  : (${ASM_DEPS})"
echo "  Label contigs      : ${JOB_LABEL}"
echo "  Run 14 tools       : ${JOB_TOOLS}"
echo ""
echo "  Monitor: squeue -u \$USER"
echo "  Outputs: ${RES_OUT_DIR}/trackB_<bg>_cov10/<tool>/"
echo "============================================================"
