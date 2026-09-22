#!/bin/bash
# Master orchestrator: submits VIRGO2 + VISTA mgCST pipeline with SLURM dependencies.
#
# Assigns metagenomic Community State Types (mgCSTs) to all 13 shotgun samples
# from Happel et al. 2021 (PRJNA767784) using VIRGO2 gene catalog mapping and
# the VISTA mgCST classifier (Williams et al. 2026).
#
# VIRGO2 and VISTA are reused from the existing installation at:
#   $VISTA_TOOLS_DIR (see config/paths.example.sh)
#
# Usage: bash scripts/08_mgcst/run_mgcst_pipeline.sh

set -euo pipefail

PROJ="${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}"
SCRIPTDIR="${PROJ}/scripts/08_mgcst"

cd "${PROJ}"
mkdir -p logs

echo "============================================"
echo "  mgCST Pipeline -- viral_bench (13 samples)"
echo "============================================"
echo ""

# Verify sample list exists
SAMPLE_LIST="${SCRIPTDIR}/00_sample_list.txt"
if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "ERROR: Sample list not found: ${SAMPLE_LIST}"
    exit 1
fi
NSAMPLES=$(wc -l < "${SAMPLE_LIST}")
echo "Samples: ${NSAMPLES} (from ${SAMPLE_LIST})"
echo ""

# Step 1: Prepare reads (concatenate R1+R2 for VIRGO2 single-end input)
JOB_PREP=$(sbatch --parsable "${SCRIPTDIR}/01_prepare_reads.sh")
echo "[1/4] Prepare reads:    Job ${JOB_PREP}"

# Step 2: VIRGO2 map (array job, depends on step 1)
JOB_MAP=$(sbatch --parsable \
    --dependency=afterok:${JOB_PREP} \
    --array=1-${NSAMPLES} \
    "${SCRIPTDIR}/02_virgo2_map.sh")
echo "[2/4] VIRGO2 map array: Job ${JOB_MAP} (${NSAMPLES} tasks, after ${JOB_PREP})"

# Step 3: VIRGO2 compile (depends on all map jobs completing)
JOB_COMPILE=$(sbatch --parsable \
    --dependency=afterok:${JOB_MAP} \
    "${SCRIPTDIR}/03_virgo2_compile.sh")
echo "[3/4] VIRGO2 compile:   Job ${JOB_COMPILE} (after ${JOB_MAP})"

# Step 4: Run VISTA (depends on compile)
JOB_VISTA=$(sbatch --parsable \
    --dependency=afterok:${JOB_COMPILE} \
    "${SCRIPTDIR}/04_run_vista.sh")
echo "[4/4] Run VISTA:        Job ${JOB_VISTA} (after ${JOB_COMPILE})"

echo ""
echo "============================================"
echo "  Pipeline submitted successfully!"
echo "============================================"
echo ""
echo "Dependency chain:"
echo "  Prepare(${JOB_PREP}) -> Map[1-${NSAMPLES}](${JOB_MAP}) -> Compile(${JOB_COMPILE}) -> VISTA(${JOB_VISTA})"
echo ""
echo "Monitor:  squeue -u \$USER"
echo "Details:  sacct -j ${JOB_PREP},${JOB_MAP},${JOB_COMPILE},${JOB_VISTA}"
echo ""
echo "Expected outputs (after ~2.5 hours):"
echo "  results/mgcst/mgCSTs.csv           -- mgCST assignment + theta score per sample"
echo "  results/mgcst/relabund_w_mgCSTs.csv -- relative abundances with mgCST labels"
echo "  results/mgcst/vista/               -- all VISTA output files + heatmap"
