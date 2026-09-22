#!/usr/bin/env bash
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -J all_gpu_tools
#SBATCH -o logs/all_gpu_tools_%j.out
#SBATCH -e logs/all_gpu_tools_%j.err
set -euo pipefail

# ============================================================================
# Run all 4 GPU-based virus ID tools on real metagenome contigs.
#
# Adapted from run_gpu_tools.sh (spike-in benchmark).
# Produces native output formats for evaluate_metagenome.py.
#
# Usage:
#   sbatch run_secondary_gpu_tools.sh SAMPLE_ID
#   e.g. sbatch run_secondary_gpu_tools.sh UC028_V2
# ============================================================================

# --- Sample ID from argument ------------------------------------------------
SAMPLE="${1:?Usage: sbatch $0 SAMPLE_ID (e.g. UC028_V2)}"

PROJ="${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}"
CPUS="${SLURM_CPUS_PER_TASK:-8}"

FASTA="${PROJ}/results/test_real/spades/${SAMPLE}_contigs.fasta"
OUTDIR="${PROJ}/results/test_real/full_run/${SAMPLE}"
mkdir -p "${OUTDIR}" logs

export TMPDIR="${SNIC_TMP:-/tmp}/all_gpu_tools_$$"
mkdir -p "${TMPDIR}"

# Bind the project directory (VBENCH_BIND) into the containers
export APPTAINER_BIND="${VBENCH_BIND:?set VBENCH_BIND (see config/paths.example.sh)},${TMPDIR}"

# Validate input
if [[ ! -f "${FASTA}" ]]; then
    echo "ERROR: Input FASTA not found: ${FASTA}" >&2
    exit 1
fi

N_CONTIGS=$(grep -c '^>' "${FASTA}")
echo "================================================================"
echo "All GPU tools: ${SAMPLE} real metagenome"
echo "Input: ${FASTA} (${N_CONTIGS} contigs)"
echo "Output: ${OUTDIR}"
echo "CPUs: ${CPUS}, GPU: ${SLURM_GPUS:-unknown}"
echo "Started: $(date)"
echo "================================================================"

run_tool() {
    local tool="$1"
    echo ""
    echo "── ${tool} ──────────────────────────────────"
    echo "[$(date)] Starting ${tool}"
}

skip_if_done() {
    local tool="$1" outfile="$2"
    if [[ -s "${outfile}" ]]; then
        echo "[$(date)] ${tool} SKIPPED — output already exists: ${outfile}"
        return 0
    fi
    return 1
}

skip_if_done_glob() {
    local tool="$1" pattern="$2"
    local match
    match=$(compgen -G "${pattern}" 2>/dev/null | head -1) || true
    if [[ -n "${match}" && -s "${match}" ]]; then
        echo "[$(date)] ${tool} SKIPPED — output already exists: ${match}"
        return 0
    fi
    return 1
}

# ═══════════════════════════════════════════════════════════════════
# 1. HVSeeker
# ═══════════════════════════════════════════════════════════════════
run_tool "HVSeeker"
HV_OUT="${OUTDIR}/hvseeker"
if ! skip_if_done "HVSeeker" "${HV_OUT}/hvseeker_results.tsv"; then
mkdir -p "${HV_OUT}"
apptainer exec --nv \
    "${PROJ}/containers/hvseeker/hvseeker-dna.sif" \
    /opt/conda/envs/HVSeekerDNA/bin/python "${PROJ}/scripts/06_tool_execution/wrappers/hvseeker_wrapper.py" \
        --input "${FASTA}" \
        --output "${HV_OUT}/hvseeker_results.tsv" \
        --model "${PROJ}/containers/hvseeker/models/model_best_acc2_test_model.pt" \
        --hvseeker-dir "${PROJ}/containers/hvseeker/HVSeeker-DNA/HVSeeker-DNA" \
        --threshold 0.5 \
        --min-length 500 \
        --max-length 500 \
&& echo "[$(date)] HVSeeker complete" || echo "[$(date)] HVSeeker FAILED (exit $?)"
fi

# ═══════════════════════════════════════════════════════════════════
# 2. Jaeger
# ═══════════════════════════════════════════════════════════════════
run_tool "Jaeger"
JG_OUT="${OUTDIR}/jaeger"
if ! skip_if_done_glob "Jaeger" "${JG_OUT}/*.tsv"; then
mkdir -p "${JG_OUT}"
apptainer exec --nv \
    "${PROJ}/containers/jaeger/jaeger.sif" \
    jaeger run \
        -i "${FASTA}" \
        -o "${JG_OUT}" \
        --batch 96 \
        --workers "${CPUS}" \
&& echo "[$(date)] Jaeger complete" || echo "[$(date)] Jaeger FAILED (exit $?)"
fi

# ═══════════════════════════════════════════════════════════════════
# 3. TransGINmer
# ═══════════════════════════════════════════════════════════════════
run_tool "TransGINmer"
TG_OUT="${OUTDIR}/transginmer"
if ! skip_if_done "TransGINmer" "${TG_OUT}/transginmer_results.tsv"; then
mkdir -p "${TG_OUT}"

# TransGINmer writes contig.fasta + tokenized data into its install dir.
# --writable-tmpfs puts this in RAM which overflows on large inputs (12k+ contigs).
# Fix: copy the TransGINmer dir to scratch and redirect all large writes there.
TG_WORK="${TMPDIR}/transginmer_scratch"
mkdir -p "${TG_WORK}/tg_dir" "${TG_WORK}/pytmp"

echo "  Copying TransGINmer install to scratch..."
apptainer exec \
    "${PROJ}/containers/transginmer/transginmer.sif" \
    cp -r /opt/transginmer/. "${TG_WORK}/tg_dir/"

echo "  Running TransGINmer from scratch (${TG_WORK})..."
apptainer exec --nv \
    --writable-tmpfs \
    --env TMPDIR="${TG_WORK}/pytmp" \
    "${PROJ}/containers/transginmer/transginmer.sif" \
    /opt/transginmer-env/bin/python "${PROJ}/scripts/06_tool_execution/wrappers/transginmer_wrapper.py" \
        --input "${FASTA}" \
        --output "${TG_OUT}/transginmer_results.tsv" \
        --transginmer-dir "${TG_WORK}/tg_dir" \
        --threshold 0.5 \
        --min-length 100 \
&& echo "[$(date)] TransGINmer complete" || echo "[$(date)] TransGINmer FAILED (exit $?)"

rm -rf "${TG_WORK}"
fi

# ═══════════════════════════════════════════════════════════════════
# 4. ViraLM
# ═══════════════════════════════════════════════════════════════════
run_tool "ViraLM"
VLM_OUT="${OUTDIR}/viralm"
if ! skip_if_done "ViraLM" "${VLM_OUT}/viralm_results.tsv"; then
mkdir -p "${VLM_OUT}"
apptainer exec --nv \
    "${PROJ}/containers/viralm/viralm.sif" bash -c "
    export HF_HOME='${TMPDIR}/hf_cache'
    export TRANSFORMERS_CACHE='${TMPDIR}/hf_cache'
    export HF_MODULES_CACHE='${TMPDIR}/hf_modules'
    mkdir -p \${HF_HOME} \${HF_MODULES_CACHE}
    viralm predict \
        --input ${FASTA} \
        --output ${VLM_OUT}/viralm_results.tsv \
        --threads ${CPUS}
" && echo "[$(date)] ViraLM complete" || echo "[$(date)] ViraLM FAILED (exit $?)"
fi

# ═══════════════════════════════════════════════════════════════════
echo ""
echo "================================================================"
echo "All GPU tools finished: $(date)"
echo "Results: ${OUTDIR}"
echo ""
echo "Next steps:"
echo "  python scripts/07_evaluation/evaluate_metagenome.py \\"
echo "    --ground-truth results/test_real/ground_truth/${SAMPLE}/ground_truth_with_kraken2.tsv \\"
echo "    --results-dir results/test_real/full_run/${SAMPLE}/ \\"
echo "    --output-dir results/test_real/secondary_benchmark/${SAMPLE}/ \\"
echo "    --sample-id ${SAMPLE} --lactobacillus-test"
echo "================================================================"

rm -rf "${TMPDIR}"
