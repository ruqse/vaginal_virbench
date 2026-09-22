#!/usr/bin/env bash
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -J gpu_tools_benchmark
#SBATCH -o logs/gpu_tools_%j.out
#SBATCH -e logs/gpu_tools_%j.err
set -euo pipefail

# ── Run 4 GPU-based virus ID tools on spike-in fragment FASTAs ──
# Usage: sbatch run_gpu_tools.sh <fragment_fasta>
# Example: sbatch run_gpu_tools.sh data/spike_in/fragments/L1500_fragments.fasta

FASTA="${1:?Usage: sbatch run_gpu_tools.sh <fragment_fasta>}"
FASTA="$(realpath "${FASTA}")"
BASENAME="$(basename "${FASTA}" .fasta)"

PROJ="${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}"
CPUS="${SLURM_CPUS_PER_TASK:-8}"

OUTDIR="${2:-${PROJ}/results/spike_in_benchmark/${BASENAME}}"
OUTDIR="$(realpath "${OUTDIR}")"
mkdir -p "${OUTDIR}" logs

export TMPDIR="${SNIC_TMP:-/tmp}/gpu_tools_$$"
mkdir -p "${TMPDIR}"

# Bind the project directory (VBENCH_BIND) into the containers
export APPTAINER_BIND="${VBENCH_BIND:?set VBENCH_BIND (see config/paths.example.sh)},${TMPDIR}"

# ── Resource tracking setup ──────────────────────────────────────
RESOURCE_TSV="${OUTDIR}/resource_usage_gpu.tsv"
TIME_LOG="${TMPDIR}/time_log.txt"
GPU_LOG="${TMPDIR}/gpu_monitor.csv"
printf "tool\twall_time_sec\tpeak_rss_kb\tgpu_mem_peak_mb\tgpu_util_max_pct\texit_code\tskipped\n" > "${RESOURCE_TSV}"

echo "================================================================"
echo "GPU tools benchmark: ${BASENAME}"
echo "Input: ${FASTA} ($(grep -c '^>' "${FASTA}") fragments)"
echo "Output: ${OUTDIR}"
echo "CPUs: ${CPUS}, GPU: ${SLURM_GPUS:-unknown}"
echo "Resource log: ${RESOURCE_TSV}"
echo "Started: $(date)"
echo "================================================================"

run_tool() {
    local tool="$1"
    echo ""
    echo "── ${tool} ──────────────────────────────────"
    echo "[$(date)] Starting ${tool}"
}

# Run a command with /usr/bin/time -v and record resource usage
run_timed() {
    local tool="$1"; shift
    local t_start=${SECONDS}
    /usr/bin/time -v "$@" 2> "${TIME_LOG}"; local rc=$?
    local t_end=${SECONDS}
    local wall=$((t_end - t_start))
    local peak_rss
    peak_rss=$(grep "Maximum resident set size" "${TIME_LOG}" 2>/dev/null | awk '{print $NF}') || true
    # Return values for caller to append GPU columns
    _TIMED_WALL=${wall}
    _TIMED_RSS=${peak_rss:-0}
    _TIMED_RC=${rc}
    echo "[$(date)] ${tool} resources: wall=${wall}s peak_rss=${peak_rss:-0}KB exit=${rc}"
    return ${rc}
}

# GPU monitoring: background nvidia-smi loop
start_gpu_monitor() {
    rm -f "${GPU_LOG}"
    (
        while true; do
            nvidia-smi --query-gpu=memory.used,utilization.gpu \
                --format=csv,noheader,nounits >> "${GPU_LOG}" 2>/dev/null
            sleep 5
        done
    ) &
    GPU_MON_PID=$!
}

stop_gpu_monitor() {
    kill ${GPU_MON_PID} 2>/dev/null
    wait ${GPU_MON_PID} 2>/dev/null || true
    local peak_mem=0 max_util=0
    if [[ -s "${GPU_LOG}" ]]; then
        peak_mem=$(awk -F', ' '{if($1+0>max) max=$1+0} END{print int(max)}' "${GPU_LOG}")
        max_util=$(awk -F', ' '{if($2+0>max) max=$2+0} END{print int(max)}' "${GPU_LOG}")
    fi
    _GPU_MEM_PEAK=${peak_mem}
    _GPU_UTIL_MAX=${max_util}
}

# Record a skipped tool
record_skip() {
    local tool="$1"
    printf "%s\t0\t0\t0\t0\t0\tyes\n" "${tool}" >> "${RESOURCE_TSV}"
}

# Record a completed tool with all metrics
record_gpu_tool() {
    local tool="$1"
    printf "%s\t%d\t%s\t%s\t%s\t%d\tno\n" \
        "${tool}" "${_TIMED_WALL}" "${_TIMED_RSS}" \
        "${_GPU_MEM_PEAK}" "${_GPU_UTIL_MAX}" "${_TIMED_RC}" >> "${RESOURCE_TSV}"
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
if skip_if_done "HVSeeker" "${HV_OUT}/hvseeker_results.tsv"; then
    record_skip "HVSeeker"
else
mkdir -p "${HV_OUT}"
start_gpu_monitor
run_timed "HVSeeker" apptainer exec --nv \
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
stop_gpu_monitor
record_gpu_tool "HVSeeker"
echo "  GPU: peak_mem=${_GPU_MEM_PEAK}MB max_util=${_GPU_UTIL_MAX}%"
fi

# ═══════════════════════════════════════════════════════════════════
# 2. Jaeger
# ═══════════════════════════════════════════════════════════════════
run_tool "Jaeger"
JG_OUT="${OUTDIR}/jaeger"
if skip_if_done_glob "Jaeger" "${JG_OUT}/*.tsv"; then
    record_skip "Jaeger"
else
mkdir -p "${JG_OUT}"
start_gpu_monitor
run_timed "Jaeger" apptainer exec --nv \
    "${PROJ}/containers/jaeger/jaeger.sif" \
    jaeger run \
        -i "${FASTA}" \
        -o "${JG_OUT}" \
        --batch 96 \
        --workers "${CPUS}" \
&& echo "[$(date)] Jaeger complete" || echo "[$(date)] Jaeger FAILED (exit $?)"
# Note: Jaeger silently skips sequences <2048 bp. For L500/L1000/L1500 bins,
# expect empty or no output — this is a valid benchmark result.
stop_gpu_monitor
record_gpu_tool "Jaeger"
echo "  GPU: peak_mem=${_GPU_MEM_PEAK}MB max_util=${_GPU_UTIL_MAX}%"
fi

# ═══════════════════════════════════════════════════════════════════
# 3. TransGINmer
# ═══════════════════════════════════════════════════════════════════
run_tool "TransGINmer"
TG_OUT="${OUTDIR}/transginmer"
if skip_if_done "TransGINmer" "${TG_OUT}/transginmer_results.tsv"; then
    record_skip "TransGINmer"
else
mkdir -p "${TG_OUT}"
start_gpu_monitor
run_timed "TransGINmer" apptainer exec --nv \
    --writable-tmpfs \
    "${PROJ}/containers/transginmer/transginmer.sif" \
    /opt/transginmer-env/bin/python "${PROJ}/scripts/06_tool_execution/wrappers/transginmer_wrapper.py" \
        --input "${FASTA}" \
        --output "${TG_OUT}/transginmer_results.tsv" \
        --threshold 0.5 \
        --min-length 100 \
&& echo "[$(date)] TransGINmer complete" || echo "[$(date)] TransGINmer FAILED (exit $?)"
stop_gpu_monitor
record_gpu_tool "TransGINmer"
echo "  GPU: peak_mem=${_GPU_MEM_PEAK}MB max_util=${_GPU_UTIL_MAX}%"
fi

# ═══════════════════════════════════════════════════════════════════
# 4. ViraLM
# ═══════════════════════════════════════════════════════════════════
run_tool "ViraLM"
VLM_OUT="${OUTDIR}/viralm"
if skip_if_done "ViraLM" "${VLM_OUT}/viralm_results.tsv"; then
    record_skip "ViraLM"
else
mkdir -p "${VLM_OUT}"
start_gpu_monitor
run_timed "ViraLM" apptainer exec --nv \
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
stop_gpu_monitor
record_gpu_tool "ViraLM"
echo "  GPU: peak_mem=${_GPU_MEM_PEAK}MB max_util=${_GPU_UTIL_MAX}%"
fi

# ═══════════════════════════════════════════════════════════════════
echo ""
echo "================================================================"
echo "All GPU tools finished: $(date)"
echo "Results: ${OUTDIR}"
echo ""
echo "Resource usage summary:"
column -t -s $'\t' "${RESOURCE_TSV}" 2>/dev/null || cat "${RESOURCE_TSV}"
echo "================================================================"

rm -rf "${TMPDIR}"
