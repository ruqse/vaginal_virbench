#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --mem=72G
#SBATCH -t 24:00:00
#SBATCH -J cpu_tools_benchmark
#SBATCH -o logs/cpu_tools_%j.out
#SBATCH -e logs/cpu_tools_%j.err
set -euo pipefail

# ── Run 10 CPU-based virus ID tools on spike-in fragment FASTAs ──
# Usage: sbatch run_cpu_tools.sh <fragment_fasta>
# Example: sbatch run_cpu_tools.sh data/spike_in/fragments/L1500_fragments.fasta

FASTA="${1:?Usage: sbatch run_cpu_tools.sh <fragment_fasta>}"
FASTA="$(realpath "${FASTA}")"
BASENAME="$(basename "${FASTA}" .fasta)"

PROJ="${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}"
WTP="${PROJ}/What_the_Phage"
IMG="${WTP}/singularity_images"
DBS="${WTP}/nextflow-autodownload-databases"
CPUS="${SLURM_CPUS_PER_TASK:-12}"

OUTDIR="${2:-${PROJ}/results/spike_in_benchmark/${BASENAME}}"
OUTDIR="$(realpath "${OUTDIR}")"
mkdir -p "${OUTDIR}" logs

# Temp directory for tools needing writable HOME/cache
export TMPDIR="${SNIC_TMP:-/tmp}/cpu_tools_$$"
mkdir -p "${TMPDIR}"

# Bind the project directory (VBENCH_BIND) into the containers
export APPTAINER_BIND="${VBENCH_BIND:?set VBENCH_BIND (see config/paths.example.sh)},${TMPDIR}"

# ── Resource tracking setup ──────────────────────────────────────
RESOURCE_TSV="${OUTDIR}/resource_usage.tsv"
TIME_LOG="${TMPDIR}/time_log.txt"
printf "tool\twall_time_sec\tpeak_rss_kb\texit_code\tskipped\n" > "${RESOURCE_TSV}"

echo "================================================================"
echo "CPU tools benchmark: ${BASENAME}"
echo "Input: ${FASTA} ($(grep -c '^>' "${FASTA}") fragments)"
echo "Output: ${OUTDIR}"
echo "CPUs: ${CPUS}, TMPDIR: ${TMPDIR}"
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
    printf "%s\t%d\t%s\t%d\tno\n" "${tool}" "${wall}" "${peak_rss:-0}" "${rc}" >> "${RESOURCE_TSV}"
    echo "[$(date)] ${tool} resources: wall=${wall}s peak_rss=${peak_rss:-0}KB exit=${rc}"
    return ${rc}
}

# Record a skipped tool
record_skip() {
    local tool="$1"
    printf "%s\t0\t0\t0\tyes\n" "${tool}" >> "${RESOURCE_TSV}"
}

# Skip a tool if its output file already exists and is non-empty.
# Usage: skip_if_done "ToolName" "output_file" && return-or-continue
# For glob patterns, use: skip_if_done_glob "ToolName" "glob_pattern"
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
# 1. DeepVirFinder
# ═══════════════════════════════════════════════════════════════════
run_tool "DeepVirFinder"
DVF_OUT="${OUTDIR}/deepvirfinder"
if skip_if_done_glob "DeepVirFinder" "${DVF_OUT}/*_dvfpred.txt"; then
    record_skip "DeepVirFinder"
else
mkdir -p "${DVF_OUT}"
run_timed "DeepVirFinder" apptainer exec "${IMG}/deepvirfinder_0.1.img" bash -c "
    export THEANO_FLAGS='base_compiledir=${TMPDIR}/theano'
    export KERAS_HOME='${TMPDIR}/keras'
    export HOME='${TMPDIR}'
    cd /DeepVirFinder
    python dvf.py -i ${FASTA} -o ${DVF_OUT} -c ${CPUS}
" && echo "[$(date)] DeepVirFinder complete" || echo "[$(date)] DeepVirFinder FAILED (exit $?)"
fi

# ═══════════════════════════════════════════════════════════════════
# 2. VirSorter v1
# ═══════════════════════════════════════════════════════════════════
run_tool "VirSorter"
VS1_OUT="${OUTDIR}/virsorter"
if skip_if_done_glob "VirSorter" "${VS1_OUT}/work/VIRSorter_global-phage-signal.csv"; then
    record_skip "VirSorter"
else
mkdir -p "${VS1_OUT}"
VS1_DB="${DBS}/virsorter/virsorter-data"
if [[ -d "${VS1_DB}" ]]; then
    run_timed "VirSorter" apptainer exec \
        "${IMG}/virsorter_0.1.2.img" \
        wrapper_phage_contigs_sorter_iPlant.pl \
            -f "${FASTA}" \
            -db 2 \
            --wdir "${VS1_OUT}/work" \
            --ncpu "${CPUS}" \
            --data-dir "${VS1_DB}" \
    && echo "[$(date)] VirSorter complete" || echo "[$(date)] VirSorter FAILED (exit $?)"
else
    echo "[$(date)] VirSorter SKIPPED — database not found at ${VS1_DB}"
    record_skip "VirSorter"
fi
fi

# ═══════════════════════════════════════════════════════════════════
# 3. VirSorter2
# ═══════════════════════════════════════════════════════════════════
run_tool "VirSorter2"
VS2_OUT="${OUTDIR}/virsorter2"
if skip_if_done "VirSorter2" "${VS2_OUT}/final-viral-score.tsv"; then
    record_skip "VirSorter2"
else
VS2_DB="${DBS}/virsorter2-db/db"
rm -rf "${VS2_OUT}"
mkdir -p "${VS2_OUT}"
if [[ -d "${VS2_DB}" ]]; then
    # NOTE: --provirus-off was used for the originally published runs; this
    # disabled VirSorter2's provirus boundary classifier and biased the H2
    # marker-vs-sequence prophage contrast. The flag is now removed so that
    # provirus detection runs at developer default. Re-running this step is
    # required to refresh Track A/B/C and multi-evidence VirSorter2 metrics.
    run_timed "VirSorter2" apptainer exec \
        "${IMG}/virsorter-2_2.2.1--fa935f8.img" bash -c "
        export XDG_CONFIG_HOME='${TMPDIR}/vs2_config'
        export TMPDIR='${TMPDIR}/vs2_tmp'
        export HOME='${TMPDIR}'
        mkdir -p \${XDG_CONFIG_HOME} \${TMPDIR}
        virsorter run \
            -d ${VS2_DB} \
            -w ${VS2_OUT} \
            -i ${FASTA} \
            -j ${CPUS}
    " && echo "[$(date)] VirSorter2 complete" || echo "[$(date)] VirSorter2 FAILED (exit $?)"
else
    echo "[$(date)] VirSorter2 SKIPPED — database not found at ${VS2_DB}"
    record_skip "VirSorter2"
fi
fi

# ═══════════════════════════════════════════════════════════════════
# 4. VirFinder
# ═══════════════════════════════════════════════════════════════════
run_tool "VirFinder"
VF_OUT="${OUTDIR}/virfinder"
if skip_if_done "VirFinder" "${VF_OUT}/virfinder_results.tsv"; then
    record_skip "VirFinder"
else
mkdir -p "${VF_OUT}"
run_timed "VirFinder" apptainer exec "${IMG}/virfinder_0.2.img" Rscript -e "
    library(VirFinder)
    predResult <- VF.pred('${FASTA}')
    results <- predResult[order(predResult\$pvalue),]
    write.table(results, file='${VF_OUT}/virfinder_results.tsv',
                sep='\t', row.names=FALSE, quote=FALSE)
" && echo "[$(date)] VirFinder complete" || echo "[$(date)] VirFinder FAILED (exit $?)"
fi

# ═══════════════════════════════════════════════════════════════════
# 5. PPR-Meta
# ═══════════════════════════════════════════════════════════════════
run_tool "PPR-Meta"
PPR_OUT="${OUTDIR}/pprmeta"
if skip_if_done "PPR-Meta" "${PPR_OUT}/pprmeta_results.csv"; then
    record_skip "PPR-Meta"
else
mkdir -p "${PPR_OUT}"
PPR_DEPS="${DBS}/pprmeta/PPR-Meta"
if [[ -d "${PPR_DEPS}" ]]; then
    run_timed "PPR-Meta" apptainer exec "${IMG}/ppr-meta_0.3.1.img" bash -c "
        export MCR_CACHE_ROOT='${TMPDIR}/mcr_cache'
        export HOME='${TMPDIR}/mcr_home'
        export MATLAB_PREFDIR='${TMPDIR}/mcr_home/.matlab'
        mkdir -p \${MCR_CACHE_ROOT} \${HOME} \${MATLAB_PREFDIR}
        cd ${PPR_OUT}
        cp ${PPR_DEPS}/* . 2>/dev/null || true
        ./PPR_Meta ${FASTA} pprmeta_results.csv
    " && echo "[$(date)] PPR-Meta complete" || echo "[$(date)] PPR-Meta FAILED (exit $?)"
else
    echo "[$(date)] PPR-Meta SKIPPED — deps not found at ${PPR_DEPS}"
    record_skip "PPR-Meta"
fi
fi

# ═══════════════════════════════════════════════════════════════════
# 6. VIBRANT
# ═══════════════════════════════════════════════════════════════════
run_tool "VIBRANT"
VIB_OUT="${OUTDIR}/vibrant"
if skip_if_done_glob "VIBRANT" "${VIB_OUT}/VIBRANT_results_*/VIBRANT_phages_*/*.phages_combined.txt"; then
    record_skip "VIBRANT"
else
mkdir -p "${VIB_OUT}"
VIB_DB="${DBS}/Vibrant/database.tar.gz"
if [[ -f "${VIB_DB}" ]]; then
    run_timed "VIBRANT" apptainer exec "${IMG}/vibrant_0.5.img" bash -c "
        cd ${VIB_OUT}
        tar xzf ${VIB_DB}
        VIBRANT_run.py -i ${FASTA} -t ${CPUS} \
            -k database/KEGG_profiles_prokaryotes.HMM \
            -p database/Pfam-A_v32.HMM \
            -v database/VOGDB94_phage.HMM \
            -e database/Pfam-A_plasmid_v32.HMM \
            -a database/Pfam-A_phage_v32.HMM \
            -c database/VIBRANT_categories.tsv \
            -n database/VIBRANT_names.tsv \
            -s database/VIBRANT_KEGG_pathways_summary.tsv \
            -m database/VIBRANT_machine_model.sav \
            -g database/VIBRANT_AMGs.tsv
    " && echo "[$(date)] VIBRANT complete" || echo "[$(date)] VIBRANT FAILED (exit $?)"
else
    echo "[$(date)] VIBRANT SKIPPED — database not found at ${VIB_DB}"
    record_skip "VIBRANT"
fi
fi

# ═══════════════════════════════════════════════════════════════════
# 7. Seeker
# ═══════════════════════════════════════════════════════════════════
run_tool "Seeker"
SKR_OUT="${OUTDIR}/seeker"
if skip_if_done "Seeker" "${SKR_OUT}/seeker_results.tsv"; then
    record_skip "Seeker"
else
mkdir -p "${SKR_OUT}"
run_timed "Seeker" apptainer exec "${IMG}/seeker_0.1.img" bash -c "
    predict-metagenome ${FASTA} > ${SKR_OUT}/seeker_results.tsv
" && echo "[$(date)] Seeker complete" || echo "[$(date)] Seeker FAILED (exit $?)"
fi

# ═══════════════════════════════════════════════════════════════════
# 8. MetaPhinder
# ═══════════════════════════════════════════════════════════════════
run_tool "MetaPhinder"
MPH_OUT="${OUTDIR}/metaphinder"
if skip_if_done_glob "MetaPhinder" "${MPH_OUT}/output.txt"; then
    record_skip "MetaPhinder"
else
mkdir -p "${MPH_OUT}"
run_timed "MetaPhinder" apptainer exec "${IMG}/metaphinder_0.1.img" \
    MetaPhinder.py \
        -i "${FASTA}" \
        -o "${MPH_OUT}" \
        -d /MetaPhinder/database/ALL_140821_hr \
&& echo "[$(date)] MetaPhinder complete" || echo "[$(date)] MetaPhinder FAILED (exit $?)"
fi

# ═══════════════════════════════════════════════════════════════════
# 9. Sourmash
# ═══════════════════════════════════════════════════════════════════
run_tool "Sourmash"
SM_OUT="${OUTDIR}/sourmash"
if skip_if_done "Sourmash" "${SM_OUT}/sourmash_results.csv"; then
    record_skip "Sourmash"
else
SM_SPLIT="${SM_OUT}/split_fastas"
mkdir -p "${SM_SPLIT}"
SM_DB="${DBS}/sourmash/phages.sbt.zip"

# Split multi-FASTA into per-contig files
awk '/^>/{if(f) close(f); f=sprintf("'"${SM_SPLIT}"'/%s.fa", substr($1,2)); print > f; next} {print >> f}' "${FASTA}"

if [[ -f "${SM_DB}" ]]; then
    run_timed "Sourmash" apptainer exec "${IMG}/sourmash_4.5.0--e12a57a.img" bash -c "
        cd ${SM_OUT}
        for fa in split_fastas/*.fa; do
            sourmash sketch dna -p k=21,scaled=100 \"\${fa}\" -o \"\${fa}.sig\" 2>/dev/null || true
        done
        : > sourmash_results.csv
        for sig in split_fastas/*.sig; do
            sourmash search -k 21 \"\${sig}\" ${SM_DB} -o \"\${sig}.tmp\" 2>/dev/null || true
            if [[ -f \"\${sig}.tmp\" ]]; then
                hits=\$(grep -v 'similarity,' \"\${sig}.tmp\" | wc -l)
                if [[ \${hits} -gt 0 ]]; then
                    name=\$(basename \"\${sig}\" .fa.sig)
                    score=\$(grep -v 'similarity,' \"\${sig}.tmp\" | sort -rn -k1 | head -1 | cut -d',' -f1)
                    echo \"\${name},\${score}\" >> sourmash_results.csv
                fi
            fi
        done
    " && echo "[$(date)] Sourmash complete" || echo "[$(date)] Sourmash FAILED (exit $?)"
else
    echo "[$(date)] Sourmash SKIPPED — database not found at ${SM_DB}"
    record_skip "Sourmash"
fi
fi

# ═══════════════════════════════════════════════════════════════════
# 10. geNomad
# ═══════════════════════════════════════════════════════════════════
run_tool "geNomad"
GN_OUT="${OUTDIR}/genomad"
if skip_if_done_glob "geNomad" "${GN_OUT}/*_summary/*_virus_summary.tsv"; then
    record_skip "geNomad"
else
GN_SIF="${PROJ}/containers/genomad/genomad.sif"
GN_DB="${PROJ}/databases/genomad_db/genomad_db"
mkdir -p "${GN_OUT}"
if [[ -f "${GN_SIF}" && -d "${GN_DB}" ]]; then
    run_timed "geNomad" apptainer exec \
        "${GN_SIF}" \
        genomad end-to-end \
            "${FASTA}" \
            "${GN_OUT}" \
            "${GN_DB}" \
            --cleanup \
            --threads "${CPUS}" \
    && echo "[$(date)] geNomad complete" || echo "[$(date)] geNomad FAILED (exit $?)"
else
    echo "[$(date)] geNomad SKIPPED — container or database not found"
    echo "  Container: ${GN_SIF} (exists: $(test -f "${GN_SIF}" && echo yes || echo no))"
    echo "  Database: ${GN_DB} (exists: $(test -d "${GN_DB}" && echo yes || echo no))"
    record_skip "geNomad"
fi
fi

# ═══════════════════════════════════════════════════════════════════
echo ""
echo "================================================================"
echo "All CPU tools finished: $(date)"
echo "Results: ${OUTDIR}"
echo ""
echo "Resource usage summary:"
column -t -s $'\t' "${RESOURCE_TSV}" 2>/dev/null || cat "${RESOURCE_TSV}"
echo "================================================================"

# Cleanup temp
rm -rf "${TMPDIR}"
