#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --mem=72G
#SBATCH -t 24:00:00
#SBATCH -J all_cpu_tools
#SBATCH -o logs/all_cpu_tools_%j.out
#SBATCH -e logs/all_cpu_tools_%j.err
set -euo pipefail
# ============================================================================
# Run all 10 CPU-based virus ID tools on real metagenome contigs.
#
# Adapted from run_cpu_tools.sh (spike-in benchmark).
# Produces native output formats (not WtP .list files) so that
# evaluate_metagenome.py can parse results with full coverage of all contigs.
#
# Usage:
#   sbatch run_secondary_cpu_tools.sh SAMPLE_ID
#   e.g. sbatch run_secondary_cpu_tools.sh UC028_V2
# ============================================================================

# --- Sample ID from argument ------------------------------------------------
SAMPLE="${1:?Usage: sbatch $0 SAMPLE_ID (e.g. UC028_V2)}"

PROJ="${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}"
WTP="${PROJ}/What_the_Phage"
IMG="${WTP}/singularity_images"
DBS="${WTP}/nextflow-autodownload-databases"
CPUS="${SLURM_CPUS_PER_TASK:-12}"

FASTA="${PROJ}/results/test_real/spades/${SAMPLE}_contigs.fasta"
OUTDIR="${PROJ}/results/test_real/full_run/${SAMPLE}"
mkdir -p "${OUTDIR}" logs

# Temp directory for tools needing writable HOME/cache
export TMPDIR="${SNIC_TMP:-/tmp}/all_cpu_tools_$$"
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
echo "All CPU tools: ${SAMPLE} real metagenome"
echo "Input: ${FASTA} (${N_CONTIGS} contigs)"
echo "Output: ${OUTDIR}"
echo "CPUs: ${CPUS}, TMPDIR: ${TMPDIR}"
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
# 1. DeepVirFinder
# ═══════════════════════════════════════════════════════════════════
run_tool "DeepVirFinder"
DVF_OUT="${OUTDIR}/deepvirfinder"
if ! skip_if_done_glob "DeepVirFinder" "${DVF_OUT}/*_dvfpred.txt"; then
mkdir -p "${DVF_OUT}"
apptainer exec "${IMG}/deepvirfinder_0.1.img" bash -c "
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
if ! skip_if_done_glob "VirSorter" "${VS1_OUT}/work/VIRSorter_global-phage-signal.csv"; then
mkdir -p "${VS1_OUT}"
VS1_DB="${DBS}/virsorter/virsorter-data"
if [[ -d "${VS1_DB}" ]]; then
    apptainer exec \
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
fi
fi

# ═══════════════════════════════════════════════════════════════════
# 3. VirSorter2
# ═══════════════════════════════════════════════════════════════════
run_tool "VirSorter2"
VS2_OUT="${OUTDIR}/virsorter2"
if ! skip_if_done "VirSorter2" "${VS2_OUT}/final-viral-score.tsv"; then
VS2_DB="${DBS}/virsorter2-db/db"
rm -rf "${VS2_OUT}"
mkdir -p "${VS2_OUT}"
if [[ -d "${VS2_DB}" ]]; then
    # VirSorter2's internal Snakemake hardcodes /tmp for hmmsearch scratch IO
    # (mktemp -d /tmp/vs2-XXX). Bind SLURM scratch to container /tmp so
    # hmmsearch temp files use real disk instead of limited container tmpfs.
    VS2_TMP="${TMPDIR}/vs2_scratch"
    mkdir -p "${VS2_TMP}/tmp" "${VS2_TMP}/config" "${VS2_TMP}/home"
    apptainer exec \
        --bind "${VS2_TMP}/tmp:/tmp" \
        "${IMG}/virsorter-2_2.2.1--fa935f8.img" bash -c "
        export XDG_CONFIG_HOME='${VS2_TMP}/config'
        export TMPDIR='${VS2_TMP}/tmp'
        export HOME='${VS2_TMP}/home'
        virsorter run \
            -d ${VS2_DB} \
            -w ${VS2_OUT} \
            -i ${FASTA} \
            -j ${CPUS} \
            --provirus-off
    " && echo "[$(date)] VirSorter2 complete" || echo "[$(date)] VirSorter2 FAILED (exit $?)"
else
    echo "[$(date)] VirSorter2 SKIPPED — database not found at ${VS2_DB}"
fi
fi

# ═══════════════════════════════════════════════════════════════════
# 4. VirFinder
# ═══════════════════════════════════════════════════════════════════
run_tool "VirFinder"
VF_OUT="${OUTDIR}/virfinder"
if ! skip_if_done "VirFinder" "${VF_OUT}/virfinder_results.tsv"; then
mkdir -p "${VF_OUT}"
apptainer exec "${IMG}/virfinder_0.2.img" Rscript -e "
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
if ! skip_if_done "PPR-Meta" "${PPR_OUT}/pprmeta_results.csv"; then
mkdir -p "${PPR_OUT}"
PPR_DEPS="${DBS}/pprmeta/PPR-Meta"
if [[ -d "${PPR_DEPS}" ]]; then
    apptainer exec "${IMG}/ppr-meta_0.3.1.img" bash -c "
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
fi
fi

# ═══════════════════════════════════════════════════════════════════
# 6. VIBRANT
# ═══════════════════════════════════════════════════════════════════
run_tool "VIBRANT"
VIB_OUT="${OUTDIR}/vibrant"
if ! skip_if_done_glob "VIBRANT" "${VIB_OUT}/VIBRANT_results_*/VIBRANT_phages_*/*.phages_combined.txt"; then
mkdir -p "${VIB_OUT}"
VIB_DB="${DBS}/Vibrant/database.tar.gz"
if [[ -f "${VIB_DB}" ]]; then
    apptainer exec "${IMG}/vibrant_0.5.img" bash -c "
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
fi
fi

# ═══════════════════════════════════════════════════════════════════
# 7. Seeker
# ═══════════════════════════════════════════════════════════════════
run_tool "Seeker"
SKR_OUT="${OUTDIR}/seeker"
if ! skip_if_done "Seeker" "${SKR_OUT}/seeker_results.tsv"; then
mkdir -p "${SKR_OUT}"
apptainer exec "${IMG}/seeker_0.1.img" bash -c "
    predict-metagenome ${FASTA} > ${SKR_OUT}/seeker_results.tsv
" && echo "[$(date)] Seeker complete" || echo "[$(date)] Seeker FAILED (exit $?)"
fi

# ═══════════════════════════════════════════════════════════════════
# 8. MetaPhinder
# ═══════════════════════════════════════════════════════════════════
run_tool "MetaPhinder"
MPH_OUT="${OUTDIR}/metaphinder"
if ! skip_if_done_glob "MetaPhinder" "${MPH_OUT}/output.txt"; then
mkdir -p "${MPH_OUT}"
apptainer exec "${IMG}/metaphinder_0.1.img" \
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
if ! skip_if_done "Sourmash" "${SM_OUT}/sourmash_results.csv"; then
SM_SPLIT="${SM_OUT}/split_fastas"
mkdir -p "${SM_SPLIT}"
SM_DB="${DBS}/sourmash/phages.sbt.zip"

# Split multi-FASTA into per-contig files
awk '/^>/{if(f) close(f); f=sprintf("'"${SM_SPLIT}"'/%s.fa", substr($1,2)); print > f; next} {print >> f}' "${FASTA}"

if [[ -f "${SM_DB}" ]]; then
    apptainer exec "${IMG}/sourmash_4.5.0--e12a57a.img" bash -c "
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
fi
fi

# ═══════════════════════════════════════════════════════════════════
# 10. geNomad
# ═══════════════════════════════════════════════════════════════════
run_tool "geNomad"
GN_OUT="${OUTDIR}/genomad"
if ! skip_if_done_glob "geNomad" "${GN_OUT}/*_summary/*_virus_summary.tsv"; then
GN_SIF="${PROJ}/containers/genomad/genomad.sif"
GN_DB="${PROJ}/databases/genomad_db/genomad_db"
mkdir -p "${GN_OUT}"
if [[ -f "${GN_SIF}" && -d "${GN_DB}" ]]; then
    apptainer exec \
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
fi
fi

# ═══════════════════════════════════════════════════════════════════
echo ""
echo "================================================================"
echo "All CPU tools finished: $(date)"
echo "Results: ${OUTDIR}"
echo ""
echo "Next steps:"
echo "  1. Run GPU tools: sbatch scripts/06_tool_execution/run_secondary_gpu_tools.sh ${SAMPLE}"
echo "  2. Run secondary benchmark:"
echo "     python scripts/07_evaluation/evaluate_metagenome.py \\"
echo "       --ground-truth results/test_real/ground_truth/${SAMPLE}/ground_truth_with_kraken2.tsv \\"
echo "       --results-dir results/test_real/full_run/${SAMPLE}/ \\"
echo "       --output-dir results/test_real/secondary_benchmark/${SAMPLE}/ \\"
echo "       --sample-id ${SAMPLE} --lactobacillus-test"
echo "================================================================"

# Cleanup temp
rm -rf "${TMPDIR}"
