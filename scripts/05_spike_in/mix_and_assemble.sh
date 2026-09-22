#!/usr/bin/env bash
# =============================================================================
# mix_and_assemble.sh — Mix Spike-In Reads with Background + Co-Assemble
# =============================================================================
# Combines pooled simulated viral reads with real host-removed background reads,
# then co-assembles with metaSPAdes at each coverage level. Also produces
# assembly-only controls (viral reads without background).
#
# Now parameterized for multiple background samples:
#   UC115_V2 — CST-I  (96.7% L. crispatus, low-diversity)
#   UC093_V3 — CST-IV (47.8% Gardnerella + diverse, high-diversity)
#
# Depends on:
#   - simulate_reads.sh output (pooled_cov{X}_R{1,2}.fastq.gz)
#   - Host-removed reads in results/test_real/samtools/
#
# Usage:
#   sbatch mix_and_assemble.sh UC115_V2          # Single background
#   sbatch mix_and_assemble.sh UC093_V3          # Single background
#   sbatch mix_and_assemble.sh control            # Viral-only controls
#   sbatch mix_and_assemble.sh all                # All backgrounds + controls
# =============================================================================

#SBATCH --cpus-per-task=12
#SBATCH -t 2-00:00:00
#SBATCH --mem=72G
#SBATCH -J spikein_assemble
#SBATCH -o logs/spikein_assemble_%j.out
#SBATCH -e logs/spikein_assemble_%j.err

set -euo pipefail

# --- Configuration -----------------------------------------------------------
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
SPIKE_IN_DIR="${SPIKE_IN_DIR:-${PROJ_DIR}/data/spike_in}"
READS_DIR="${READS_DIR:-${SPIKE_IN_DIR}/simulated_reads}"
OUTDIR="${OUTDIR:-${SPIKE_IN_DIR}/assemblies}"
COVERAGE_LEVELS="${COVERAGE_LEVELS:-0.1 0.5 1 5 10 50}"
MIN_CONTIG_LENGTH="${MIN_CONTIG_LENGTH:-1500}"
THREADS="${SLURM_CPUS_PER_TASK:-12}"
MEMORY="${MEMORY:-60}"  # GB for SPAdes

# Background sample to process (positional argument)
BACKGROUND="${1:-all}"

# Host-removed reads directory (override for expansion-cohort backgrounds)
HOST_REMOVED_DIR="${HOST_REMOVED_DIR:-${PROJ_DIR}/results/test_real/samtools}"

mkdir -p "${OUTDIR}" logs

# --- Load SPAdes module ------------------------------------------------------
module load SPAdes/4.2.0-GCC-13.3.0

# =============================================================================
# Helper: filter contigs by minimum length
# =============================================================================
filter_contigs() {
    local asm_dir="$1"
    if [[ -f "${asm_dir}/contigs.fasta" ]]; then
        python3 -c "
min_len = ${MIN_CONTIG_LENGTH}
infile = '${asm_dir}/contigs.fasta'
outfile = '${asm_dir}/contigs_filtered.fasta'
kept = 0
total = 0
current_id = None
current_seq = []
with open(infile) as fin, open(outfile, 'w') as fout:
    for line in fin:
        line = line.strip()
        if line.startswith('>'):
            if current_id and len(''.join(current_seq)) >= min_len:
                fout.write(f'>{current_id}\n')
                fout.write(''.join(current_seq) + '\n')
                kept += 1
            total += 1
            current_id = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
    if current_id and len(''.join(current_seq)) >= min_len:
        fout.write(f'>{current_id}\n')
        fout.write(''.join(current_seq) + '\n')
        kept += 1
    total += 1
print(f'  Filtered {infile}: {kept}/{total} contigs >= {min_len} bp')
"
    fi
}

# =============================================================================
# Helper: run co-assembly for a given background sample
# =============================================================================
run_background_assembly() {
    local sample="$1"
    local bg_r1="${HOST_REMOVED_DIR}/${sample}_host_removed_R1.fastq.gz"
    local bg_r2="${HOST_REMOVED_DIR}/${sample}_host_removed_R2.fastq.gz"
    local bg_outdir="${OUTDIR}/${sample}"

    # Validate background reads
    if [[ ! -f "${bg_r1}" ]] || [[ ! -f "${bg_r2}" ]]; then
        echo "ERROR: Host-removed reads not found for ${sample}:" >&2
        echo "  Expected: ${bg_r1}" >&2
        echo "  Expected: ${bg_r2}" >&2
        echo "  Run run_host_removal.sh first" >&2
        return 1
    fi

    mkdir -p "${bg_outdir}"

    local bg_reads
    bg_reads=$(zcat "${bg_r1}" | awk 'END{print NR/4}')
    echo ""
    echo "################################################################"
    echo "[$(date)] Background: ${sample} (${bg_reads} read pairs)"
    echo "################################################################"
    echo ""

    for cov in ${COVERAGE_LEVELS}; do
        echo "================================================================"
        echo "[$(date)] ${sample} — Coverage: ${cov}x"
        echo "================================================================"

        local pooled_r1="${READS_DIR}/pooled_cov${cov}_R1.fastq.gz"
        local pooled_r2="${READS_DIR}/pooled_cov${cov}_R2.fastq.gz"

        if [[ ! -f "${pooled_r1}" ]]; then
            echo "  WARNING: Pooled reads not found for cov${cov}, skipping" >&2
            continue
        fi

        # --- Mix reads ---
        local mixed_dir="${bg_outdir}/mixed_cov${cov}"
        mkdir -p "${mixed_dir}"

        echo "  Mixing reads..."
        cat "${pooled_r1}" "${bg_r1}" > "${mixed_dir}/mixed_R1.fastq.gz"
        cat "${pooled_r2}" "${bg_r2}" > "${mixed_dir}/mixed_R2.fastq.gz"

        local mixed_reads
        mixed_reads=$(zcat "${mixed_dir}/mixed_R1.fastq.gz" | awk 'END{print NR/4}')
        echo "  Mixed read pairs: ${mixed_reads}"

        # --- Co-assemble (mixed) ---
        local asm_dir="${bg_outdir}/assembly_cov${cov}"
        if [[ -f "${asm_dir}/contigs.fasta" ]]; then
            echo "  Assembly exists, skipping co-assembly"
        else
            echo "  Co-assembling with metaSPAdes..."
            spades.py --meta -k 33,55,77 \
                -1 "${mixed_dir}/mixed_R1.fastq.gz" \
                -2 "${mixed_dir}/mixed_R2.fastq.gz" \
                -o "${asm_dir}" \
                -t "${THREADS}" \
                -m "${MEMORY}" \
                2>&1 | tail -5
        fi

        # --- Filter contigs ---
        filter_contigs "${asm_dir}"

        # --- Clean up mixed reads (large files) ---
        rm -f "${mixed_dir}/mixed_R1.fastq.gz" "${mixed_dir}/mixed_R2.fastq.gz"
        rmdir "${mixed_dir}" 2>/dev/null || true

        echo ""
    done

    # --- Summary for this background ---
    echo "[$(date)] ${sample} assemblies complete."
    for cov in ${COVERAGE_LEVELS}; do
        local asm="${bg_outdir}/assembly_cov${cov}/contigs_filtered.fasta"
        local n_asm
        n_asm=$(grep -c "^>" "${asm}" 2>/dev/null || echo "N/A")
        printf "  cov%5s:  %s contigs\n" "${cov}" "${n_asm}"
    done
    echo ""
}

# =============================================================================
# Helper: run viral-only control assemblies (no background)
# =============================================================================
run_control_assembly() {
    local ctl_outdir="${OUTDIR}/control"
    mkdir -p "${ctl_outdir}"

    echo ""
    echo "################################################################"
    echo "[$(date)] Control assemblies (viral reads only, no background)"
    echo "################################################################"
    echo ""

    for cov in ${COVERAGE_LEVELS}; do
        echo "================================================================"
        echo "[$(date)] Control — Coverage: ${cov}x"
        echo "================================================================"

        local pooled_r1="${READS_DIR}/pooled_cov${cov}_R1.fastq.gz"
        local pooled_r2="${READS_DIR}/pooled_cov${cov}_R2.fastq.gz"

        if [[ ! -f "${pooled_r1}" ]]; then
            echo "  WARNING: Pooled reads not found for cov${cov}, skipping" >&2
            continue
        fi

        local ctl_dir="${ctl_outdir}/assembly_control_cov${cov}"
        if [[ -f "${ctl_dir}/contigs.fasta" ]]; then
            echo "  Control assembly exists, skipping"
        else
            echo "  Control assembly (viral reads only)..."
            spades.py --meta -k 33,55,77 \
                -1 "${pooled_r1}" \
                -2 "${pooled_r2}" \
                -o "${ctl_dir}" \
                -t "${THREADS}" \
                -m "${MEMORY}" \
                2>&1 | tail -5
        fi

        filter_contigs "${ctl_dir}"
        echo ""
    done

    # --- Summary ---
    echo "[$(date)] Control assemblies complete."
    for cov in ${COVERAGE_LEVELS}; do
        local ctl="${ctl_outdir}/assembly_control_cov${cov}/contigs_filtered.fasta"
        local n_ctl
        n_ctl=$(grep -c "^>" "${ctl}" 2>/dev/null || echo "N/A")
        printf "  cov%5s:  %s contigs\n" "${cov}" "${n_ctl}"
    done
    echo ""
}

# =============================================================================
# Main dispatch
# =============================================================================
echo "[$(date)] Spike-in mix & co-assembly"
echo "  Background mode: ${BACKGROUND}"
echo "  Coverage levels: ${COVERAGE_LEVELS}"
echo "  Min contig length: ${MIN_CONTIG_LENGTH} bp"
echo "  SPAdes: $(spades.py --version 2>&1 || echo 'loaded')"
echo ""

case "${BACKGROUND}" in
    all)
        run_background_assembly "UC115_V2"
        run_background_assembly "UC093_V3"
        run_control_assembly
        ;;
    control)
        run_control_assembly
        ;;
    UC115_V2|UC093_V3)
        run_background_assembly "${BACKGROUND}"
        ;;
    *)
        # Any other token is treated as an explicit background sample name —
        # required for the expansion-cohort pipeline (scripts/09_expansion/),
        # which submits arbitrary SRR* run accessions as backgrounds. The
        # host-removed reads must exist at
        #   ${HOST_REMOVED_DIR}/${sample}_host_removed_R{1,2}.fastq.gz
        run_background_assembly "${BACKGROUND}"
        ;;
esac

echo ""
echo "[$(date)] All requested assemblies complete."
echo "  Next step: python label_spike_in_contigs.py"
