#!/usr/bin/env bash
# =============================================================================
# simulate_reads.sh — InSilicoSeq Read Simulation for Spike-In Benchmark
# =============================================================================
# Generates simulated Illumina NovaSeq reads from each viral spike-in genome
# at 6 coverage depths (0.1x-50x). Outputs per-genome and pooled read sets.
#
# Container: InSilicoSeq via Singularity (auto-pulled from BioContainers)
#
# Usage: sbatch simulate_reads.sh
#        COVERAGE_LEVELS="1 5 10" sbatch simulate_reads.sh  # override coverages
# =============================================================================

#SBATCH --cpus-per-task=4
#SBATCH -t 02:00:00
#SBATCH --mem=16G
#SBATCH -J iss_simulate
#SBATCH -o logs/iss_simulate_%j.out
#SBATCH -e logs/iss_simulate_%j.err

set -euo pipefail

# --- Configuration -----------------------------------------------------------
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
SPIKE_IN_DIR="${PROJ_DIR}/data/spike_in"
OUTDIR="${SPIKE_IN_DIR}/simulated_reads"
CONTAINER="${PROJ_DIR}/containers/insilicoseq/insilicoseq.sif"
READ_LENGTH=150
COVERAGE_LEVELS="${COVERAGE_LEVELS:-0.1 0.5 1 5 10 50}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

mkdir -p "${OUTDIR}" logs

# Bind the project directory (VBENCH_BIND) into the containers
export APPTAINER_BIND="${VBENCH_BIND:?set VBENCH_BIND (see config/paths.example.sh)},${OUTDIR}"

# --- Pull container if needed ------------------------------------------------
if [[ ! -f "${CONTAINER}" ]]; then
    echo "[$(date)] Pulling InSilicoSeq container..."
    mkdir -p "$(dirname "${CONTAINER}")"
    singularity pull "${CONTAINER}" https://depot.galaxyproject.org/singularity/insilicoseq%3A2.0.1--pyh7cba7a3_0
fi

# --- Collect viral genome FASTAs ---------------------------------------------
# Only viral genomes, not negative controls
# Includes 4 MetaVR v5 phages added in 20-genome panel expansion
GENOME_DIRS="lactobacillus_phage gardnerella_phage megasphaera_phage fannyhessea_phage sneathia_phage hpv anellovirus herpesvirus"
GENOME_FILES=()

for dir in ${GENOME_DIRS}; do
    for f in "${SPIKE_IN_DIR}/${dir}"/*.fasta; do
        [[ -f "$f" ]] && GENOME_FILES+=("$f")
    done
done

echo "[$(date)] Simulating reads for ${#GENOME_FILES[@]} viral genomes"
echo "  Coverage levels: ${COVERAGE_LEVELS}"
echo "  Read length:     ${READ_LENGTH} bp PE"
echo "  Error model:     NovaSeq"
echo "  Output:          ${OUTDIR}"
echo ""

# =============================================================================
# Simulate reads per genome × coverage
# =============================================================================
for genome_fasta in "${GENOME_FILES[@]}"; do
    accession=$(basename "${genome_fasta}" .fasta)
    genome_dir="${OUTDIR}/${accession}"
    mkdir -p "${genome_dir}"

    # Get genome size
    genome_size=$(grep -v "^>" "${genome_fasta}" | tr -d '\n' | wc -c)

    echo "  ${accession} (${genome_size} bp):"

    for cov in ${COVERAGE_LEVELS}; do
        prefix="${genome_dir}/${accession}_cov${cov}"

        # Skip if already exists
        if [[ -f "${prefix}_R1.fastq.gz" ]] && [[ -s "${prefix}_R1.fastq.gz" ]]; then
            echo "    ${cov}x: SKIP (exists)"
            continue
        fi

        # Calculate number of reads
        # n_reads = (genome_size * coverage) / (2 * read_length)
        # Use awk for float arithmetic
        n_reads=$(awk "BEGIN {n = int((${genome_size} * ${cov}) / (2 * ${READ_LENGTH})); if (n < 1) n = 1; print n}")

        echo "    ${cov}x: ${n_reads} read pairs"

        # Run InSilicoSeq (single-threaded for reliability; genomes are small)
        singularity exec \
            --bind "${PROJ_DIR}:${PROJ_DIR}" \
            --bind /tmp:/tmp \
            --pwd "${PROJ_DIR}" \
            "${CONTAINER}" \
            iss generate \
                --genomes "${genome_fasta}" \
                --model novaseq \
                --n_reads "${n_reads}" \
                --cpus 1 \
                --output "${prefix}" \
            2>/dev/null || true

        # ISS outputs uncompressed FASTQ — gzip them if both exist
        if [[ -f "${prefix}_R1.fastq" ]] && [[ -f "${prefix}_R2.fastq" ]]; then
            gzip -f "${prefix}_R1.fastq"
            gzip -f "${prefix}_R2.fastq"
        elif [[ -f "${prefix}_R1.fastq" ]]; then
            # ISS sometimes fails on very low read counts — clean up partial
            rm -f "${prefix}_R1.fastq"
            echo "         ✗ ISS produced partial output (R1 only), skipping"
            continue
        fi

        # Clean up ISS temp files
        rm -f "${prefix}".iss.tmp.* "${prefix}"_abundance.txt 2>/dev/null || true

        # Verify
        if [[ -f "${prefix}_R1.fastq.gz" ]] && [[ -s "${prefix}_R1.fastq.gz" ]]; then
            actual_reads=$(zcat "${prefix}_R1.fastq.gz" | awk 'END{print NR/4}')
            echo "         ✓ ${actual_reads} reads generated"
        else
            echo "         ✗ FAILED" >&2
        fi
    done
done

# =============================================================================
# Pool reads per coverage level (all genomes combined)
# =============================================================================
echo ""
echo "[$(date)] Pooling reads per coverage level..."

# Remove old pooled files first (may contain stale data from old 10-genome panel)
rm -f "${OUTDIR}"/pooled_cov*_R?.fastq.gz

for cov in ${COVERAGE_LEVELS}; do
    pooled_r1="${OUTDIR}/pooled_cov${cov}_R1.fastq.gz"
    pooled_r2="${OUTDIR}/pooled_cov${cov}_R2.fastq.gz"

    # Concatenate all per-genome reads at this coverage
    find "${OUTDIR}" -name "*_cov${cov}_R1.fastq.gz" -not -name "pooled_*" -exec cat {} + > "${pooled_r1}"
    find "${OUTDIR}" -name "*_cov${cov}_R2.fastq.gz" -not -name "pooled_*" -exec cat {} + > "${pooled_r2}"

    n_reads_pooled=$(zcat "${pooled_r1}" 2>/dev/null | awk 'NR%4==1' | wc -l)
    echo "  cov${cov}: ${n_reads_pooled} pooled read pairs"
done

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "[$(date)] Read simulation complete."
echo ""
echo "  Per-genome reads: ${OUTDIR}/{accession}/{accession}_cov{X}_R{1,2}.fastq.gz"
echo "  Pooled reads:     ${OUTDIR}/pooled_cov{X}_R{1,2}.fastq.gz"
echo ""
echo "  Next step: sbatch mix_and_assemble.sh"
