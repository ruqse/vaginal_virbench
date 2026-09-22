#!/usr/bin/env bash
# =============================================================================
# download_databases.sh — Download and build databases for ground truth evidence
# =============================================================================
# Usage: sbatch download_databases.sh
# Or:    bash download_databases.sh
#
# Downloads:
#   1. NCBI RefSeq viral proteins → DIAMOND database
#   2. MetaVR v5 (IMG/VR v5) nucleotide sequences → BLAST database
# =============================================================================

#SBATCH --cpus-per-task=4
#SBATCH -t 24:00:00
#SBATCH --mem=128G
#SBATCH -J download_dbs
#SBATCH -o logs/download_databases_%j.out
#SBATCH -e logs/download_databases_%j.err

set -euo pipefail

# --- Configuration -----------------------------------------------------------
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
DB_DIR="${PROJ_DIR}/databases"
REFSEQ_DIR="${DB_DIR}/refseq_viral_prot"
METAVR_DIR="${DB_DIR}/metavr_v5"
NCBI_FTP="https://ftp.ncbi.nlm.nih.gov/refseq/release/viral"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

mkdir -p "${REFSEQ_DIR}" "${METAVR_DIR}" logs

# --- Load modules ------------------------------------------------------------
module load DIAMOND/2.1.11-GCC-13.3.0 2>/dev/null || true

# =============================================================================
# 1. RefSeq Viral Proteins
# =============================================================================
echo "[$(date)] Downloading RefSeq viral protein sequences..."

cd "${REFSEQ_DIR}"

# Download all viral protein FASTA files from RefSeq release
for i in $(seq 1 5); do
    URL="${NCBI_FTP}/viral.${i}.protein.faa.gz"
    OUTFILE="viral.${i}.protein.faa.gz"
    if [[ ! -f "${OUTFILE}" ]]; then
        echo "  Downloading ${OUTFILE}..."
        wget -q --retry-connrefused --waitretry=5 --tries=3 "${URL}" -O "${OUTFILE}" || {
            echo "  WARNING: ${OUTFILE} not found (may not exist for current release). Skipping."
            rm -f "${OUTFILE}"
        }
    else
        echo "  ${OUTFILE} already exists, skipping download."
    fi
done

# Concatenate into single FASTA
echo "[$(date)] Concatenating viral protein FASTAs..."
zcat viral.*.protein.faa.gz > refseq_viral_proteins.faa 2>/dev/null || {
    echo "ERROR: No viral protein files downloaded. Check NCBI FTP availability."
    exit 1
}

NSEQS=$(grep -c "^>" refseq_viral_proteins.faa)
echo "  Total viral protein sequences: ${NSEQS}"

# Build DIAMOND database
echo "[$(date)] Building DIAMOND database..."
diamond makedb \
    --in refseq_viral_proteins.faa \
    --db refseq_viral \
    --threads "${THREADS}"

echo "  DIAMOND database: ${REFSEQ_DIR}/refseq_viral.dmnd"
ls -lh refseq_viral.dmnd

# =============================================================================
# 2. MetaVR v5 (IMG/VR v5) nucleotide sequences
# =============================================================================
METAVR_URL="https://meta-virome.org/Data/Downloads/IMGVR5_UViG.fna.gz"
METAVR_FASTA="${METAVR_DIR}/IMGVR5_UViG.fna"
METAVR_DB="${METAVR_DIR}/metavr_v5"

if [[ ! -f "${METAVR_DB}.ndb" ]] && [[ ! -f "${METAVR_DB}.nsq" ]]; then
    echo "[$(date)] Downloading MetaVR v5 (IMG/VR v5) nucleotide sequences..."
    echo "  Source: ${METAVR_URL}"
    echo "  WARNING: ~77 GB compressed download. This will take a while."

    if [[ ! -f "${METAVR_FASTA}" ]]; then
        METAVR_GZ="${METAVR_DIR}/IMGVR5_UViG.fna.gz"
        MAX_ATTEMPTS=20
        BACKOFF=30  # initial wait (seconds); doubles each attempt

        for attempt in $(seq 1 ${MAX_ATTEMPTS}); do
            echo "[$(date)] Download attempt ${attempt}/${MAX_ATTEMPTS}..."

            # Capture curl exit code separately to prevent set -e from killing the script
            curl_exit=0
            http_code=$(curl -C - -L -o "${METAVR_GZ}" \
                -w '%{http_code}' \
                --retry 3 --retry-delay 15 --retry-max-time 300 \
                --connect-timeout 30 --max-time 0 \
                -H 'User-Agent: Mozilla/5.0 (compatible; viral_bench/1.0; academic research)' \
                "${METAVR_URL}" 2>>"${METAVR_DIR}/download.log") || curl_exit=$?

            echo "  curl exit=${curl_exit}, HTTP ${http_code}"

            # Check file size regardless of exit code — partial downloads are still useful
            actual_size=$(stat -c%s "${METAVR_GZ}" 2>/dev/null || echo 0)

            if [[ "${http_code}" == "200" ]] || [[ "${http_code}" == "206" ]]; then
                if (( actual_size >= 82000000000 )); then
                    echo "[$(date)] Download complete (${actual_size} bytes)."
                    break
                fi
            fi

            echo "  Partial download (${actual_size} bytes), will resume..."

            if (( attempt == MAX_ATTEMPTS )); then
                echo "ERROR: Download failed after ${MAX_ATTEMPTS} attempts." >&2
                exit 1
            fi

            sleep_time=$(( BACKOFF * (2 ** (attempt - 1)) ))
            # Cap at 30 minutes
            (( sleep_time > 1800 )) && sleep_time=1800
            echo "  Sleeping ${sleep_time}s before retry..."
            sleep "${sleep_time}"
        done

        echo "[$(date)] Decompressing..."
        gunzip "${METAVR_GZ}"
    fi

    echo "[$(date)] Building BLAST database..."
    module load BLAST+/2.17.0-gompi-2024a 2>/dev/null || true
    makeblastdb \
        -in "${METAVR_FASTA}" \
        -dbtype nucl \
        -out "${METAVR_DB}" \
        -title "MetaVR v5 (IMG/VR v5)"
else
    echo "  MetaVR v5 BLAST DB already exists: ${METAVR_DB}"
fi

cd "${PROJ_DIR}"

echo ""
echo "[$(date)] Database setup complete."
echo "  RefSeq viral DIAMOND: ${REFSEQ_DIR}/refseq_viral.dmnd"
echo "  MetaVR v5 BLAST DB:   ${METAVR_DB}"
