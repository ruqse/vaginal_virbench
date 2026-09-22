#!/usr/bin/env bash
set -euo pipefail

# ── Setup geNomad container and database ──
# RUN ON LOGIN NODE (compute nodes lack internet access):
#   bash scripts/06_tool_execution/setup_genomad.sh

PROJ="${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}"
CONTAINER_DIR="${PROJ}/containers/genomad"
DB_DIR="${PROJ}/databases/genomad_db"

mkdir -p "${CONTAINER_DIR}" "${DB_DIR}"

# ── 1. Pull container (skip if already done) ──
if [[ ! -f "${CONTAINER_DIR}/genomad.sif" ]]; then
    echo "[$(date)] Pulling geNomad container..."
    #module load apptainer/1.3.4
    apptainer pull \
        "${CONTAINER_DIR}/genomad.sif" \
        "docker://quay.io/biocontainers/genomad:1.11.2--pyhdfd78af_0"
    echo "[$(date)] Container pulled."
else
    echo "[$(date)] Container already exists: ${CONTAINER_DIR}/genomad.sif"
fi

# ── 2. Download database (~3 GB) ──
# Uses wget instead of `genomad download-database` to bypass
# SSL CERTIFICATE_VERIFY_FAILED inside the BioContainers image.
if [[ ! -f "${DB_DIR}/genomad_db/names.dmp" ]]; then
    DB_VERSION="1.9"
    DB_URL="https://portal.nersc.gov/genomad/__data__/genomad_db_v${DB_VERSION}.tar.gz"
    echo "[$(date)] Downloading geNomad database v${DB_VERSION} via wget..."
    wget -q "${DB_URL}" -O "${DB_DIR}/genomad_db_v${DB_VERSION}.tar.gz"
    tar xzf "${DB_DIR}/genomad_db_v${DB_VERSION}.tar.gz" -C "${DB_DIR}"
    rm "${DB_DIR}/genomad_db_v${DB_VERSION}.tar.gz"
    echo "[$(date)] Database downloaded to: ${DB_DIR}/genomad_db/"
else
    echo "[$(date)] Database already exists: ${DB_DIR}/genomad_db/"
fi

echo "[$(date)] geNomad setup complete."
echo "  Container: ${CONTAINER_DIR}/genomad.sif"
echo "  Database:  ${DB_DIR}/genomad_db/"
