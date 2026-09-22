#!/usr/bin/env bash
# =============================================================================
# download_spike_in_genomes.sh — Download Spike-In Genomes from NCBI
# =============================================================================
# Downloads all viral and bacterial genomes specified in spike_in_genomes.tsv
# using NCBI datasets CLI (for assemblies) and efetch (for nucleotide accessions).
#
# Output: data/spike_in/{category}/{accession}.fasta
#         data/spike_in/all_viral_genomes.fasta   (concatenated viral panel)
#         data/spike_in/all_negative_controls.fasta (concatenated negatives)
#
# Usage: bash download_spike_in_genomes.sh
#        sbatch download_spike_in_genomes.sh   (for SLURM)
# =============================================================================

#SBATCH --cpus-per-task=2
#SBATCH -t 02:00:00
#SBATCH --mem=16G
#SBATCH -J download_spikein
#SBATCH -o logs/download_spikein_%j.out
#SBATCH -e logs/download_spikein_%j.err

set -euo pipefail

# --- Configuration -----------------------------------------------------------
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
SCRIPT_DIR="${PROJ_DIR}/scripts/05_spike_in"
MANIFEST="${SCRIPT_DIR}/spike_in_genomes.tsv"
OUTDIR="${PROJ_DIR}/data/spike_in"
METAVR_FASTA="${PROJ_DIR}/databases/metavr_v5/IMGVR5_UViG.fna"

mkdir -p "${OUTDIR}" logs

# --- Validate manifest -------------------------------------------------------
if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: Manifest not found: ${MANIFEST}" >&2
    exit 1
fi

# --- Check tools -------------------------------------------------------------
# Try to load Entrez Direct (efetch) and NCBI datasets
# the HPC cluster: no bioinfo-tools needed, load modules directly
module load entrez-direct 2>/dev/null || true

if ! command -v efetch &>/dev/null; then
    echo "WARNING: efetch not found. Trying curl fallback for nucleotide accessions."
    USE_EFETCH=false
else
    USE_EFETCH=true
fi

# =============================================================================
# Download each genome
# =============================================================================
echo "[$(date)] Downloading spike-in genomes..."
echo "  Manifest: ${MANIFEST}"
echo "  Output:   ${OUTDIR}"
echo ""

TOTAL=0
SUCCESS=0
FAILED=0

while IFS=$'\t' read -r accession organism category genome_size genome_type source notes; do
    # Skip comments and header
    [[ "${accession}" =~ ^#.*$ ]] && continue
    [[ -z "${accession}" ]] && continue

    TOTAL=$((TOTAL + 1))
    CATDIR="${OUTDIR}/${category}"
    mkdir -p "${CATDIR}"

    OUTFILE="${CATDIR}/${accession}.fasta"

    if [[ -f "${OUTFILE}" ]] && [[ -s "${OUTFILE}" ]]; then
        echo "  [SKIP] ${accession} (${organism}) — already exists"
        SUCCESS=$((SUCCESS + 1))
        continue
    fi

    echo "  [GET]  ${accession} (${organism}) → ${category}/"

    if [[ "${accession}" == IMGVR_* ]]; then
        # MetaVR v5 UViG — extract from local FASTA (254 GB)
        if [[ ! -f "${METAVR_FASTA}" ]]; then
            echo "         ✗ FAILED — MetaVR FASTA not found: ${METAVR_FASTA}"
            FAILED=$((FAILED + 1))
            continue
        fi

        # Index if needed (first run only; takes ~15 min for 254 GB)
        if [[ ! -f "${METAVR_FASTA}.fai" ]]; then
            echo "         Indexing MetaVR FASTA (first run only, ~15 min)..."
            module load samtools 2>/dev/null || module load SAMtools/1.22-GCC-13.3.0 2>/dev/null || true
            samtools faidx "${METAVR_FASTA}"
        fi

        module load samtools 2>/dev/null || module load SAMtools/1.22-GCC-13.3.0 2>/dev/null || true

        # FASTA headers include pipe-separated metadata after the UViG ID
        # (e.g., IMGVR_UViG_xxx|taxon_oid|contig_id). Look up the full header
        # from the .fai index so samtools faidx extracts the full sequence.
        FULL_HEADER=$(grep "^${accession}" "${METAVR_FASTA}.fai" | cut -f1)
        if [[ -z "${FULL_HEADER}" ]]; then
            echo "         ✗ FAILED — ${accession} not found in MetaVR FASTA index"
            FAILED=$((FAILED + 1))
            continue
        fi
        samtools faidx "${METAVR_FASTA}" "${FULL_HEADER}" > "${OUTFILE}" 2>/dev/null

        # Rename FASTA header to match convention: >ACCESSION Organism, complete genome
        sed -i "1s|^>.*|>${accession} ${organism}, complete genome|" "${OUTFILE}"

    elif [[ "${accession}" == GCF_* ]] || [[ "${accession}" == GCA_* ]]; then
        # Assembly accession — use NCBI datasets CLI or direct FTP
        # Try datasets first, fall back to curl
        if command -v datasets &>/dev/null; then
            datasets download genome accession "${accession}" \
                --include genome \
                --filename "${CATDIR}/${accession}.zip" 2>/dev/null

            if [[ -f "${CATDIR}/${accession}.zip" ]]; then
                unzip -o -j "${CATDIR}/${accession}.zip" "ncbi_dataset/data/${accession}/*.fna" \
                    -d "${CATDIR}/" 2>/dev/null
                # Rename to standard name
                mv "${CATDIR}"/*.fna "${OUTFILE}" 2>/dev/null || true
                rm -f "${CATDIR}/${accession}.zip"
            fi
        else
            # Fallback: direct NCBI FTP download
            # Construct FTP path: GCF_NNNNNNNNN.V → GCF/NNN/NNN/NNN/
            ACC_NUM=$(echo "${accession}" | sed 's/GCF_//' | sed 's/\..*//')
            P1=${ACC_NUM:0:3}
            P2=${ACC_NUM:3:3}
            P3=${ACC_NUM:6:3}
            # List the FTP directory to find the exact assembly folder name
            FTP_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${P1}/${P2}/${P3}"
            ASM_DIR=$(curl -sS "${FTP_BASE}/" 2>/dev/null | grep -oP "href=\"${accession}[^\"]*\"" | head -1 | tr -d '"' | sed 's/href=//')
            if [[ -n "${ASM_DIR}" ]]; then
                FTP_URL="${FTP_BASE}/${ASM_DIR}${ASM_DIR%/}_genomic.fna.gz"
                curl -sS -o "${OUTFILE}.gz" "${FTP_URL}" 2>/dev/null && \
                    gunzip -f "${OUTFILE}.gz" 2>/dev/null || true
            else
                echo "         ✗ FAILED — could not resolve FTP path for ${accession}"
            fi
        fi
    else
        # Nucleotide accession — use efetch or curl
        if [[ "${USE_EFETCH}" == true ]]; then
            efetch -db nucleotide -id "${accession}" -format fasta > "${OUTFILE}" 2>/dev/null
        else
            # Fallback: NCBI E-utilities API
            curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${accession}&rettype=fasta&retmode=text" \
                > "${OUTFILE}" 2>/dev/null
            # Rate limit: NCBI allows 3 requests/second without API key
            sleep 0.4
        fi
    fi

    # Validate download
    if [[ -f "${OUTFILE}" ]] && [[ -s "${OUTFILE}" ]] && grep -q "^>" "${OUTFILE}"; then
        SIZE=$(wc -c < "${OUTFILE}")
        NSEQS=$(grep -c "^>" "${OUTFILE}")
        echo "         ✓ ${SIZE} bytes, ${NSEQS} sequence(s)"
        SUCCESS=$((SUCCESS + 1))
    else
        echo "         ✗ FAILED — empty or invalid FASTA"
        rm -f "${OUTFILE}"
        FAILED=$((FAILED + 1))
    fi

done < "${MANIFEST}"

# =============================================================================
# Concatenate genomes
# =============================================================================
echo ""
echo "[$(date)] Creating concatenated files..."

# All viral genomes (everything except negative_control)
cat "${OUTDIR}"/lactobacillus_phage/*.fasta \
    "${OUTDIR}"/gardnerella_phage/*.fasta \
    "${OUTDIR}"/megasphaera_phage/*.fasta \
    "${OUTDIR}"/fannyhessea_phage/*.fasta \
    "${OUTDIR}"/sneathia_phage/*.fasta \
    "${OUTDIR}"/hpv/*.fasta \
    "${OUTDIR}"/anellovirus/*.fasta \
    "${OUTDIR}"/herpesvirus/*.fasta \
    2>/dev/null > "${OUTDIR}/all_viral_genomes.fasta" || true

# All negative controls
cat "${OUTDIR}"/negative_control/*.fasta \
    2>/dev/null > "${OUTDIR}/all_negative_controls.fasta" || true

NVIRAL=$(grep -c "^>" "${OUTDIR}/all_viral_genomes.fasta" 2>/dev/null || echo "0")
NNEG=$(grep -c "^>" "${OUTDIR}/all_negative_controls.fasta" 2>/dev/null || echo "0")

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "[$(date)] Download complete."
echo "  Total genomes:    ${TOTAL}"
echo "  Successful:       ${SUCCESS}"
echo "  Failed:           ${FAILED}"
echo ""
echo "  Viral panel:      ${NVIRAL} sequences → ${OUTDIR}/all_viral_genomes.fasta"
echo "  Negative controls: ${NNEG} sequences → ${OUTDIR}/all_negative_controls.fasta"
echo ""
echo "  Per-category directories:"
for d in "${OUTDIR}"/*/; do
    [[ -d "$d" ]] || continue
    N=$(find "$d" -name "*.fasta" | wc -l)
    echo "    $(basename "$d"): ${N} genomes"
done
echo ""
echo "  Next steps:"
echo "    1. Verify accessions: check_training_overlap.py"
echo "    2. Simulate reads:    simulate_spike_in_reads.sh"
echo "    3. Fragment genomes:  fragment_genomes.py"
