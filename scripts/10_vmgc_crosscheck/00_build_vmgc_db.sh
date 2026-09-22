#!/usr/bin/env bash
# =============================================================================
# 00_build_vmgc_db.sh — Prepare VMGC reference DBs for the cross-check analysis
# =============================================================================
# Prerequisite for the VMGC cross-check (Plan: VMGC vaginal-catalogue cross-check).
#   1. Extract BOTH VMGC viral FASTAs (full catalogue + vOTU reps).
#   2. Assert genome counts (14,224 all / 4,263 vOTU).
#   3. Build two BLAST nucleotide DBs: vmgc_all (Part 1), vmgc_votu (Part 2).
#   4. Download VMGC_virus.info metadata + assert the schema (incl. required
#      `viral_sequences` column).
#   5. Build + assert a member->vOTU lookup (every all-genome ID maps to exactly
#      one vOTU).
#
# VMGC is used as an ANNOTATION layer only; it is NOT a ground-truth evidence
# line (it was built with DeepVirFinder/VIBRANT/CheckV, all panel tools).
#
# Usage:  sbatch scripts/10_vmgc_crosscheck/00_build_vmgc_db.sh
#    or:  bash   scripts/10_vmgc_crosscheck/00_build_vmgc_db.sh   (it is light)
# =============================================================================

#SBATCH --cpus-per-task=8
#SBATCH -t 02:00:00
#SBATCH --mem=32G
#SBATCH -J vmgc_build
#SBATCH -o logs/vmgc_build_%j.out
#SBATCH -e logs/vmgc_build_%j.err

set -euo pipefail

PROJ="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
VMGC_TARBALL="${VMGC_TARBALL:?set VMGC_TARBALL (VMGC_virus.tar.gz, Zenodo 10.5281/zenodo.10457006)}"
VMGC_DIR="${PROJ}/databases/vmgc_virus"
DATA_DIR="${PROJ}/data/vmgc"
INFO_URL="https://zenodo.org/records/10457006/files/VMGC_virus.info?download=1"
INFO_TSV="${DATA_DIR}/vmgc_virus_info.tsv"
MEMBER_MAP="${DATA_DIR}/member_to_votu.tsv"
ALL_FA="${VMGC_DIR}/VMGC_virus_all.fa"
VOTU_FA="${VMGC_DIR}/VMGC_virus_vOTU.fa"

mkdir -p "${VMGC_DIR}" "${DATA_DIR}" "${PROJ}/logs"
cd "${PROJ}"

module load BLAST+/2.17.0-gompi-2024a 2>/dev/null || module load BLAST+/2.16.0-gompi-2024a 2>/dev/null || true
module load SAMtools/1.21 2>/dev/null || true

echo "================================================================"
echo "VMGC cross-check DB build"
echo "  Tarball: ${VMGC_TARBALL}"
echo "  Out:     ${VMGC_DIR}"
echo "  Started: $(date)"
echo "================================================================"

# --- 1. Extract both FASTAs --------------------------------------------------
if [[ ! -f "${VOTU_FA}" || ! -f "${ALL_FA}" ]]; then
    echo "[$(date)] Extracting VMGC viral FASTAs..."
    tar -xzf "${VMGC_TARBALL}" -C "${VMGC_DIR}" VMGC_virus_vOTU.fa VMGC_virus_all.fa
else
    echo "[$(date)] FASTAs already extracted — skipping."
fi

# --- 2. Index + assert counts ------------------------------------------------
echo "[$(date)] Indexing + counting..."
samtools faidx "${VOTU_FA}"
samtools faidx "${ALL_FA}"
N_VOTU=$(grep -c '^>' "${VOTU_FA}")
N_ALL=$(grep -c '^>' "${ALL_FA}")
echo "  vOTU reps: ${N_VOTU} (expect 4263)"
echo "  all genomes: ${N_ALL} (expect 14224)"
[[ "${N_VOTU}" -eq 4263  ]] || { echo "ERROR: vOTU count ${N_VOTU} != 4263" >&2; exit 1; }
[[ "${N_ALL}"  -eq 14224 ]] || { echo "ERROR: all-genome count ${N_ALL} != 14224" >&2; exit 1; }

# --- 3. Build two BLAST DBs --------------------------------------------------
# Note: -parse_seqids omitted (the working project scripts omit it too); the
# member ID we need is the first token of each header == BLAST outfmt6 sseqid.
echo "[$(date)] Building BLAST DB: vmgc_all (Part 1, full catalogue)..."
makeblastdb -in "${ALL_FA}" -dbtype nucl -out "${VMGC_DIR}/vmgc_all" \
    -title "VMGC viral all-genomes (Huang 2024, 14224)"

echo "[$(date)] Building BLAST DB: vmgc_votu (Part 2, species reps)..."
makeblastdb -in "${VOTU_FA}" -dbtype nucl -out "${VMGC_DIR}/vmgc_votu" \
    -title "VMGC viral vOTU reps (Huang 2024, 4263)"

# --- 4. Download metadata + assert schema ------------------------------------
if [[ ! -s "${INFO_TSV}" ]]; then
    echo "[$(date)] Downloading VMGC_virus.info metadata..."
    curl -L --retry 3 --retry-delay 10 -o "${INFO_TSV}" "${INFO_URL}"
else
    echo "[$(date)] Metadata already present — skipping download."
fi

# --- 5. Build + assert member->vOTU map (also asserts schema) ----------------
echo "[$(date)] Building member->vOTU lookup + asserting schema/coverage..."
python3 - "${INFO_TSV}" "${ALL_FA}.fai" "${MEMBER_MAP}" << 'PYEOF'
import csv, sys
info_tsv, all_fai, out_map = sys.argv[1], sys.argv[2], sys.argv[3]

with open(info_tsv, newline='') as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    cols = reader.fieldnames or []
    # --- schema assertions ---
    required = ['vOTU_ID', 'viral_sequences']
    missing = [c for c in required if c not in cols]
    if missing:
        sys.exit(f"ERROR: VMGC_virus.info missing required column(s): {missing}\n"
                 f"  columns present: {cols}")
    rows = list(reader)

n_rows = len(rows)
print(f"  metadata rows: {n_rows} (expect 4263)")
if n_rows != 4263:
    sys.exit(f"ERROR: metadata row count {n_rows} != 4263")

# explode viral_sequences -> member_to_votu
member_to_votu = {}
dups = 0
for r in rows:
    votu = r['vOTU_ID'].strip()
    members = [m.strip() for m in r['viral_sequences'].split(',') if m.strip()]
    for m in members:
        if m in member_to_votu and member_to_votu[m] != votu:
            dups += 1
        member_to_votu[m] = votu

with open(out_map, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['member_id', 'vOTU_ID'])
    for m, v in sorted(member_to_votu.items()):
        w.writerow([m, v])
print(f"  member IDs mapped: {len(member_to_votu)} (dups across vOTUs: {dups})")

# assert every all-genome FASTA ID maps to exactly one vOTU
fasta_ids = []
with open(all_fai) as fh:
    for line in fh:
        fasta_ids.append(line.split('\t')[0])
orphans = [i for i in fasta_ids if i not in member_to_votu]
print(f"  all-genome FASTA IDs: {len(fasta_ids)} ; orphans (no vOTU): {len(orphans)}")
if orphans:
    ex = ', '.join(orphans[:5])
    sys.exit(f"ERROR: {len(orphans)} all-genome IDs have no vOTU mapping "
             f"(header-ID vs viral_sequences format mismatch). e.g. {ex}")
if dups:
    sys.exit(f"ERROR: {dups} member IDs map to >1 vOTU — not a clean partition")
print("  OK: every all-genome ID maps to exactly one vOTU (14224/14224).")
PYEOF

echo ""
echo "[$(date)] DONE. Built:"
echo "  ${VMGC_DIR}/vmgc_all.*   (Part 1)"
echo "  ${VMGC_DIR}/vmgc_votu.*  (Part 2)"
echo "  ${INFO_TSV}"
echo "  ${MEMBER_MAP}"
