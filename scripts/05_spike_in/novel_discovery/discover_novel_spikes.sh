#!/usr/bin/env bash
# =============================================================================
# discover_novel_spikes.sh — Discover novel phage contigs from internal vaginal
#                            assemblies for training-leakage-free spike-ins
# =============================================================================
# Pipeline: Extract SPAdes assemblies -> filter >= 1500 bp -> CheckV -> DIAMOND
#           -> skani novelty filter -> candidate summary
#
# Uses only non-benchmarked tools (CheckV, DIAMOND, skani) to avoid circularity.
#
# Usage:
#   sbatch discover_novel_spikes.sh [--pilot N]   # N = number of assemblies (default: all)
#   sbatch discover_novel_spikes.sh --pilot 5      # quick pilot on 5 assemblies
# =============================================================================

#SBATCH --cpus-per-task=16
#SBATCH -t 12:00:00
#SBATCH --mem=72G
#SBATCH -J novel_spikes
#SBATCH -o logs/novel_spikes_%j.out
#SBATCH -e logs/novel_spikes_%j.err

set -euo pipefail

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
PILOT_N=0  # 0 = all assemblies
while [[ $# -gt 0 ]]; do
    case "$1" in
        --pilot) PILOT_N="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
TARBALL="${MITCH_ASSEMBLY_TARBALL:?set MITCH_ASSEMBLY_TARBALL (see config/paths.example.sh)}"
OUTDIR="${PROJ_DIR}/data/novel_spike_discovery"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
MIN_CONTIG_LEN=1500

# CheckV setup (same as scripts/04_ground_truth/run_checkv.sh)
CHECKV_DB="${PROJ_DIR}/REFs/checkv-db/checkv-db-v1.5"
CHECKV_SIF="${PROJ_DIR}/singularity-images/checkv-1.0.1--pyhdfd78af_0.sif"
CHECKV_CONTAINER="https://depot.galaxyproject.org/singularity/checkv:1.0.1--pyhdfd78af_0"

# DIAMOND setup (same as scripts/04_ground_truth/run_diamond_blastx.sh)
DIAMOND_DB="${PROJ_DIR}/databases/refseq_viral_prot/refseq_viral.dmnd"

# skani reference DBs for novelty filter
REFSEQ_VIRAL_GENOMES="${PROJ_DIR}/databases/refseq_viral_genomes"  # will build sketch if needed

# Bind the project directory (VBENCH_BIND) into the containers
export APPTAINER_BIND="${VBENCH_BIND:?set VBENCH_BIND (see config/paths.example.sh)},${TMPDIR:-/tmp}"

mkdir -p "${OUTDIR}"/{internal_assemblies,checkv_output,diamond_hits,candidates} logs

echo "============================================================"
echo "[$(date)] Novel Spike-In Discovery Pipeline"
echo "============================================================"
echo "  Tarball:     ${TARBALL}"
echo "  Output:      ${OUTDIR}"
echo "  Threads:     ${THREADS}"
echo "  Min contig:  ${MIN_CONTIG_LEN} bp"
echo "  Pilot mode:  ${PILOT_N:-all assemblies}"
echo ""

# =============================================================================
# Step 1: Extract SPAdes assemblies from tar.gz
# =============================================================================
echo "[$(date)] Step 1: Extracting SPAdes assemblies..."

# Get list of SPAdes assemblies
ASSEMBLY_LIST="${OUTDIR}/spades_assembly_list.txt"
tar tzf "${TARBALL}" | grep "^Assembly/SPAdes-" > "${ASSEMBLY_LIST}"
TOTAL_ASSEMBLIES=$(wc -l < "${ASSEMBLY_LIST}")
echo "  Found ${TOTAL_ASSEMBLIES} SPAdes assemblies in tarball"

# Subset for pilot mode
if [[ ${PILOT_N} -gt 0 ]]; then
    head -n "${PILOT_N}" "${ASSEMBLY_LIST}" > "${ASSEMBLY_LIST}.pilot"
    mv "${ASSEMBLY_LIST}.pilot" "${ASSEMBLY_LIST}"
    echo "  Pilot mode: extracting first ${PILOT_N} assemblies"
fi

N_EXTRACT=$(wc -l < "${ASSEMBLY_LIST}")
echo "  Extracting ${N_EXTRACT} assemblies..."

# Extract the selected assemblies one at a time
# (--files-from is very slow on 30 GB tar.gz since it scans entire archive per file)
echo "  Extracting assemblies one at a time..."
N_DONE=0
while IFS= read -r entry; do
    tar xzf "${TARBALL}" -C "${OUTDIR}/internal_assemblies/" "${entry}" 2>/dev/null && {
        N_DONE=$((N_DONE + 1))
        echo "    [${N_DONE}/${N_EXTRACT}] ${entry}"
    } || {
        echo "    WARNING: Failed to extract ${entry}" >&2
    }
done < "${ASSEMBLY_LIST}"

echo "  Extracted ${N_DONE} assemblies to: ${OUTDIR}/internal_assemblies/Assembly/"

# =============================================================================
# Step 2: Concatenate and filter contigs >= MIN_CONTIG_LEN bp
# =============================================================================
echo ""
echo "[$(date)] Step 2: Concatenating contigs >= ${MIN_CONTIG_LEN} bp..."

CONCAT_FASTA="${OUTDIR}/all_contigs_filtered.fasta"
> "${CONCAT_FASTA}"  # truncate

N_CONTIGS_TOTAL=0
N_CONTIGS_PASS=0

for gz_file in "${OUTDIR}/internal_assemblies/Assembly"/SPAdes-*.contigs.fa.gz; do
    sample_name=$(basename "${gz_file}" .contigs.fa.gz)

    # Decompress and filter by length, adding sample prefix to contig IDs
    zcat "${gz_file}" | awk -v min_len="${MIN_CONTIG_LEN}" -v prefix="${sample_name}" '
    /^>/ {
        if (seq != "" && length(seq) >= min_len) {
            print ">" prefix "__" substr(header, 2)
            print seq
            pass++
        }
        header = $0
        seq = ""
        total++
        next
    }
    { seq = seq $0 }
    END {
        if (seq != "" && length(seq) >= min_len) {
            print ">" prefix "__" substr(header, 2)
            print seq
            pass++
        }
        total++
        # Print stats to stderr
        printf "  %s: %d total contigs, %d >= %d bp\n", prefix, total, pass, min_len > "/dev/stderr"
    }
    ' >> "${CONCAT_FASTA}" 2>&1
done

N_CONTIGS_PASS=$(grep -c "^>" "${CONCAT_FASTA}" || echo 0)
echo "  Total contigs >= ${MIN_CONTIG_LEN} bp: ${N_CONTIGS_PASS}"
echo "  Output: ${CONCAT_FASTA}"

if [[ ${N_CONTIGS_PASS} -eq 0 ]]; then
    echo "ERROR: No contigs >= ${MIN_CONTIG_LEN} bp found. Check assemblies." >&2
    exit 1
fi

# =============================================================================
# Step 3: Run CheckV end_to_end
# =============================================================================
echo ""
echo "[$(date)] Step 3: Running CheckV end_to_end..."

# Pull CheckV container if needed
if [[ ! -f "${CHECKV_SIF}" ]]; then
    echo "  Pulling CheckV container..."
    mkdir -p "$(dirname "${CHECKV_SIF}")"
    singularity pull "${CHECKV_SIF}" "${CHECKV_CONTAINER}"
fi

# Build DIAMOND index if missing
DMND_INDEX="${CHECKV_DB}/genome_db/checkv_reps.dmnd"
if [[ ! -f "${DMND_INDEX}" ]]; then
    echo "  Building DIAMOND index for CheckV..."
    singularity exec \
        "${CHECKV_SIF}" \
        diamond makedb \
            --in "${CHECKV_DB}/genome_db/checkv_reps.faa" \
            --db "${CHECKV_DB}/genome_db/checkv_reps"
fi

CHECKV_OUT="${OUTDIR}/checkv_output"
singularity exec \
    --env TMPDIR="${CHECKV_OUT}/tmp" \
    --pwd "${CHECKV_OUT}" \
    "${CHECKV_SIF}" \
    checkv end_to_end \
        "${CONCAT_FASTA}" \
        "${CHECKV_OUT}" \
        -d "${CHECKV_DB}" \
        -t "${THREADS}"

echo "  CheckV output: ${CHECKV_OUT}/quality_summary.tsv"

# --- Filter for high-quality viral contigs ---
echo ""
echo "[$(date)] Step 3b: Filtering for high-quality viral contigs..."

# CheckV quality_summary.tsv columns:
# contig_id, contig_length, provirus, proviral_length, gene_count, viral_genes,
# host_genes, checkv_quality, miuvig_quality, completeness, completeness_method, ...

HQ_CONTIGS="${OUTDIR}/candidates/hq_viral_contigs.txt"
python3 -c "
import csv, sys
with open('${CHECKV_OUT}/quality_summary.tsv') as f:
    reader = csv.DictReader(f, delimiter='\t')
    hq = []
    for row in reader:
        quality = row.get('checkv_quality', '')
        viral_genes = int(row.get('viral_genes', '0') or '0')
        completeness = float(row.get('completeness', '0') or '0')
        length = int(row.get('contig_length', '0') or '0')
        is_provirus = row.get('provirus', '') == 'Yes'

        # Selection criteria:
        # - High-quality or Complete (CheckV quality tiers)
        # - At least 1 viral gene detected
        # - >= 50% completeness (relaxed from 90% for pilot to see more candidates)
        if quality in ('High-quality', 'Complete') and viral_genes >= 1 and completeness >= 50:
            hq.append({
                'contig_id': row['contig_id'],
                'length': length,
                'quality': quality,
                'completeness': completeness,
                'viral_genes': viral_genes,
                'host_genes': int(row.get('host_genes', '0') or '0'),
                'is_provirus': is_provirus,
            })

    print(f'High-quality viral contigs: {len(hq)}', file=sys.stderr)

    # Write contig IDs
    with open('${HQ_CONTIGS}', 'w') as out:
        for h in hq:
            out.write(h['contig_id'] + '\n')

    # Write summary
    with open('${OUTDIR}/candidates/checkv_hq_summary.tsv', 'w') as out:
        writer = csv.DictWriter(out, fieldnames=['contig_id','length','quality','completeness','viral_genes','host_genes','is_provirus'], delimiter='\t')
        writer.writeheader()
        writer.writerows(hq)

    # Print top candidates
    hq.sort(key=lambda x: (-x['completeness'], -x['length']))
    print(f'\nTop candidates:', file=sys.stderr)
    for h in hq[:20]:
        print(f'  {h[\"contig_id\"]}: {h[\"length\"]:,} bp, {h[\"quality\"]}, {h[\"completeness\"]:.0f}% complete, {h[\"viral_genes\"]} viral genes, provirus={h[\"is_provirus\"]}', file=sys.stderr)
"

N_HQ=$(wc -l < "${HQ_CONTIGS}" || echo 0)
echo "  HQ viral contigs: ${N_HQ}"

if [[ ${N_HQ} -eq 0 ]]; then
    echo "WARNING: No high-quality viral contigs found. Try relaxing filters or processing more assemblies."
    echo "  CheckV quality distribution:"
    awk -F'\t' 'NR>1 {print $8}' "${CHECKV_OUT}/quality_summary.tsv" | sort | uniq -c | sort -rn
    exit 0
fi

# Extract HQ contig sequences (use proviruses.fna for trimmed proviral regions)
HQ_FASTA="${OUTDIR}/candidates/hq_viral_contigs.fasta"
python3 -c "
import sys

# Load HQ contig IDs
with open('${HQ_CONTIGS}') as f:
    hq_ids = set(l.strip() for l in f)

# First collect from proviruses.fna (trimmed proviral regions)
provirus_ids = set()
provirus_file = '${CHECKV_OUT}/proviruses.fna'
try:
    with open(provirus_file) as f:
        for line in f:
            if line.startswith('>'):
                cid = line[1:].strip().split()[0]
                if cid in hq_ids:
                    provirus_ids.add(cid)
except FileNotFoundError:
    pass

# Collect sequences
# Proviruses from proviruses.fna, others from original FASTA
found = set()
with open('${HQ_FASTA}', 'w') as out:
    # First: proviruses (trimmed)
    if provirus_ids:
        with open(provirus_file) as f:
            writing = False
            for line in f:
                if line.startswith('>'):
                    cid = line[1:].strip().split()[0]
                    writing = cid in hq_ids
                    if writing:
                        found.add(cid)
                if writing:
                    out.write(line)

    # Second: non-proviral from original FASTA
    remaining = hq_ids - found
    if remaining:
        with open('${CONCAT_FASTA}') as f:
            writing = False
            for line in f:
                if line.startswith('>'):
                    cid = line[1:].strip().split()[0]
                    writing = cid in remaining
                    if writing:
                        found.add(cid)
                if writing:
                    out.write(line)

print(f'Extracted {len(found)}/{len(hq_ids)} HQ sequences ({len(provirus_ids)} proviruses trimmed)', file=sys.stderr)
"

# =============================================================================
# Step 4: DIAMOND BLASTx confirmation
# =============================================================================
echo ""
echo "[$(date)] Step 4: DIAMOND BLASTx vs RefSeq viral proteins..."

if [[ ! -f "${DIAMOND_DB}" ]]; then
    echo "ERROR: DIAMOND database not found: ${DIAMOND_DB}" >&2
    echo "  Run scripts/04_ground_truth/run_diamond_blastx.sh first to set up the DB"
    exit 1
fi

DIAMOND_OUT="${OUTDIR}/diamond_hits/hq_viral_blastx.tsv"
module load DIAMOND/2.1.11-GCC-13.3.0 2>/dev/null || true

diamond blastx \
    --query "${HQ_FASTA}" \
    --db "${DIAMOND_DB}" \
    --out "${DIAMOND_OUT}" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
    --evalue 1e-10 \
    --max-target-seqs 5 \
    --threads "${THREADS}" \
    --block-size 4 \
    --index-chunks 1 \
    --quiet

# Filter: contigs with at least 1 significant hit
DIAMOND_CONFIRMED="${OUTDIR}/candidates/diamond_confirmed.txt"
awk -F'\t' '{print $1}' "${DIAMOND_OUT}" | sort -u > "${DIAMOND_CONFIRMED}"
N_CONFIRMED=$(wc -l < "${DIAMOND_CONFIRMED}")
echo "  Contigs with viral protein hits: ${N_CONFIRMED}/${N_HQ}"

# =============================================================================
# Step 5: skani novelty filter vs public databases
# =============================================================================
echo ""
echo "[$(date)] Step 5: skani novelty filter..."

module load skani/0.3.1 2>/dev/null || true

# Extract only DIAMOND-confirmed contigs for skani
CONFIRMED_FASTA="${OUTDIR}/candidates/diamond_confirmed.fasta"
python3 -c "
with open('${DIAMOND_CONFIRMED}') as f:
    ids = set(l.strip() for l in f)
with open('${HQ_FASTA}') as fin, open('${CONFIRMED_FASTA}', 'w') as fout:
    writing = False
    for line in fin:
        if line.startswith('>'):
            cid = line[1:].strip().split()[0]
            writing = cid in ids
        if writing:
            fout.write(line)
"

# Run skani dist against RefSeq viral genomes
# If we don't have a local RefSeq viral genome collection, use the DIAMOND DB's
# source FASTAs or skip this step
SKANI_OUT="${OUTDIR}/candidates/skani_dist.tsv"

# Check if we have reference viral genomes for skani comparison
# For now, compare against the spike-in genomes themselves as a sanity check
# + any available viral reference genomes
SPIKE_IN_DIR="${PROJ_DIR}/data/spike_in"
SPIKE_ALL_VIRAL="${SPIKE_IN_DIR}/all_viral_genomes.fasta"
if [[ -f "${SPIKE_ALL_VIRAL}" ]]; then
    echo "  Running skani dist against spike-in reference genomes (sanity check)..."

    # Use pre-concatenated all_viral_genomes.fasta + all_negative_controls.fasta
    SPIKE_REF="${OUTDIR}/candidates/spike_in_reference.fasta"
    cat "${SPIKE_ALL_VIRAL}" "${SPIKE_IN_DIR}/all_negative_controls.fasta" 2>/dev/null > "${SPIKE_REF}" || true

    if [[ -s "${SPIKE_REF}" ]]; then
        skani dist \
            --ql "${CONFIRMED_FASTA}" \
            --rl "${SPIKE_REF}" \
            -o "${SKANI_OUT}" \
            -t "${THREADS}" \
            --ci \
            -n 5 \
            2>/dev/null || true

        if [[ -s "${SKANI_OUT}" ]]; then
            echo "  Matches to spike-in genomes (should be rejected as non-novel):"
            awk -F'\t' 'NR>1 && $3 >= 95 {printf "    %s -> %s (ANI=%.1f%%, AF=%.1f%%)\n", $1, $2, $3, $4}' "${SKANI_OUT}"
        fi
    fi
fi

# For full novelty check against all public DBs, we would need:
# - RefSeq viral genomes (all)
# - IMG/VR v4 / MetaVR v5 UViG sequences
# This requires large downloads. For the pilot, we flag contigs that DON'T match
# our spike-in genomes as "potentially novel" and note that full novelty verification
# is needed before final panel selection.

echo ""
echo "  NOTE: Full novelty verification against RefSeq/IMG-VR requires additional"
echo "  reference downloads. For pilot, candidates not matching spike-in genomes"
echo "  are flagged as potentially novel."

# =============================================================================
# Step 6: Generate candidate summary
# =============================================================================
echo ""
echo "[$(date)] Step 6: Generating candidate summary..."

python3 -c "
import csv, sys

# Load CheckV quality data
checkv_data = {}
with open('${CHECKV_OUT}/quality_summary.tsv') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        checkv_data[row['contig_id']] = row

# Load DIAMOND-confirmed IDs
with open('${DIAMOND_CONFIRMED}') as f:
    confirmed_ids = set(l.strip() for l in f)

# Load skani matches (contigs matching known genomes at >= 95% ANI)
known_matches = {}
try:
    with open('${SKANI_OUT}') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ani = float(row.get('ANI', 0) or 0)
            af_q = float(row.get('Align_fraction_query', 0) or row.get('AF_query', 0) or 0)
            if ani >= 95:
                qname = row.get('Query_file', row.get('Ref_file', ''))
                known_matches[qname] = {'ani': ani, 'af': af_q, 'ref': row.get('Ref_file', '')}
except (FileNotFoundError, KeyError):
    pass

# Build candidate summary
candidates = []
for cid in sorted(confirmed_ids):
    cv = checkv_data.get(cid, {})
    is_known = cid in known_matches
    candidates.append({
        'contig_id': cid,
        'contig_length': cv.get('contig_length', ''),
        'checkv_quality': cv.get('checkv_quality', ''),
        'completeness': cv.get('completeness', ''),
        'viral_genes': cv.get('viral_genes', ''),
        'host_genes': cv.get('host_genes', ''),
        'provirus': cv.get('provirus', ''),
        'taxonomy': cv.get('taxonomy', ''),
        'diamond_confirmed': 'yes',
        'matches_known_genome': 'yes' if is_known else 'no',
        'novelty_status': 'known' if is_known else 'potentially_novel',
    })

# Write summary
outpath = '${OUTDIR}/candidates/candidate_summary.tsv'
with open(outpath, 'w') as f:
    fieldnames = ['contig_id', 'contig_length', 'checkv_quality', 'completeness',
                  'viral_genes', 'host_genes', 'provirus', 'taxonomy',
                  'diamond_confirmed', 'matches_known_genome', 'novelty_status']
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(candidates)

# Print summary
n_novel = sum(1 for c in candidates if c['novelty_status'] == 'potentially_novel')
n_known = sum(1 for c in candidates if c['novelty_status'] == 'known')
print(f'')
print(f'Candidate summary: {len(candidates)} total')
print(f'  Potentially novel:  {n_novel}')
print(f'  Matches known:      {n_known}')
print(f'')
print(f'Top potentially novel candidates (by length):')
novel = [c for c in candidates if c['novelty_status'] == 'potentially_novel']
novel.sort(key=lambda x: -int(x['contig_length'] or 0))
for c in novel[:15]:
    print(f\"  {c['contig_id']}: {int(c['contig_length']):,} bp, {c['checkv_quality']}, {c['completeness']}% complete, {c['viral_genes']} viral genes, taxonomy={c['taxonomy']}\")
print(f'')
print(f'Output: {outpath}')
"

echo ""
echo "============================================================"
echo "[$(date)] Pipeline complete!"
echo "============================================================"
echo ""
echo "Next steps:"
echo "  1. Review: ${OUTDIR}/candidates/candidate_summary.tsv"
echo "  2. Select panel: python select_novel_spikes.py --candidates candidate_summary.tsv"
echo "  3. If too few candidates, rerun with more assemblies (--pilot N or no --pilot for all)"
echo ""
