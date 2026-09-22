#!/usr/bin/env bash
#SBATCH --cpus-per-task=16
#SBATCH -t 04:00:00
#SBATCH --mem=64G
#SBATCH -J merge_char
#SBATCH -o logs/merge_characterize_%j.out
#SBATCH -e logs/merge_characterize_%j.err

# =============================================================================
# Step 3: Merge CheckV results + DIAMOND + BLASTn MetaVR v5 + candidate selection
# =============================================================================
set -euo pipefail

PROJ="${PROJ:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
OUTDIR="${OUTDIR:-${PROJ}/data/novel_spike_discovery/full_scale}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

CHECKV_DB="${PROJ}/REFs/checkv-db/checkv-db-v1.5"
CHECKV_SIF="${PROJ}/singularity-images/checkv-1.0.1--pyhdfd78af_0.sif"
DIAMOND_DB="${PROJ}/databases/refseq_viral_prot/refseq_viral.dmnd"
METAVR_DB="${PROJ}/databases/metavr_v5/metavr_v5"

export APPTAINER_BIND="${VBENCH_BIND:?set VBENCH_BIND (see config/paths.example.sh)},${SNIC_TMP:-/tmp}"

echo "============================================================"
echo "[$(date)] Step 3: Merge + Characterize"
echo "============================================================"

# --- 3a: Merge all CheckV quality_summary.tsv ---
echo ""
echo "[$(date)] Merging CheckV results..."
MERGED="${OUTDIR}/merged_quality_summary.tsv"
FIRST=true
for qsv in "${OUTDIR}"/checkv/*/quality_summary.tsv; do
    [[ ! -s "${qsv}" ]] && continue
    if ${FIRST}; then
        cat "${qsv}" > "${MERGED}"
        FIRST=false
    else
        tail -n +2 "${qsv}" >> "${MERGED}"
    fi
done

N_TOTAL=$(( $(wc -l < "${MERGED}") - 1 ))
echo "  Merged: ${N_TOTAL} contigs from $(ls "${OUTDIR}"/checkv/*/quality_summary.tsv | wc -l) samples"

# --- 3b: Filter HQ/Complete candidates ---
echo ""
echo "[$(date)] Filtering high-quality viral candidates..."

python3 -c "
import csv, sys

with open('${MERGED}') as f:
    reader = csv.DictReader(f, delimiter='\t')
    hq = []
    for row in reader:
        quality = row.get('checkv_quality', '')
        viral_genes = int(row.get('viral_genes', '0') or '0')
        host_genes = int(row.get('host_genes', '0') or '0')
        comp_str = row.get('completeness', '0') or '0'
        completeness = float(comp_str) if comp_str != 'NA' else 0.0
        length = int(row.get('contig_length', '0') or '0')
        provirus = row.get('provirus', '') == 'Yes'

        # Selection: Complete/HQ, viral_genes >= 5, length >= 10kb
        if quality in ('High-quality', 'Complete') and viral_genes >= 5 and length >= 10000:
            hq.append(row)

print(f'  HQ viral candidates: {len(hq)}')

# Write candidate IDs
with open('${OUTDIR}/candidates/hq_contig_ids.txt', 'w') as f:
    for h in hq:
        f.write(h['contig_id'] + '\n')

# Write summary
with open('${OUTDIR}/candidates/checkv_hq_summary.tsv', 'w') as f:
    writer = csv.DictWriter(f, fieldnames=reader.fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(hq)

# Quality distribution
from collections import Counter
quals = Counter(h['checkv_quality'] for h in hq)
for q, n in sorted(quals.items()):
    print(f'    {q}: {n}')
"

N_HQ=$(wc -l < "${OUTDIR}/candidates/hq_contig_ids.txt")
echo "  HQ candidates: ${N_HQ}"

if [[ ${N_HQ} -eq 0 ]]; then
    echo "WARNING: No HQ candidates found."
    exit 0
fi

# --- 3c: Extract HQ sequences ---
echo ""
echo "[$(date)] Extracting HQ candidate sequences..."

HQ_FASTA="${OUTDIR}/candidates/hq_viral_contigs.fasta"
python3 -c "
import sys

with open('${OUTDIR}/candidates/hq_contig_ids.txt') as f:
    ids = set(l.strip() for l in f)

# Check proviruses first, then filtered FASTAs
found = set()
with open('${HQ_FASTA}', 'w') as out:
    # Proviruses (trimmed regions)
    import glob
    for prov_file in sorted(glob.glob('${OUTDIR}/checkv/*/proviruses.fna')):
        writing = False
        for line in open(prov_file):
            if line.startswith('>'):
                cid = line[1:].strip().split()[0]
                writing = cid in ids and cid not in found
                if writing:
                    found.add(cid)
            if writing:
                out.write(line)

    # Non-proviral from filtered FASTAs
    remaining = ids - found
    if remaining:
        for filt_file in sorted(glob.glob('${OUTDIR}/filtered/*.fasta')):
            writing = False
            for line in open(filt_file):
                if line.startswith('>'):
                    cid = line[1:].strip().split()[0]
                    writing = cid in remaining and cid not in found
                    if writing:
                        found.add(cid)
                if writing:
                    out.write(line)

print(f'  Extracted {len(found)}/{len(ids)} sequences')
"

# --- 3d: DIAMOND BLASTx ---
echo ""
echo "[$(date)] DIAMOND BLASTx vs RefSeq viral proteins..."
module load DIAMOND/2.1.11-GCC-13.3.0 2>/dev/null || true

diamond blastx \
    --query "${HQ_FASTA}" \
    --db "${DIAMOND_DB}" \
    --out "${OUTDIR}/candidates/diamond_blastx.tsv" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
    --evalue 1e-10 --max-target-seqs 5 --threads "${THREADS}" \
    --block-size 4 --index-chunks 1 --quiet

N_DIAMOND=$(awk -F'\t' '{print $1}' "${OUTDIR}/candidates/diamond_blastx.tsv" | sort -u | wc -l)
echo "  Contigs with viral protein hits: ${N_DIAMOND}"

# --- 3e: BLASTn against MetaVR v5 (95% ANI filter) ---
echo ""
echo "[$(date)] BLASTn vs MetaVR v5 (novelty filter)..."
module load BLAST+/2.17.0-gompi-2024a 2>/dev/null || true

blastn \
    -query "${HQ_FASTA}" \
    -db "${METAVR_DB}" \
    -outfmt "6 qseqid sseqid pident length qlen slen qcovs evalue" \
    -evalue 1e-5 -max_target_seqs 3 -num_threads "${THREADS}" \
    -out "${OUTDIR}/candidates/blastn_metavr_v5.tsv"

echo "  BLASTn hits: $(wc -l < "${OUTDIR}/candidates/blastn_metavr_v5.tsv")"

# --- 3f: Generate final candidate summary with novelty status ---
echo ""
echo "[$(date)] Generating candidate summary..."

python3 -c "
import csv, sys
from collections import defaultdict

# Load CheckV HQ summary
checkv = {}
with open('${OUTDIR}/candidates/checkv_hq_summary.tsv') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        checkv[row['contig_id']] = row

# Load DIAMOND confirmed IDs + best hit
diamond_hits = {}
with open('${OUTDIR}/candidates/diamond_blastx.tsv') as f:
    for line in f:
        fields = line.strip().split('\t')
        qid = fields[0]
        if qid not in diamond_hits:
            diamond_hits[qid] = fields[6] if len(fields) > 6 else fields[1]  # stitle or sseqid

# Load BLASTn MetaVR v5: best hit per query
metavr_best = {}
with open('${OUTDIR}/candidates/blastn_metavr_v5.tsv') as f:
    for line in f:
        fields = line.strip().split('\t')
        qid, sid, pident, alen = fields[0], fields[1], float(fields[2]), int(fields[3])
        qcovs = float(fields[6]) if len(fields) > 6 else 0
        if qid not in metavr_best or pident > metavr_best[qid]['pident']:
            metavr_best[qid] = {'subject': sid, 'pident': pident, 'qcovs': qcovs}

# Build summary
candidates = []
for cid in sorted(checkv.keys()):
    cv = checkv[cid]
    mv = metavr_best.get(cid, {})
    mv_pident = mv.get('pident', 0)
    mv_qcovs = mv.get('qcovs', 0)

    # Novelty: <95% pident to MetaVR = novel species
    if mv_pident < 95:
        novelty = 'novel_below_95pct'
    elif mv_pident >= 99:
        novelty = 'same_species_in_metavr'
    else:
        novelty = 'close_relative_in_metavr'

    candidates.append({
        'contig_id': cid,
        'length': cv.get('contig_length', ''),
        'checkv_quality': cv.get('checkv_quality', ''),
        'completeness': cv.get('completeness', ''),
        'viral_genes': cv.get('viral_genes', ''),
        'host_genes': cv.get('host_genes', ''),
        'provirus': cv.get('provirus', ''),
        'diamond_hit': diamond_hits.get(cid, 'none'),
        'metavr_best_pident': f'{mv_pident:.1f}' if mv_pident > 0 else 'no_hit',
        'metavr_best_qcovs': f'{mv_qcovs:.0f}' if mv_qcovs > 0 else '',
        'metavr_best_subject': mv.get('subject', ''),
        'novelty': novelty,
    })

# Write
outpath = '${OUTDIR}/candidates/candidate_summary.tsv'
with open(outpath, 'w') as f:
    writer = csv.DictWriter(f, fieldnames=candidates[0].keys(), delimiter='\t')
    writer.writeheader()
    writer.writerows(candidates)

# Summary
from collections import Counter
novelty_counts = Counter(c['novelty'] for c in candidates)
print(f'')
print(f'Candidate summary: {len(candidates)} total')
for status, n in sorted(novelty_counts.items()):
    print(f'  {status}: {n}')
print(f'')

novel = [c for c in candidates if c['novelty'] == 'novel_below_95pct']
print(f'Truly novel (<95% ANI to MetaVR v5): {len(novel)}')
for c in sorted(novel, key=lambda x: -int(x['length']))[:20]:
    print(f'  {c[\"contig_id\"]}: {int(c[\"length\"]):,} bp, {c[\"checkv_quality\"]}, {c[\"completeness\"]}% complete, {c[\"viral_genes\"]} viral genes')
    print(f'    DIAMOND: {c[\"diamond_hit\"][:80]}')
    print(f'    MetaVR: {c[\"metavr_best_pident\"]}% id, {c[\"metavr_best_qcovs\"]}% qcov')
print(f'')
print(f'Output: {outpath}')
"

echo ""
echo "============================================================"
echo "[$(date)] Full-scale discovery complete!"
echo "============================================================"
