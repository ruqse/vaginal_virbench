#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# BLAST the 48 Track C "dark matter" contigs against the 15 novel-panel phage
# genomes (including the 2 "dark phages" with zero MetaVR v5 match).
#
# Purpose: test whether any Track C dark-matter contigs match genomes from the
# novel panel used in Track A. A hit would be a gold-standard positive control
# for Track C: the R > 10 enrichment filter recovers novel material that the
# database-based evidence lines (E1a/E1b/CheckV/Kraken2) cannot detect.
#
# Inputs:
#   - 48 dark-matter contig IDs: rca_dark_matter_contig_ids.txt
#     (patient-prefixed, live in all_master_contigs.fasta)
#   - 15 novel-panel genomes: data/novel_spike_discovery/full_scale/selected/
#     novel_viral_genomes.fasta (13 has_relatives + 2 dark phages)
#
# Outputs (to results/test_real/coassembly/dark_matter_novel_blast/):
#   - dark_matter_contigs.fasta       (48 extracted sequences)
#   - novel_panel.db.{nhr,nin,nsq,...}
#   - dark_matter_vs_novel_panel.tsv  (BLAST hits, tabular)
#   - hits_summary.tsv                (hits at >= 80% identity, >= 50% qcov)
#   - README.md                       (interpretation)
#
# Requires: BLAST+ 2.17.0 (or 2.16.0), samtools (for fasta indexing).
# -----------------------------------------------------------------------------

set -euo pipefail

# --- Paths --------------------------------------------------------------------
PROJ=${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}
MASTER_FASTA="${PROJ}/results/test_real/coassembly/all_master_contigs.fasta"
DARK_IDS="${PROJ}/results/test_real/coassembly/rca_dark_matter_contig_ids.txt"
NOVEL_FASTA="${PROJ}/data/novel_spike_discovery/full_scale/selected/novel_viral_genomes.fasta"
NOVEL_ID_MAP="${PROJ}/data/novel_spike_discovery/full_scale/selected/id_mapping.tsv"
OUTDIR="${PROJ}/results/test_real/coassembly/dark_matter_novel_blast"
mkdir -p "${OUTDIR}"

# --- Module environment -------------------------------------------------------
# the HPC cluster: module load bioinfo-tools does NOT work; load directly.
module load BLAST+/2.17.0-gompi-2024a 2>/dev/null || module load BLAST+/2.16.0-gompi-2024a
module load SAMtools/1.21 2>/dev/null || true

# --- 1. Extract 48 dark-matter contigs ---------------------------------------
echo "[$(date)] Extracting dark-matter contigs from master assembly..."
DARK_FASTA="${OUTDIR}/dark_matter_contigs.fasta"

# Index master fasta if needed (samtools faidx is idempotent)
[[ -f "${MASTER_FASTA}.fai" ]] || samtools faidx "${MASTER_FASTA}"

# Use samtools faidx --region-file for bulk extraction (fast, ordered)
# Strip any trailing whitespace from IDs just in case
awk 'NF' "${DARK_IDS}" > "${OUTDIR}/.dark_ids.txt"
samtools faidx "${MASTER_FASTA}" -r "${OUTDIR}/.dark_ids.txt" > "${DARK_FASTA}"

n_extracted=$(grep -c "^>" "${DARK_FASTA}")
n_expected=$(wc -l < "${OUTDIR}/.dark_ids.txt")
echo "  Extracted ${n_extracted}/${n_expected} dark-matter contigs"
if [[ "${n_extracted}" -ne "${n_expected}" ]]; then
    echo "  WARNING: count mismatch — some IDs not found in ${MASTER_FASTA}" >&2
fi

# --- 2. Build BLAST db from novel-panel genomes ------------------------------
echo "[$(date)] Building BLAST database from 15 novel-panel genomes..."
DB_PREFIX="${OUTDIR}/novel_panel.db"
makeblastdb -in "${NOVEL_FASTA}" -dbtype nucl -out "${DB_PREFIX}" -title "novel_panel_15"

# --- 3. BLASTn dark matter against novel panel -------------------------------
echo "[$(date)] Running BLASTn..."
HITS="${OUTDIR}/dark_matter_vs_novel_panel.tsv"
blastn \
    -query "${DARK_FASTA}" \
    -db "${DB_PREFIX}" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
    -evalue 1e-5 \
    -num_threads 4 \
    -perc_identity 70 \
    -max_target_seqs 5 \
    > "${HITS}"

n_hits=$(wc -l < "${HITS}")
echo "  Raw BLAST hits (perc_identity >= 70, e-value <= 1e-5): ${n_hits}"

# --- 4. Filter hits at publication thresholds --------------------------------
# Keep hits at >= 80% identity AND >= 50% query coverage
SUMMARY="${OUTDIR}/hits_summary.tsv"
{
    echo -e "query_contig\tsubject_genome\tpct_identity\taln_length\tqlen\tslen\tqcov_pct\tevalue\tbitscore"
    awk -F'\t' -v OFS='\t' '
        $3 >= 80 && $15 >= 50 {
            print $1, $2, $3, $4, $13, $14, $15, $11, $12
        }
    ' "${HITS}"
} > "${SUMMARY}"

n_strong=$(( $(wc -l < "${SUMMARY}") - 1 ))
echo "  Strong hits (>=80% ID, >=50% qcov): ${n_strong}"

# Also summarise unique dark-matter contigs with at least one strong hit
uniq_hits="${OUTDIR}/hits_unique_dark_contigs.tsv"
{
    echo -e "dark_contig\tn_novel_panel_hits\ttop_subject\ttop_pct_identity"
    awk -F'\t' -v OFS='\t' 'NR > 1 {print $1, $2, $3}' "${SUMMARY}" \
      | sort -k1,1 -k3,3nr \
      | awk -F'\t' -v OFS='\t' '
            {
                if (prev != $1) {
                    n_hits[$1] = 0
                    top_subj[$1] = $2
                    top_pid[$1] = $3
                }
                n_hits[$1]++
                prev = $1
            }
            END {
                for (c in n_hits) print c, n_hits[c], top_subj[c], top_pid[c]
            }
      '
} > "${uniq_hits}"

n_uniq=$(( $(wc -l < "${uniq_hits}") - 1 ))
echo "  Unique dark-matter contigs with >=1 strong novel-panel hit: ${n_uniq}/${n_extracted}"

# --- 5. Report ---------------------------------------------------------------
cat > "${OUTDIR}/README.md" <<REPORT
# Track C dark-matter vs Track A novel-panel BLAST

Generated: $(date)

## Inputs
- Dark-matter contigs: \`${DARK_IDS}\`  (n=${n_extracted})
- Novel-panel genomes: \`${NOVEL_FASTA}\`  (n=15: 13 has_relatives + 2 dark phages)
- ID mapping (novel-panel patient labels): \`${NOVEL_ID_MAP}\`

## Command
BLASTn, perc_identity >= 70, e-value <= 1e-5, max 5 targets per query.

## Filtering
Strong hits retained at:
- percent identity >= 80
- query coverage >= 50%

## Results
- Raw hits: ${n_hits}
- Strong hits: ${n_strong}
- Unique dark-matter contigs with >=1 strong hit: ${n_uniq} / ${n_extracted}

See \`hits_summary.tsv\` for the hit table and \`hits_unique_dark_contigs.tsv\`
for a one-row-per-contig summary.

## Interpretation
- Any hit to one of the 2 *dark phages* (which had zero MetaVR v5 match) is a
  gold-standard positive control: it shows Track C's R > 10 filter recovers
  exactly the kind of novel material the database evidence lines miss.
- A hit to a *has_relatives* novel-panel genome means the dark-matter contig
  is related to recently-discovered vaginal phages that are not yet in MetaVR v5.
- Zero hits is also informative: it means Track C dark matter is capturing
  different novelty than the Track A novel panel (a distinct sequence space
  rather than reassembly of the same genomes).
REPORT

echo "[$(date)] Done. Outputs in: ${OUTDIR}/"
echo "  Summary: ${SUMMARY}"
echo "  Per-contig: ${uniq_hits}"
echo "  Report: ${OUTDIR}/README.md"
