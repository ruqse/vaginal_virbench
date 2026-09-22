#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -J blast_spikein_vs_assembly
#SBATCH -o logs/blast_spikein_%j.out
#SBATCH -e logs/blast_spikein_%j.err

# ============================================================================
# Check spike-in genome homology against a Track B background assembly
#
# Purpose: Verify whether any of the spike-in viral genomes have close
# homologs already present in a background sample's metagenome assembly.
# Close homologs would cause chimeric assemblies in Track B (spike-in +
# real background co-assembly), confounding ground truth labels.
#
# Usage:
#   sbatch check_spikein_homology.sh UC115_V2    # Check CST-I background
#   sbatch check_spikein_homology.sh UC093_V3    # Check CST-IV background
#
# Output (per sample):
#   results/spike_in/homology_check/${SAMPLE}/blast_hits.tsv
#   results/spike_in/homology_check/${SAMPLE}/filtered_hits.tsv
#   results/spike_in/homology_check/${SAMPLE}/HOMOLOGY_REPORT.md
# ============================================================================

set -euo pipefail

# --- Sample ID from argument ---
SAMPLE="${1:?Usage: sbatch $0 SAMPLE_ID (e.g. UC115_V2, UC093_V3)}"

# --- Paths ---
PROJ=${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}
QUERY="${PROJ}/data/spike_in/all_viral_genomes.fasta"
SUBJECT="${PROJ}/results/test_real/spades/${SAMPLE}_contigs.fasta"
OUTDIR="${PROJ}/results/spike_in/homology_check/${SAMPLE}"
TMPDB="${SNIC_TMP:-${TMPDIR:-/tmp}}/blastdb_${SAMPLE}"

# --- Validate inputs ---
if [[ ! -f "$QUERY" ]]; then
    echo "ERROR: Query file not found: $QUERY" >&2
    exit 1
fi
if [[ ! -f "$SUBJECT" ]]; then
    echo "ERROR: Subject file not found: $SUBJECT" >&2
    exit 1
fi

# --- Setup ---
module load BLAST+/2.17.0-gompi-2024a
mkdir -p "$OUTDIR" "$TMPDB"

echo "=== Spike-in vs ${SAMPLE} homology check ==="
echo "Query:   $QUERY"
echo "Subject: $SUBJECT"
echo "Output:  $OUTDIR"
echo "Date:    $(date)"
echo ""

# --- Build BLAST database in scratch ---
echo "Building BLAST database..."
makeblastdb -in "$SUBJECT" -dbtype nucl -out "${TMPDB}/${SAMPLE}" -title "${SAMPLE}_contigs"
echo "Database built."

# --- Run BLASTn ---
echo "Running BLASTn..."
blastn \
    -query "$QUERY" \
    -db "${TMPDB}/${SAMPLE}" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
    -evalue 1e-5 \
    -num_threads "${SLURM_CPUS_PER_TASK:-4}" \
    -max_target_seqs 50 \
    -out "${OUTDIR}/blast_hits.tsv"

N_HITS=$(wc -l < "${OUTDIR}/blast_hits.tsv")
echo "BLASTn complete: ${N_HITS} raw hits."

# --- Filter and generate report ---
echo "Generating report..."

python3 - "$OUTDIR" "$SAMPLE" << 'PYEOF'
import sys
from collections import defaultdict
from pathlib import Path

OUTDIR = Path(sys.argv[1])
SAMPLE_ID = sys.argv[2] if len(sys.argv) > 2 else "unknown"
BLAST_FILE = OUTDIR / "blast_hits.tsv"
FILTERED_FILE = OUTDIR / "filtered_hits.tsv"
REPORT_FILE = OUTDIR / "HOMOLOGY_REPORT.md"

# Column indices for outfmt 6
QSEQID, SSEQID, PIDENT, LENGTH, MISMATCH, GAPOPEN = 0, 1, 2, 3, 4, 5
QSTART, QEND, SSTART, SEND, EVALUE, BITSCORE, QLEN, SLEN = 6, 7, 8, 9, 10, 11, 12, 13

MIN_IDENTITY = 80.0  # percent
MIN_ALIGNMENT = 500  # bp

# Spike-in genome metadata (for the report)
GENOME_INFO = {
    "NC_007924.1": ("Lactobacillus phage KC5a", "lactobacillus_phage", 38239),
    "NC_011801.1": ("Lactobacillus phage Lv-1", "lactobacillus_phage", 38934),
    "MW387018.1":  ("Gardnerella phage vB_Gva_AB1", "gardnerella_phage", 50268),
    "NC_001526.4": ("HPV-16", "hpv", 7906),
    "NC_001357.1": ("HPV-18", "hpv", 7857),
    "NC_001592.1": ("HPV-52", "hpv", 7942),
    "NC_001443.1": ("HPV-58", "hpv", 7824),
    "NC_002076.2": ("Torque teno virus 1", "anellovirus", 3852),
    "NC_001806.2": ("HSV-1 strain 17", "herpesvirus", 152222),
    "NC_001798.2": ("HSV-2 strain HG52", "herpesvirus", 154675),
}

# --- Parse and filter BLAST hits ---
hits_by_genome = defaultdict(list)  # {qseqid: [(pident, length, qstart, qend, sseqid, evalue, qlen), ...]}
filtered_lines = []

with open(BLAST_FILE) as fh:
    for line in fh:
        fields = line.strip().split("\t")
        if len(fields) < 14:
            continue
        qseqid = fields[QSEQID]
        pident = float(fields[PIDENT])
        alen = int(fields[LENGTH])
        qstart = int(fields[QSTART])
        qend = int(fields[QEND])
        sseqid = fields[SSEQID]
        evalue = float(fields[EVALUE])
        qlen = int(fields[QLEN])

        if pident >= MIN_IDENTITY and alen >= MIN_ALIGNMENT:
            hits_by_genome[qseqid].append({
                "pident": pident,
                "alen": alen,
                "qstart": min(qstart, qend),
                "qend": max(qstart, qend),
                "sseqid": sseqid,
                "evalue": evalue,
                "qlen": qlen,
            })
            filtered_lines.append(line)

# Write filtered hits
with open(FILTERED_FILE, "w") as fh:
    fh.write("qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\n")
    for line in filtered_lines:
        fh.write(line)


def merge_intervals(intervals):
    """Merge overlapping intervals, return total non-overlapping coverage."""
    if not intervals:
        return 0, []
    sorted_ivs = sorted(intervals, key=lambda x: x[0])
    merged = [sorted_ivs[0]]
    for start, end in sorted_ivs[1:]:
        if start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    total = sum(e - s + 1 for s, e in merged)
    return total, merged


# --- Compute per-genome summary ---
summaries = []
for accession in GENOME_INFO:
    name, category, genome_len = GENOME_INFO[accession]
    hits = hits_by_genome.get(accession, [])

    if not hits:
        summaries.append({
            "accession": accession,
            "name": name,
            "category": category,
            "genome_len": genome_len,
            "n_hits": 0,
            "best_pident": 0.0,
            "total_aligned_bp": 0,
            "query_coverage_pct": 0.0,
            "risk": "LOW",
            "top_subject": "-",
        })
        continue

    # Merge overlapping alignment intervals on query
    intervals = [(h["qstart"], h["qend"]) for h in hits]
    total_covered, merged = merge_intervals(intervals)
    coverage_pct = 100.0 * total_covered / genome_len

    best_pident = max(h["pident"] for h in hits)
    best_hit = max(hits, key=lambda h: h["alen"])

    # Risk classification
    if coverage_pct >= 50 and best_pident >= 95:
        risk = "HIGH"
    elif coverage_pct >= 10 and best_pident >= 90:
        risk = "MODERATE"
    else:
        risk = "LOW"

    summaries.append({
        "accession": accession,
        "name": name,
        "category": category,
        "genome_len": genome_len,
        "n_hits": len(hits),
        "best_pident": best_pident,
        "total_aligned_bp": total_covered,
        "query_coverage_pct": round(coverage_pct, 1),
        "risk": risk,
        "top_subject": best_hit["sseqid"],
    })


# --- Write report ---
with open(REPORT_FILE, "w") as rpt:
    rpt.write(f"# Spike-In vs {SAMPLE_ID} Assembly Homology Report\n\n")
    rpt.write("**Purpose**: Check whether spike-in viral genomes have close homologs\n")
    rpt.write(f"already present in the {SAMPLE_ID} metagenome assembly, which would cause\n")
    rpt.write("chimeric assemblies in Track B co-assembly benchmarking.\n\n")
    rpt.write(f"**Filter thresholds**: identity >= {MIN_IDENTITY}%, alignment >= {MIN_ALIGNMENT} bp\n\n")

    # Risk summary
    risk_counts = defaultdict(int)
    for s in summaries:
        risk_counts[s["risk"]] += 1
    rpt.write("## Risk Summary\n\n")
    rpt.write(f"- **HIGH** (>= 50% query coverage, >= 95% identity): {risk_counts['HIGH']} genomes\n")
    rpt.write(f"- **MODERATE** (>= 10% coverage, >= 90% identity): {risk_counts['MODERATE']} genomes\n")
    rpt.write(f"- **LOW** (below thresholds or no hits): {risk_counts['LOW']} genomes\n\n")

    # Risk interpretation
    rpt.write("### Risk Interpretation\n\n")
    rpt.write("- **HIGH**: Native virus very similar to spike-in genome. Track A co-assembly\n")
    rpt.write("  will likely produce chimeric contigs mixing spike-in and native reads.\n")
    rpt.write("  Ground truth labels unreliable for this genome. Consider excluding from Track B.\n")
    rpt.write("- **MODERATE**: Partial homology. Some chimeric contigs possible but most labels\n")
    rpt.write("  likely correct. Interpret Track A results with caution for this genome.\n")
    rpt.write(f"- **LOW**: No close homologs in {SAMPLE_ID}. Track B reliable for this genome.\n\n")

    # Per-genome table
    rpt.write("## Per-Genome Results\n\n")
    rpt.write("| Accession | Name | Category | Genome (bp) | Hits | Best %ID | Covered (bp) | Coverage % | Risk | Top Hit |\n")
    rpt.write("|-----------|------|----------|-------------|------|----------|--------------|------------|------|----------|\n")
    for s in summaries:
        rpt.write(f"| {s['accession']} | {s['name']} | {s['category']} | "
                   f"{s['genome_len']:,} | {s['n_hits']} | {s['best_pident']:.1f} | "
                   f"{s['total_aligned_bp']:,} | {s['query_coverage_pct']}% | "
                   f"**{s['risk']}** | {s['top_subject']} |\n")

    # Recommendations
    rpt.write("\n## Recommendations\n\n")
    high_risk = [s for s in summaries if s["risk"] == "HIGH"]
    mod_risk = [s for s in summaries if s["risk"] == "MODERATE"]
    if high_risk:
        rpt.write("### HIGH-risk genomes (consider excluding from Track B)\n\n")
        for s in high_risk:
            rpt.write(f"- **{s['name']}** ({s['accession']}): {s['query_coverage_pct']}% coverage "
                       f"at {s['best_pident']:.1f}% identity. Native homolog in {SAMPLE_ID} "
                       f"({s['top_subject']}) will produce chimeric assemblies.\n")
        rpt.write("\n")
    if mod_risk:
        rpt.write("### MODERATE-risk genomes (interpret with caution)\n\n")
        for s in mod_risk:
            rpt.write(f"- **{s['name']}** ({s['accession']}): {s['query_coverage_pct']}% coverage "
                       f"at {s['best_pident']:.1f}% identity.\n")
        rpt.write("\n")
    if not high_risk and not mod_risk:
        rpt.write(f"All spike-in genomes show LOW risk — Track B co-assembly with {SAMPLE_ID} is reliable.\n\n")

    rpt.write("---\n")
    from datetime import datetime
    rpt.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")

print(f"Report written: {REPORT_FILE}")
print(f"Filtered hits:  {FILTERED_FILE}")
print(f"\nSummary:")
for s in summaries:
    print(f"  {s['risk']:8s}  {s['name']:30s}  {s['query_coverage_pct']:6.1f}% coverage  {s['best_pident']:5.1f}% identity  ({s['n_hits']} hits)")
PYEOF

echo ""
echo "=== Done ==="
