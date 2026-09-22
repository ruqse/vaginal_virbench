#!/usr/bin/env python3
"""
build_ground_truth.py — Tier Assignment from Multi-Evidence Consensus

Aggregates evidence from all ground truth scripts (E1a-E5) and assigns
each contig a confidence tier and category per the study strategy (§4.4).

Tier logic:
  count = sum(E1, E2, E3, E4, E5)   # E1 = E1a OR E1b
  Tier 1: count >= 3  → high confidence viral positive
  Tier 2: count == 2  → medium confidence viral positive
  Tier 3: count == 1  → excluded from primary benchmarking
  Tier 0: count == 0 AND bacterial taxonomy → true negative
  AMBIGUOUS: count == 0, no taxonomy assignment

Category assignment:
  E3 triggered → prophage (with coordinates from Phigaro)
  E5 triggered + eukaryotic markers → eukaryotic_virus
  E1 or E4 without E3 → free_phage
  No evidence + bacterial taxonomy → bacterial

CLI flags:
  --merge-e1-e2: Treat E1+E2 as single evidence line (correlation check)
  --exclude-evidence: Comma-separated evidence lines to exclude (e.g., "E4" or "E1b,E4")
  --evidence-dir: Path to evidence TSV directory
  --output: Path to output ground truth TSV
  --contigs: Path to input contigs FASTA (for length info)

Usage:
    SAMPLE=UC028_V2

    python build_ground_truth.py \
        --evidence-dir results/test_real/ground_truth/$SAMPLE/ \
        --contigs results/test_real/spades/${SAMPLE}_contigs.fasta \
        --output results/test_real/ground_truth/$SAMPLE/ground_truth.tsv

    # Sensitivity analysis: merge E1+E2
    python build_ground_truth.py \
        --evidence-dir results/test_real/ground_truth/$SAMPLE/ \
        --contigs results/test_real/spades/${SAMPLE}_contigs.fasta \
        --output results/test_real/ground_truth/$SAMPLE/ground_truth_merged_e1e2.tsv \
        --merge-e1-e2

    # Sensitivity analysis: exclude E4 (CRISPR)
    python build_ground_truth.py \
        --evidence-dir results/test_real/ground_truth/$SAMPLE/ \
        --contigs results/test_real/spades/${SAMPLE}_contigs.fasta \
        --output results/test_real/ground_truth/$SAMPLE/ground_truth_no_e4.tsv \
        --exclude-evidence E4

    # Sensitivity analysis: exclude E1b (IMG/VR circularity check)
    python build_ground_truth.py \
        --evidence-dir results/test_real/ground_truth/$SAMPLE/ \
        --contigs results/test_real/spades/${SAMPLE}_contigs.fasta \
        --output results/test_real/ground_truth/$SAMPLE/ground_truth_no_e1b.tsv \
        --exclude-evidence E1b
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
import re


# =============================================================================
# Eukaryotic viral family classification (§7.5.2)
# =============================================================================

EUKARYOTIC_FAMILY_KEYWORDS = {
    "papillomavirus": "Papillomaviridae",
    "papilloma virus": "Papillomaviridae",
    "herpesvirus": "Herpesviridae",
    "herpes virus": "Herpesviridae",
    "torque teno": "Anelloviridae",
    "anellovirus": "Anelloviridae",
    "alphatorquevirus": "Anelloviridae",
    "betatorquevirus": "Anelloviridae",
    "gammatorquevirus": "Anelloviridae",
    "genomovirus": "Genomoviridae",
    "gemycircularvirus": "Genomoviridae",
    "polyomavirus": "Polyomaviridae",
    "polyoma virus": "Polyomaviridae",
}


def classify_viral_family(stitle: str) -> str:
    """Extract eukaryotic viral family from DIAMOND BLASTx stitle.

    Parses bracketed organism name (NCBI convention, e.g.,
    '[Human papillomavirus type 16]') then falls back to full stitle.
    Returns family name (e.g., 'Papillomaviridae'), 'phage', or 'unknown'.
    """
    if not stitle:
        return 'unknown'

    # Try to extract organism from bracketed portion [Organism name]
    bracket_match = re.search(r'\[([^\]]+)\]', stitle)
    search_text = bracket_match.group(1) if bracket_match else stitle
    search_lower = search_text.lower()

    # Check eukaryotic family keywords (case-insensitive)
    for keyword, family in EUKARYOTIC_FAMILY_KEYWORDS.items():
        if keyword in search_lower:
            return family

    # Fallback: check full stitle if bracket didn't match
    if bracket_match:
        stitle_lower = stitle.lower()
        for keyword, family in EUKARYOTIC_FAMILY_KEYWORDS.items():
            if keyword in stitle_lower:
                return family

    # Not eukaryotic — classify as phage or unknown
    full_lower = stitle.lower()
    if 'phage' in full_lower or 'virus' in full_lower:
        return 'phage'

    return 'unknown'


# =============================================================================
# Evidence Loaders
# =============================================================================

def load_evidence_1a(filepath):
    """Load DIAMOND BLASTx evidence (E1a)."""
    evidence = {}
    if not os.path.isfile(filepath):
        print(f"  [E1a] File not found: {filepath} — skipping")
        return evidence

    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            cid = row.get('contig_id', '')
            if cid and row.get('viral_hit', '') == 'TRUE':
                stitle = row.get('stitle', '')
                evidence[cid] = {
                    'best_hit': row.get('best_hit', ''),
                    'evalue': row.get('evalue', ''),
                    'pident': row.get('pident', ''),
                    'stitle': stitle,
                    'viral_family': classify_viral_family(stitle),
                }
    print(f"  [E1a] DIAMOND BLASTx: {len(evidence)} viral hits")
    return evidence


def load_evidence_1b(filepath):
    """Load BLASTn IMG/VR evidence (E1b)."""
    evidence = {}
    if not os.path.isfile(filepath):
        print(f"  [E1b] File not found: {filepath} — skipping (data deferred)")
        return evidence

    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            cid = row.get('contig_id', '')
            if cid and row.get('viral_hit', '') == 'TRUE':
                evidence[cid] = {
                    'imgvr_hit': row.get('imgvr_hit', ''),
                    'ani': row.get('ani', ''),
                }
    print(f"  [E1b] BLASTn IMG/VR: {len(evidence)} viral hits")
    return evidence


def load_evidence_2(filepath):
    """Load CheckV structural evidence (E2)."""
    evidence = {}
    if not os.path.isfile(filepath):
        print(f"  [E2] File not found: {filepath} — skipping")
        return evidence

    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            cid = row.get('contig_id', '')
            if cid and row.get('evidence_triggered', '') == 'TRUE':
                evidence[cid] = {
                    'checkv_quality': row.get('checkv_quality', ''),
                    'completeness': row.get('completeness', ''),
                    'provirus': row.get('provirus', ''),
                    'viral_genes': row.get('viral_genes', '0'),
                    'host_genes': row.get('host_genes', '0'),
                    'trigger_reasons': row.get('trigger_reasons', ''),
                }
    print(f"  [E2] CheckV structural: {len(evidence)} triggered")
    return evidence


def load_evidence_3(filepath):
    """Load Phigaro prophage evidence (E3)."""
    evidence = {}
    if not os.path.isfile(filepath):
        print(f"  [E3] File not found: {filepath} — skipping")
        return evidence

    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            cid = row.get('contig_id', '')
            if cid and row.get('evidence_triggered', '') == 'TRUE':
                coords = f"{row.get('prophage_start', '?')}-{row.get('prophage_end', '?')}"
                if cid in evidence:
                    # Multiple prophage regions on same contig
                    evidence[cid]['coords'].append(coords)
                else:
                    evidence[cid] = {
                        'coords': [coords],
                        'n_genes': row.get('n_genes', 'NA'),
                    }
    print(f"  [E3] Phigaro prophage: {len(evidence)} contigs with prophages")
    return evidence


def load_evidence_4(filepath):
    """Load CRISPR spacer matching evidence (E4)."""
    evidence = {}
    if not os.path.isfile(filepath):
        print(f"  [E4] File not found: {filepath} — skipping")
        return evidence

    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            cid = row.get('target_contig_id', '')
            if cid and row.get('evidence_triggered', '') == 'TRUE':
                evidence[cid] = {
                    'spacer_id': row.get('spacer_id', ''),
                    'source_contig': row.get('source_array_contig', ''),
                    'pident': row.get('pident', ''),
                }
    print(f"  [E4] CRISPR spacers: {len(evidence)} target contigs")
    return evidence


def load_evidence_5(filepath, mode='assembly'):
    """Load RCA cross-validation evidence (E5).

    Supports two modes:
      - 'assembly': Original assembly-to-assembly crossmap (run_rca_crossmap.sh)
                     Columns: shotgun_contig_id, rca_contig_id, ani, aligned_frac, evidence_triggered
      - 'read':     Read-level RCA mapping (run_rca_read_mapping.sh)
                     Columns: shotgun_contig_id, contig_length, n_read_pairs,
                              breadth_pct, mean_depth, rca_status, evidence_triggered
    """
    evidence = {}
    if not os.path.isfile(filepath):
        print(f"  [E5] File not found: {filepath} — skipping")
        return evidence

    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            cid = row.get('shotgun_contig_id', '')
            if cid and row.get('evidence_triggered', '') == 'TRUE':
                if mode == 'read':
                    evidence[cid] = {
                        'n_read_pairs': row.get('n_read_pairs', ''),
                        'breadth_pct': row.get('breadth_pct', ''),
                        'mean_depth': row.get('mean_depth', ''),
                        'rca_status': row.get('rca_status', ''),
                    }
                else:
                    evidence[cid] = {
                        'rca_contig': row.get('rca_contig_id', ''),
                        'ani': row.get('ani', ''),
                        'aligned_frac': row.get('aligned_frac', ''),
                    }

    mode_label = "read-level" if mode == 'read' else "assembly-based"
    print(f"  [E5] RCA cross-validation ({mode_label}): {len(evidence)} concordant contigs")
    return evidence


def load_checkv_host_genes(evidence_2_path):
    """
    Load host_genes counts from CheckV E2 evidence for bacterial assignment.
    Used as a lightweight proxy for taxonomy alongside Kraken2.
    """
    host_genes = {}
    if not os.path.isfile(evidence_2_path):
        return host_genes

    with open(evidence_2_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            cid = row.get('contig_id', '')
            try:
                hg = int(row.get('host_genes', '0'))
                vg = int(row.get('viral_genes', '0'))
            except (ValueError, TypeError):
                hg = 0
                vg = 0
            if cid:
                host_genes[cid] = {'host_genes': hg, 'viral_genes': vg}
    return host_genes


def load_kraken2_taxonomy(filepath):
    """
    Load Kraken2 taxonomic assignments for bacterial/archaeal contig identification.

    Kraken2 output format (tab-separated, with --use-names):
      C/U  contig_id  taxid (taxon name)  length  k-mer hits

    Returns dict of contig_id → taxon_name for classified contigs,
    excluding those classified as Viruses (which should be captured by
    evidence lines E1/E2 instead).
    """
    taxonomy = {}
    n_viral_excluded = 0
    if not filepath or not os.path.isfile(filepath):
        return taxonomy

    with open(filepath) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            status = fields[0]
            contig_id = fields[1]
            taxon = fields[2]  # e.g. "Lactobacillus crispatus (taxid 47770)"

            if status != 'C':  # Only classified contigs
                continue

            # Exclude contigs classified as Viruses — these should be
            # captured by the viral evidence lines (E1/E2), not treated
            # as true negatives. The Kraken2 --use-names format puts the
            # lowest-level taxon name here, so we check for known viral
            # indicators. This is conservative: unrecognized taxa are
            # included (most will be bacterial in vaginal metagenomes).
            taxon_lower = taxon.lower()
            if any(term in taxon_lower for term in ('virus', 'viridae', 'phage')):
                n_viral_excluded += 1
                continue

            taxonomy[contig_id] = taxon

    if n_viral_excluded > 0:
        print(f"  [Kraken2] Excluded {n_viral_excluded} virus-classified contigs from Tier 0")

    return taxonomy


def get_contig_lengths(contigs_fasta):
    """Get contig lengths from FASTA file."""
    lengths = {}
    try:
        from Bio import SeqIO
        for rec in SeqIO.parse(contigs_fasta, 'fasta'):
            lengths[rec.id] = len(rec.seq)
    except ImportError:
        # Fallback: manual FASTA parsing
        current_id = None
        current_len = 0
        with open(contigs_fasta) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        lengths[current_id] = current_len
                    current_id = line[1:].split()[0]
                    current_len = 0
                else:
                    current_len += len(line)
            if current_id:
                lengths[current_id] = current_len
    return lengths


# =============================================================================
# Tier Assignment Logic
# =============================================================================

def assign_tier_and_category(contig_id, e1a, e1b, e2, e3, e4, e5,
                              host_gene_info, merge_e1_e2=False,
                              kraken2_taxonomy=None):
    """
    Assign confidence tier and category to a contig.

    Returns: (tier, category, evidence_lines, n_evidence, notes)
    """
    # Determine which evidence lines are triggered
    has_e1a = contig_id in e1a
    has_e1b = contig_id in e1b
    has_e1 = has_e1a or has_e1b
    has_e2 = contig_id in e2
    has_e3 = contig_id in e3
    has_e4 = contig_id in e4
    has_e5 = contig_id in e5

    # Build evidence list
    evidence_lines = []
    if has_e1a:
        evidence_lines.append('E1a')
    if has_e1b:
        evidence_lines.append('E1b')
    if has_e2:
        evidence_lines.append('E2')
    if has_e3:
        evidence_lines.append('E3')
    if has_e4:
        evidence_lines.append('E4')
    if has_e5:
        evidence_lines.append('E5')

    # Count independent evidence lines
    if merge_e1_e2:
        # Treat E1+E2 as a single evidence line (correlation sensitivity check)
        has_e1_or_e2 = has_e1 or has_e2
        n_evidence = sum([has_e1_or_e2, has_e3, has_e4, has_e5])
    else:
        # Standard: each evidence line is independent
        n_evidence = sum([has_e1, has_e2, has_e3, has_e4, has_e5])

    # Notes
    notes = []

    # --- Tier assignment ---
    if n_evidence >= 3:
        tier = 1
    elif n_evidence == 2:
        tier = 2
    elif n_evidence == 1:
        tier = 3  # Excluded from primary benchmarking
    else:
        # n_evidence == 0: check for bacterial taxonomy
        hg_info = host_gene_info.get(contig_id, {})
        host_g = hg_info.get('host_genes', 0)
        viral_g = hg_info.get('viral_genes', 0)

        if host_g > 0 and viral_g == 0:
            tier = 0  # True negative (bacterial via CheckV)
            notes.append(f'checkv:{viral_g}_viral_genes')
        elif kraken2_taxonomy and contig_id in kraken2_taxonomy:
            tier = 0  # True negative (bacterial via Kraken2)
            taxon = kraken2_taxonomy[contig_id]
            notes.append(f'kraken2:{taxon[:60]}')
        else:
            tier = -1  # Ambiguous
            if host_g == 0 and viral_g == 0:
                notes.append('no_gene_classification')

    # --- Category assignment ---
    if tier <= 0 and n_evidence == 0:
        if tier == 0:
            category = 'bacterial'
        else:
            category = 'ambiguous'
    elif has_e3:
        category = 'prophage'
        # Add Phigaro coordinates to notes
        if contig_id in e3:
            coords = e3[contig_id].get('coords', [])
            if coords:
                notes.append(f"phigaro:{';'.join(coords)}")
    else:
        # Check for eukaryotic virus via E1a family assignment (§7.5.2)
        e1a_info = e1a.get(contig_id, {})
        e1a_family = e1a_info.get('viral_family', 'unknown')
        is_eukaryotic = e1a_family not in ('phage', 'unknown')

        if is_eukaryotic:
            category = 'eukaryotic_virus'
            notes.append(f'family:{e1a_family}')
            if has_e5:
                notes.append('rca_validated')
        elif has_e5:
            category = 'free_phage'
            notes.append('rca_validated')
        elif has_e4:
            category = 'free_phage'
            spacer_info = e4.get(contig_id, {})
            notes.append(f"crispr_spacer_hit:{spacer_info.get('spacer_id', '')}")
        elif has_e1 or has_e2:
            category = 'free_phage'
        else:
            category = 'ambiguous'

    return (
        tier,
        category,
        ','.join(evidence_lines) if evidence_lines else '-',
        n_evidence,
        '; '.join(notes) if notes else ''
    )


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Build multi-evidence consensus ground truth with tiered confidence'
    )
    parser.add_argument(
        '--evidence-dir',
        required=True,
        help='Directory containing evidence TSV files'
    )
    parser.add_argument(
        '--contigs',
        required=True,
        help='Path to input contigs FASTA (for contig lengths)'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output ground truth TSV path'
    )
    parser.add_argument(
        '--merge-e1-e2',
        action='store_true',
        default=False,
        help='Treat E1+E2 as single evidence line (correlation sensitivity check)'
    )
    parser.add_argument(
        '--exclude-evidence',
        default='',
        help='Comma-separated evidence lines to exclude (e.g., "E4" or "E1b,E4"). '
             'Used for sensitivity analyses.'
    )
    parser.add_argument(
        '--kraken2-taxonomy',
        default=None,
        help='Path to Kraken2 output TSV (strengthens Tier 0 bacterial assignment)'
    )
    parser.add_argument(
        '--e5-mode',
        choices=['assembly', 'read'],
        default='assembly',
        help='E5 evidence source: "assembly" (run_rca_crossmap.sh, contig-to-contig) '
             'or "read" (run_rca_read_mapping.sh, read-to-contig). Read mode uses '
             'evidence_5_rca_readlevel.tsv; assembly mode uses evidence_5_rca.tsv. '
             '(default: assembly)'
    )
    args = parser.parse_args()

    edir = args.evidence_dir

    print("=" * 70)
    print("Ground Truth Construction — Multi-Evidence Consensus")
    print("=" * 70)
    if args.merge_e1_e2:
        print("  MODE: E1+E2 MERGED (sensitivity analysis)")
    else:
        print("  MODE: Standard (all evidence lines independent)")
    print(f"  E5 MODE: {args.e5_mode} ({'evidence_5_rca_readlevel.tsv' if args.e5_mode == 'read' else 'evidence_5_rca.tsv'})")
    if args.exclude_evidence:
        print(f"  EXCLUDED: {args.exclude_evidence}")
    print()

    # --- Load contig lengths ---
    print("Loading contig lengths...")
    contig_lengths = get_contig_lengths(args.contigs)
    print(f"  Contigs: {len(contig_lengths)}")
    print()

    # --- Load all evidence ---
    print("Loading evidence files...")
    e1a = load_evidence_1a(os.path.join(edir, 'evidence_1a_diamond.tsv'))
    e1b = load_evidence_1b(os.path.join(edir, 'evidence_1b_imgvr.tsv'))
    e2  = load_evidence_2(os.path.join(edir, 'evidence_2_structural.tsv'))
    e3  = load_evidence_3(os.path.join(edir, 'evidence_3_prophage.tsv'))
    e4  = load_evidence_4(os.path.join(edir, 'evidence_4_crispr.tsv'))
    if args.e5_mode == 'read':
        e5 = load_evidence_5(
            os.path.join(edir, 'evidence_5_rca_readlevel.tsv'), mode='read'
        )
    else:
        e5 = load_evidence_5(
            os.path.join(edir, 'evidence_5_rca.tsv'), mode='assembly'
        )
    print()

    # --- Apply evidence exclusion for sensitivity analyses ---
    if args.exclude_evidence:
        excluded = set(x.strip().upper() for x in args.exclude_evidence.split(','))
        print(f"Excluding evidence lines: {', '.join(sorted(excluded))}")
        if 'E1A' in excluded:
            print(f"  E1a zeroed: was {len(e1a)} contigs")
            e1a = {}
        if 'E1B' in excluded:
            print(f"  E1b zeroed: was {len(e1b)} contigs")
            e1b = {}
        if 'E2' in excluded:
            print(f"  E2 zeroed: was {len(e2)} contigs")
            e2 = {}
        if 'E3' in excluded:
            print(f"  E3 zeroed: was {len(e3)} contigs")
            e3 = {}
        if 'E4' in excluded:
            print(f"  E4 zeroed: was {len(e4)} contigs")
            e4 = {}
        if 'E5' in excluded:
            print(f"  E5 zeroed: was {len(e5)} contigs")
            e5 = {}
        print()

    # --- Load host gene info for bacterial assignment ---
    host_gene_info = load_checkv_host_genes(
        os.path.join(edir, 'evidence_2_structural.tsv')
    )

    # --- Load Kraken2 taxonomy (optional, strengthens Tier 0) ---
    kraken2_taxonomy = {}
    if args.kraken2_taxonomy:
        kraken2_taxonomy = load_kraken2_taxonomy(args.kraken2_taxonomy)
        print(f"  [Kraken2] Classified contigs: {len(kraken2_taxonomy)}")
        print()

    # --- Assign tiers ---
    print("Assigning tiers and categories...")
    results = []
    for contig_id in sorted(contig_lengths.keys()):
        tier, category, evidence_lines, n_evidence, notes = assign_tier_and_category(
            contig_id, e1a, e1b, e2, e3, e4, e5,
            host_gene_info, merge_e1_e2=args.merge_e1_e2,
            kraken2_taxonomy=kraken2_taxonomy,
        )
        # Look up viral family from E1a (for downstream §7.5.2 analysis)
        e1a_family = e1a.get(contig_id, {}).get('viral_family', '')
        results.append({
            'contig_id': contig_id,
            'tier': tier,
            'category': category,
            'evidence_lines': evidence_lines,
            'n_evidence': n_evidence,
            'length': contig_lengths.get(contig_id, 0),
            'viral_family': e1a_family,
            'notes': notes,
        })

    # --- Write output ---
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    with open(args.output, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            'contig_id', 'tier', 'category', 'evidence_lines',
            'n_evidence', 'length', 'viral_family', 'notes'
        ])
        for r in results:
            writer.writerow([
                r['contig_id'], r['tier'], r['category'],
                r['evidence_lines'], r['n_evidence'],
                r['length'], r['viral_family'], r['notes']
            ])

    # --- Summary statistics ---
    print()
    print("=" * 70)
    print("GROUND TRUTH SUMMARY")
    print("=" * 70)

    # Tier distribution
    tier_counts = defaultdict(int)
    tier_bp = defaultdict(int)
    for r in results:
        tier_counts[r['tier']] += 1
        tier_bp[r['tier']] += r['length']

    tier_labels = {
        1: 'Tier 1 (high confidence viral)',
        2: 'Tier 2 (medium confidence viral)',
        3: 'Tier 3 (single evidence, excluded)',
        0: 'Tier 0 (true negative, bacterial)',
        -1: 'Ambiguous (excluded)',
    }

    print("\nTier Distribution:")
    print(f"  {'Tier':<45} {'Contigs':>8} {'Total bp':>12}")
    print(f"  {'-'*45} {'-'*8} {'-'*12}")
    for tier in sorted(tier_counts.keys()):
        label = tier_labels.get(tier, f'Tier {tier}')
        print(f"  {label:<45} {tier_counts[tier]:>8,} {tier_bp[tier]:>12,}")
    print(f"  {'TOTAL':<45} {sum(tier_counts.values()):>8,} {sum(tier_bp.values()):>12,}")

    # Category distribution
    cat_counts = defaultdict(int)
    for r in results:
        cat_counts[r['category']] += 1

    print("\nCategory Distribution:")
    for cat in sorted(cat_counts.keys()):
        print(f"  {cat:<20} {cat_counts[cat]:>8,}")

    # Evidence line contribution
    print("\nEvidence Line Contribution (contigs with evidence):")
    e_counts = {
        'E1a (DIAMOND)': len(e1a),
        'E1b (IMG/VR)': len(e1b),
        'E2 (CheckV structural)': len(e2),
        'E3 (Phigaro prophage)': len(e3),
        'E4 (CRISPR spacers)': len(e4),
        'E5 (RCA cross-val)': len(e5),
    }
    for label, count in e_counts.items():
        print(f"  {label:<30} {count:>6}")

    # Benchmarking-ready subset
    tier12 = sum(1 for r in results if r['tier'] in (1, 2))
    tier0 = tier_counts.get(0, 0)
    print(f"\nBenchmarking-ready contigs:")
    print(f"  Viral positives (Tier 1+2): {tier12}")
    print(f"  True negatives (Tier 0):    {tier0}")
    print(f"  Total usable:               {tier12 + tier0}")
    print(f"  Excluded (Tier 3/Ambiguous): {sum(tier_counts.values()) - tier12 - tier0}")

    print(f"\nOutput: {args.output}")


if __name__ == '__main__':
    main()
