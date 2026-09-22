#!/usr/bin/env python3
"""
fragment_genomes.py — Fixed-Length Genome Fragmentation for Spike-In Benchmark

Produces non-overlapping adjacent fragments at specified lengths from spike-in
genomes, consistent with Schackart et al. 2023. Fragments are the direct input
to virus identification tools — bypasses assembly entirely, enabling
length-stratified performance analysis without assembly confounding.

Output:
  - Per-length FASTA: fragments/L{length}_fragments.fasta
  - Manifest TSV:     fragments/fragment_manifest.tsv

Usage:
    python fragment_genomes.py \
        --input data/spike_in/all_viral_genomes.fasta \
        --output data/spike_in/fragments/ \
        --lengths 500 1000 1500 3000 5000 10000

    # With negative controls:
    python fragment_genomes.py \
        --input data/spike_in/all_viral_genomes.fasta \
        --include-negatives data/spike_in/all_negative_controls.fasta \
        --output data/spike_in/fragments/ \
        --lengths 500 1000 1500 3000 5000 10000
"""

import argparse
import csv
import os
import sys
from collections import OrderedDict


# =============================================================================
# FASTA I/O
# =============================================================================

def parse_fasta(filepath):
    """Parse FASTA file into OrderedDict of {id: (description, sequence)}."""
    sequences = OrderedDict()
    current_id = None
    current_desc = None
    current_seq = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = (current_desc, ''.join(current_seq))
                header = line[1:]
                current_id = header.split()[0]
                current_desc = header  # full header line (id + description)
                current_seq = []
            elif current_id is not None:
                current_seq.append(line)
        if current_id is not None:
            sequences[current_id] = (current_desc, ''.join(current_seq))

    return sequences


def infer_category(seq_id):
    """Infer genome category from sequence ID or description."""
    sid = seq_id.lower()
    if 'papillomavirus' in sid or 'hpv' in sid:
        return 'hpv'
    if 'herpes' in sid or 'hsv' in sid:
        return 'herpesvirus'
    if 'torque' in sid or 'teno' in sid or 'anello' in sid or 'ttv' in sid:
        return 'anellovirus'
    if 'gardnerella' in sid and 'phage' in sid:
        return 'gardnerella_phage'
    if 'lactobacillus' in sid and 'phage' in sid:
        return 'lactobacillus_phage'
    if 'megasphaera' in sid and 'phage' in sid:
        return 'megasphaera_phage'
    if 'fannyhessea' in sid and 'phage' in sid:
        return 'fannyhessea_phage'
    if 'sneathia' in sid and 'phage' in sid:
        return 'sneathia_phage'
    # Bacterial negative controls
    if 'lactobacillus' in sid:
        return 'negative_control'
    if 'gardnerella' in sid:
        return 'negative_control'
    if 'megasphaera' in sid:
        return 'negative_control'
    if 'fannyhessea' in sid or 'atopobium' in sid:
        return 'negative_control'
    if 'sneathia' in sid:
        return 'negative_control'
    return 'unknown'


# =============================================================================
# Fragmentation
# =============================================================================

def fragment_sequence(seq_id, sequence, length):
    """
    Split sequence into non-overlapping adjacent fragments of exactly `length` bp.
    Terminal fragment shorter than `length` is discarded (per Schackart et al. 2023).

    Returns list of (fragment_id, fragment_sequence).
    """
    fragments = []
    seq_len = len(sequence)
    n_frags = seq_len // length

    for i in range(n_frags):
        start = i * length
        end = start + length
        frag_seq = sequence[start:end]
        frag_id = f"{seq_id}_frag_{i+1}_L{length}"
        fragments.append((frag_id, frag_seq))

    return fragments


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Fragment spike-in genomes at fixed lengths (Schackart et al. 2023)'
    )
    parser.add_argument(
        '--input', required=True,
        help='Input FASTA of viral spike-in genomes'
    )
    parser.add_argument(
        '--include-negatives', default=None,
        help='Optional FASTA of negative control (bacterial) genomes'
    )
    parser.add_argument(
        '--output', '-o', required=True,
        help='Output directory for fragment files'
    )
    parser.add_argument(
        '--lengths', nargs='+', type=int,
        default=[500, 1000, 1500, 3000, 5000, 10000],
        help='Fragment lengths in bp (default: 500 1000 1500 3000 5000 10000)'
    )
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    # --- Load genomes ---
    print("Loading viral genomes...")
    viral_seqs = parse_fasta(args.input)  # {id: (description, sequence)}
    print(f"  Viral genomes: {len(viral_seqs)}")
    for sid, (desc, seq) in viral_seqs.items():
        print(f"    {sid:<18} {len(seq):>10,} bp  {infer_category(desc):<25} {desc[:60]}")

    neg_seqs = OrderedDict()
    if args.include_negatives and os.path.isfile(args.include_negatives):
        print("Loading negative control genomes...")
        neg_seqs = parse_fasta(args.include_negatives)
        print(f"  Negative controls: {len(neg_seqs)}")
        for sid, (desc, seq) in neg_seqs.items():
            print(f"    {sid:<18} {len(seq):>10,} bp  {infer_category(desc):<25} {desc[:60]}")

    # --- Fragment at each length ---
    manifest_rows = []

    for length in sorted(args.lengths):
        fasta_path = os.path.join(args.output, f"L{length}_fragments.fasta")
        n_viral = 0
        n_neg = 0
        n_skipped = 0

        with open(fasta_path, 'w') as f:
            # Viral fragments
            for seq_id, (desc, sequence) in viral_seqs.items():
                if len(sequence) < length:
                    n_skipped += 1
                    continue
                category = infer_category(desc)
                frags = fragment_sequence(seq_id, sequence, length)
                for frag_id, frag_seq in frags:
                    f.write(f">{frag_id} source={seq_id} category={category} label=viral\n")
                    # Write sequence in 80-char lines
                    for i in range(0, len(frag_seq), 80):
                        f.write(frag_seq[i:i+80] + '\n')
                    manifest_rows.append({
                        'fragment_id': frag_id,
                        'source_genome': seq_id,
                        'category': category,
                        'label': 'viral',
                        'fragment_length': length,
                        'source_genome_length': len(sequence),
                    })
                    n_viral += 1

            # Negative control fragments
            for seq_id, (desc, sequence) in neg_seqs.items():
                if len(sequence) < length:
                    n_skipped += 1
                    continue
                category = infer_category(desc)
                frags = fragment_sequence(seq_id, sequence, length)
                for frag_id, frag_seq in frags:
                    f.write(f">{frag_id} source={seq_id} category={category} label=bacterial\n")
                    for i in range(0, len(frag_seq), 80):
                        f.write(frag_seq[i:i+80] + '\n')
                    manifest_rows.append({
                        'fragment_id': frag_id,
                        'source_genome': seq_id,
                        'category': category,
                        'label': 'bacterial',
                        'fragment_length': length,
                        'source_genome_length': len(sequence),
                    })
                    n_neg += 1

        print(f"  L{length:>5}: {n_viral:>6} viral + {n_neg:>6} bacterial = {n_viral+n_neg:>6} fragments"
              f"  ({n_skipped} genomes too short, skipped)")

    # --- Write manifest ---
    manifest_path = os.path.join(args.output, 'fragment_manifest.tsv')
    with open(manifest_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=[
            'fragment_id', 'source_genome', 'category', 'label',
            'fragment_length', 'source_genome_length'
        ])
        writer.writeheader()
        writer.writerows(manifest_rows)

    # --- Summary ---
    print(f"\n{'='*60}")
    print(f"Fragmentation Summary")
    print(f"{'='*60}")
    print(f"  Input viral genomes:     {len(viral_seqs)}")
    print(f"  Input negative controls: {len(neg_seqs)}")
    print(f"  Fragment lengths:        {args.lengths}")
    print(f"  Total fragments:         {len(manifest_rows)}")

    # Per-category breakdown
    from collections import Counter
    label_counts = Counter(r['label'] for r in manifest_rows)
    cat_counts = Counter(r['category'] for r in manifest_rows)
    print(f"\n  By label:")
    for label, count in sorted(label_counts.items()):
        print(f"    {label:<12} {count:>8}")
    print(f"\n  By category:")
    for cat, count in sorted(cat_counts.items()):
        print(f"    {cat:<25} {count:>8}")

    print(f"\n  Output directory: {args.output}")
    print(f"  Manifest:         {manifest_path}")


if __name__ == '__main__':
    main()
