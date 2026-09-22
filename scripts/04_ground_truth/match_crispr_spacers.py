#!/usr/bin/env python3
"""
match_crispr_spacers.py — Evidence 4b: CRISPR Spacer Matching

Matches extracted CRISPR spacers against shotgun contigs using BLASTn-short.
A contig targeted by spacers from a DIFFERENT contig = phage evidence.

Thresholds (strategy §4.2):
  - Identity   >= 95%
  - Coverage   >= 95% of spacer length

Requires: BLAST+ module loaded, spacers FASTA from run_crisprcasfinder.sh

Usage:
    python match_crispr_spacers.py \
        --spacers results/test_real/ground_truth/${SAMPLE}/crispr_spacers.fasta \
        --contigs results/test_real/spades/${SAMPLE}_contigs.fasta \
        --output results/test_real/ground_truth/${SAMPLE}/evidence_4_crispr.tsv

    # Or with BLAST+ already run:
    python match_crispr_spacers.py \
        --spacers crispr_spacers.fasta \
        --contigs contigs.fasta \
        --blast-output blast_spacers_raw.tsv \
        --output evidence_4_crispr.tsv
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict


def run_blastn_short(spacers_fasta, contigs_fasta, output_file, threads=4):
    """Run BLASTn-short to match spacers against contigs."""
    cmd = [
        'blastn',
        '-task', 'blastn-short',
        '-query', spacers_fasta,
        '-subject', contigs_fasta,
        '-out', output_file,
        '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen',
        '-evalue', '1',
        '-dust', 'yes',
        '-num_threads', str(threads),
    ]

    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"WARNING: BLASTn exited with code {result.returncode}", file=sys.stderr)
        if result.stderr:
            print(f"  STDERR: {result.stderr[:500]}", file=sys.stderr)

    return result.returncode


def get_spacer_source_contig(spacer_id):
    """
    Extract the source contig from a spacer ID.

    Spacer IDs are formatted as: {contig_id}_CRISPR_{n}_spacer_{m}
    We need to recover the original contig_id, which may itself contain underscores.
    """
    # The pattern is: {contig}_{CRISPR}_{n}_{spacer}_{m}
    # Find the last occurrence of _CRISPR_ to split
    idx = spacer_id.rfind('_CRISPR_')
    if idx > 0:
        return spacer_id[:idx]
    # Fallback: return everything before last 3 underscore-separated parts
    parts = spacer_id.rsplit('_', 4)
    if len(parts) >= 5:
        return '_'.join(parts[:-4])
    return spacer_id


def parse_spacer_lengths(spacers_fasta):
    """Get spacer lengths from FASTA file."""
    lengths = {}
    current_id = None
    current_seq = []

    with open(spacers_fasta) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id and current_seq:
                    lengths[current_id] = len(''.join(current_seq))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id and current_seq:
            lengths[current_id] = len(''.join(current_seq))

    return lengths


def main():
    parser = argparse.ArgumentParser(
        description='Match CRISPR spacers to contigs for Evidence 4'
    )
    parser.add_argument(
        '--spacers',
        required=True,
        help='Path to spacers FASTA (from run_crisprcasfinder.sh)'
    )
    parser.add_argument(
        '--contigs',
        required=True,
        help='Path to shotgun contigs FASTA'
    )
    parser.add_argument(
        '--blast-output',
        default=None,
        help='Pre-computed BLASTn output (skip running BLAST)'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output evidence TSV path'
    )
    parser.add_argument(
        '--min-identity',
        type=float,
        default=95.0,
        help='Minimum percent identity (default: 95)'
    )
    parser.add_argument(
        '--min-coverage',
        type=float,
        default=95.0,
        help='Minimum spacer coverage percent (default: 95)'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of threads for BLASTn (default: 4)'
    )
    parser.add_argument(
        '--min-spacer-length',
        type=int,
        default=25,
        help='Minimum spacer length in bp (default: 25; real CRISPR spacers are 26-72 bp)'
    )
    parser.add_argument(
        '--max-targets-per-spacer',
        type=int,
        default=50,
        help='Max target contigs per spacer; higher = promiscuous/conserved sequence (default: 50)'
    )
    args = parser.parse_args()

    # Validate inputs
    if not os.path.isfile(args.spacers):
        print(f"ERROR: Spacers FASTA not found: {args.spacers}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(args.contigs):
        print(f"ERROR: Contigs FASTA not found: {args.contigs}", file=sys.stderr)
        sys.exit(1)

    # Check if spacers file is empty
    spacer_lengths_all = parse_spacer_lengths(args.spacers)
    if not spacer_lengths_all:
        print("WARNING: No spacers found in FASTA. Creating empty evidence file.")
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        with open(args.output, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow([
                'target_contig_id', 'spacer_id', 'source_array_contig',
                'pident', 'coverage', 'evidence_triggered'
            ])
        sys.exit(0)

    # Filter spacers by minimum length (removes microsatellite-length fragments)
    spacer_lengths = {k: v for k, v in spacer_lengths_all.items()
                      if v >= args.min_spacer_length}
    spacers_filtered_count = len(spacer_lengths_all) - len(spacer_lengths)

    print(f"Spacers in FASTA:                 {len(spacer_lengths_all)}")
    print(f"Spacers after length filter (>= {args.min_spacer_length} bp): {len(spacer_lengths)}")
    if spacers_filtered_count > 0:
        print(f"  Removed {spacers_filtered_count} short spacers")
    print(f"Thresholds: identity >= {args.min_identity}%, coverage >= {args.min_coverage}%")

    if not spacer_lengths:
        print("WARNING: No spacers pass length filter. Creating empty evidence file.")
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        with open(args.output, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow([
                'target_contig_id', 'spacer_id', 'source_array_contig',
                'pident', 'coverage', 'evidence_triggered'
            ])
        sys.exit(0)

    # Run BLASTn if needed
    if args.blast_output and os.path.isfile(args.blast_output):
        blast_file = args.blast_output
        print(f"Using pre-computed BLAST output: {blast_file}")
    else:
        blast_file = os.path.join(os.path.dirname(args.output), 'blast_spacers_raw.tsv')
        print(f"\nRunning BLASTn-short...")
        ret = run_blastn_short(args.spacers, args.contigs, blast_file, args.threads)
        if ret != 0:
            print("WARNING: BLASTn may have failed. Continuing with available output.")

    # Parse BLAST output
    # Columns: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
    hits = []
    if os.path.isfile(blast_file):
        with open(blast_file) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 14:
                    continue

                spacer_id = fields[0]
                target_contig = fields[1]
                pident = float(fields[2])
                alen = int(fields[3])
                qlen = int(fields[12])

                # Skip spacers that didn't pass length filter
                if spacer_id not in spacer_lengths:
                    continue

                # Compute spacer coverage
                coverage = (alen / qlen) * 100 if qlen > 0 else 0

                # Get source contig of the spacer
                source_contig = get_spacer_source_contig(spacer_id)

                # CRITICAL: Only count hits to DIFFERENT contigs
                # A spacer hitting its own contig is the CRISPR array itself
                if target_contig == source_contig:
                    continue

                # Apply thresholds
                if pident >= args.min_identity and coverage >= args.min_coverage:
                    hits.append({
                        'target_contig_id': target_contig,
                        'spacer_id': spacer_id,
                        'source_array_contig': source_contig,
                        'pident': pident,
                        'coverage': coverage,
                    })

    print(f"\nFiltered hits (different contig, identity >= {args.min_identity}%, "
          f"coverage >= {args.min_coverage}%): {len(hits)}")

    # Filter promiscuous spacers (conserved sequences that hit too many targets)
    spacer_target_counts = Counter(h['spacer_id'] for h in hits)
    promiscuous = {s for s, c in spacer_target_counts.items()
                   if c > args.max_targets_per_spacer}
    hits_before_promiscuity = len(hits)
    if promiscuous:
        hits = [h for h in hits if h['spacer_id'] not in promiscuous]
        print(f"  Removed {len(promiscuous)} promiscuous spacers "
              f"(>{args.max_targets_per_spacer} targets each)")
        print(f"  Hits after promiscuity filter: {len(hits)} "
              f"(was {hits_before_promiscuity})")

    # Group by target contig — keep best hit per target
    best_per_target = {}
    for hit in hits:
        tc = hit['target_contig_id']
        if tc not in best_per_target or hit['pident'] > best_per_target[tc]['pident']:
            best_per_target[tc] = hit

    # Write evidence TSV
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    with open(args.output, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            'target_contig_id', 'spacer_id', 'source_array_contig',
            'pident', 'coverage', 'evidence_triggered'
        ])

        for target_contig in sorted(best_per_target.keys()):
            hit = best_per_target[target_contig]
            writer.writerow([
                hit['target_contig_id'],
                hit['spacer_id'],
                hit['source_array_contig'],
                f"{hit['pident']:.1f}",
                f"{hit['coverage']:.1f}",
                'TRUE'
            ])

    print(f"\n{'='*60}")
    print(f"Evidence 4 — Diagnostic Summary")
    print(f"{'='*60}")
    print(f"  Spacers in FASTA (total):        {len(spacer_lengths_all)}")
    print(f"  Spacers after length filter:     {len(spacer_lengths)}")
    print(f"  Spacers removed (short):         {spacers_filtered_count}")
    print(f"  Spacers removed (promiscuous):   {len(promiscuous)}")
    print(f"  BLAST hits (threshold-passing):  {hits_before_promiscuity}")
    print(f"  Hits after promiscuity filter:   {len(hits)}")
    print(f"  Target contigs with spacer hits: {len(best_per_target)}")
    print(f"  Unique source array contigs:     {len(set(h['source_array_contig'] for h in best_per_target.values()))}")
    print(f"\nOutput: {args.output}")


if __name__ == '__main__':
    main()
