#!/usr/bin/env python3
"""
label_spike_in_contigs.py — Absolute Ground Truth Labeling for Spike-In Assemblies

Maps assembled contigs back to spike-in reference genomes via minimap2 to assign
absolute TP/TN labels. This replaces the multi-evidence consensus as the primary
ground truth for headline benchmarking metrics.

Labels:
  viral     — contig maps to a spike-in viral genome (>=95% id, >=80% AF)
  bacterial — contig maps to a negative control genome (>=95% id, >=80% AF)
  background — contig from real metagenome background (maps to neither)

Usage:
    python label_spike_in_contigs.py \
        --contigs data/spike_in/assemblies/assembly_cov10/contigs_filtered.fasta \
        --viral-genomes data/spike_in/all_viral_genomes.fasta \
        --negative-genomes data/spike_in/all_negative_controls.fasta \
        --output data/spike_in/assemblies/assembly_cov10/ground_truth.tsv

    # Batch mode: label all coverage levels
    for cov in 0.1 0.5 1 5 10 50; do
        python label_spike_in_contigs.py \
            --contigs data/spike_in/assemblies/assembly_cov${cov}/contigs_filtered.fasta \
            --viral-genomes data/spike_in/all_viral_genomes.fasta \
            --negative-genomes data/spike_in/all_negative_controls.fasta \
            --output data/spike_in/assemblies/assembly_cov${cov}/ground_truth.tsv
    done
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile
from collections import defaultdict


def get_contig_lengths(fasta_path):
    """Get contig lengths from FASTA file."""
    lengths = {}
    current_id = None
    current_len = 0
    with open(fasta_path) as f:
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


def run_minimap2(query_fasta, ref_fasta, threads=4):
    """
    Run minimap2 in asm5 mode (>=95% identity assembly mapping).
    Returns list of alignment dicts.
    """
    cmd = [
        'minimap2',
        '-x', 'asm5',
        '-t', str(threads),
        '--secondary=no',
        '-c',  # output CIGAR in PAF
        ref_fasta,
        query_fasta,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"WARNING: minimap2 exited with code {result.returncode}", file=sys.stderr)
        if result.stderr:
            print(f"  STDERR: {result.stderr[:300]}", file=sys.stderr)

    alignments = []
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        fields = line.split('\t')
        if len(fields) < 12:
            continue

        qname = fields[0]
        qlen = int(fields[1])
        qstart = int(fields[2])
        qend = int(fields[3])
        tname = fields[5]
        tlen = int(fields[6])
        nmatch = int(fields[9])
        block_len = int(fields[10])
        mapq = int(fields[11])

        # Compute alignment identity and query coverage
        pident = (nmatch / block_len * 100) if block_len > 0 else 0
        query_cov = ((qend - qstart) / qlen * 100) if qlen > 0 else 0

        alignments.append({
            'query': qname,
            'query_len': qlen,
            'target': tname,
            'target_len': tlen,
            'pident': pident,
            'query_cov': query_cov,
            'mapq': mapq,
        })

    return alignments


def infer_category(genome_id):
    """Infer category from genome accession/description."""
    gid = genome_id.lower()
    if 'papillomavirus' in gid or 'hpv' in gid:
        return 'hpv'
    if 'herpes' in gid or 'hsv' in gid:
        return 'herpesvirus'
    if 'torque' in gid or 'teno' in gid or 'anello' in gid:
        return 'anellovirus'
    if 'gardnerella' in gid and 'phage' in gid:
        return 'gardnerella_phage'
    if 'lactobacillus' in gid and 'phage' in gid:
        return 'lactobacillus_phage'
    if 'megasphaera' in gid and 'phage' in gid:
        return 'megasphaera_phage'
    if 'fannyhessea' in gid and 'phage' in gid:
        return 'fannyhessea_phage'
    if 'sneathia' in gid and 'phage' in gid:
        return 'sneathia_phage'
    # Accession-based fallback
    acc = genome_id.split('.')[0] if '.' in genome_id else genome_id
    # Handle IMGVR accessions (MetaVR v5 phages)
    imgvr_map = {
        'IMGVR_UViG_3300045988_117760': 'megasphaera_phage',
        'IMGVR_UViG_3300045988_083571': 'megasphaera_phage',
        'IMGVR_UViG_3300007272_000043': 'fannyhessea_phage',
        'IMGVR_UViG_2636415511_000001': 'sneathia_phage',
    }
    if genome_id in imgvr_map:
        return imgvr_map[genome_id]
    acc_map = {
        'NC_001526': 'hpv', 'NC_001357': 'hpv', 'NC_001592': 'hpv', 'NC_001443': 'hpv',
        'NC_001806': 'herpesvirus', 'NC_001798': 'herpesvirus',
        'NC_002076': 'anellovirus',
        'NC_011801': 'lactobacillus_phage', 'NC_007924': 'lactobacillus_phage',
        'MW387018': 'gardnerella_phage',
        'CP072197': 'negative_control', 'NC_014644': 'negative_control',
    }
    if acc in acc_map:
        return acc_map[acc]
    if acc.startswith('NZ_'):
        return 'negative_control'
    return 'unknown'


def main():
    parser = argparse.ArgumentParser(
        description='Label spike-in assembled contigs with absolute ground truth'
    )
    parser.add_argument(
        '--contigs', required=True,
        help='Assembled contigs FASTA (from co-assembly)'
    )
    parser.add_argument(
        '--viral-genomes', required=True,
        help='Spike-in viral genomes FASTA'
    )
    parser.add_argument(
        '--negative-genomes', default=None,
        help='Negative control bacterial genomes FASTA'
    )
    parser.add_argument(
        '--min-identity', type=float, default=95.0,
        help='Minimum alignment identity %% (default: 95)'
    )
    parser.add_argument(
        '--min-coverage', type=float, default=80.0,
        help='Minimum query coverage %% (default: 80)'
    )
    parser.add_argument(
        '--threads', type=int, default=4,
        help='Threads for minimap2 (default: 4)'
    )
    parser.add_argument(
        '--exclude-genomes', default=None,
        help='Comma-separated accessions to exclude from viral labeling '
             '(e.g., MW387018.1 for vB_Gva_AB1 on UC093_V3 due to native homolog). '
             'Contigs mapping to excluded genomes are labeled "excluded_chimeric_risk".'
    )
    parser.add_argument(
        '--output', '-o', required=True,
        help='Output ground truth TSV'
    )
    args = parser.parse_args()

    # Parse excluded genomes
    excluded_genomes = set()
    if args.exclude_genomes:
        excluded_genomes = {g.strip() for g in args.exclude_genomes.split(',')}
        print(f"Excluding {len(excluded_genomes)} genomes from viral labeling: {excluded_genomes}")

    # Validate inputs
    if not os.path.isfile(args.contigs):
        print(f"ERROR: Contigs not found: {args.contigs}", file=sys.stderr)
        sys.exit(1)

    # Get contig lengths
    contig_lengths = get_contig_lengths(args.contigs)
    print(f"Contigs: {len(contig_lengths)}")

    # --- Map contigs to viral genomes ----------------------------------------
    print(f"\nMapping contigs to viral genomes...")
    viral_alns = run_minimap2(args.contigs, args.viral_genomes, args.threads)
    print(f"  Raw alignments: {len(viral_alns)}")

    # Filter by identity and coverage, keep best hit per contig
    viral_hits = {}
    for aln in viral_alns:
        if aln['pident'] >= args.min_identity and aln['query_cov'] >= args.min_coverage:
            qname = aln['query']
            if qname not in viral_hits or aln['pident'] > viral_hits[qname]['pident']:
                viral_hits[qname] = aln
    print(f"  Viral hits (>={args.min_identity}% id, >={args.min_coverage}% cov): {len(viral_hits)}")

    # --- Map contigs to negative genomes -------------------------------------
    neg_hits = {}
    if args.negative_genomes and os.path.isfile(args.negative_genomes):
        print(f"\nMapping contigs to negative control genomes...")
        neg_alns = run_minimap2(args.contigs, args.negative_genomes, args.threads)
        print(f"  Raw alignments: {len(neg_alns)}")

        for aln in neg_alns:
            if aln['pident'] >= args.min_identity and aln['query_cov'] >= args.min_coverage:
                qname = aln['query']
                if qname not in neg_hits or aln['pident'] > neg_hits[qname]['pident']:
                    neg_hits[qname] = aln
        print(f"  Bacterial hits: {len(neg_hits)}")

    # --- Assign labels -------------------------------------------------------
    results = []
    for contig_id in sorted(contig_lengths.keys()):
        length = contig_lengths[contig_id]

        if contig_id in viral_hits:
            hit = viral_hits[contig_id]
            target_acc = hit['target'].split()[0]
            category = infer_category(hit['target'])

            # Check if this genome is excluded (chimeric risk from homology check)
            is_excluded = any(target_acc.startswith(ex) for ex in excluded_genomes)
            label = 'excluded_chimeric_risk' if is_excluded else 'viral'

            results.append({
                'contig_id': contig_id,
                'label': label,
                'source_genome': hit['target'],
                'source_category': category,
                'pident': f"{hit['pident']:.1f}",
                'query_coverage': f"{hit['query_cov']:.1f}",
                'length': length,
            })
        elif contig_id in neg_hits:
            hit = neg_hits[contig_id]
            results.append({
                'contig_id': contig_id,
                'label': 'bacterial',
                'source_genome': hit['target'],
                'source_category': 'negative_control',
                'pident': f"{hit['pident']:.1f}",
                'query_coverage': f"{hit['query_cov']:.1f}",
                'length': length,
            })
        else:
            results.append({
                'contig_id': contig_id,
                'label': 'background',
                'source_genome': '-',
                'source_category': '-',
                'pident': '-',
                'query_coverage': '-',
                'length': length,
            })

    # --- Write output --------------------------------------------------------
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=[
            'contig_id', 'label', 'source_genome', 'source_category',
            'pident', 'query_coverage', 'length'
        ])
        writer.writeheader()
        writer.writerows(results)

    # --- Summary -------------------------------------------------------------
    from collections import Counter
    label_counts = Counter(r['label'] for r in results)
    cat_counts = Counter(r['source_category'] for r in results if r['label'] == 'viral')

    print(f"\n{'='*60}")
    print(f"Spike-In Ground Truth Summary")
    print(f"{'='*60}")
    print(f"  Total contigs:  {len(results)}")
    print(f"\n  By label:")
    for label in ['viral', 'bacterial', 'background', 'excluded_chimeric_risk']:
        count = label_counts.get(label, 0)
        if count > 0 or label != 'excluded_chimeric_risk':
            print(f"    {label:<25} {count:>6}")
    if cat_counts:
        print(f"\n  Viral by category:")
        for cat, count in sorted(cat_counts.items()):
            print(f"    {cat:<25} {count:>6}")
    print(f"\n  Output: {args.output}")


if __name__ == '__main__':
    main()
