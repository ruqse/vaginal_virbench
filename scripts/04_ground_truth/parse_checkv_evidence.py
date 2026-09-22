#!/usr/bin/env python3
"""
parse_checkv_evidence.py — Evidence 2: Structural Signatures from CheckV

Parses CheckV quality_summary.tsv output to extract structural evidence
for viral classification of contigs.

Trigger conditions (strategy §4.2):
  - checkv_quality in {Medium-quality, High-quality, Complete}
  - provirus == "Yes" (DTR/ITR indicator)
  - viral_genes >= 1

Output: evidence_2_structural.tsv

Usage:
    python parse_checkv_evidence.py \
        --checkv-dir results/test_real/checkv/ \
        --output results/test_real/ground_truth/evidence_2_structural.tsv
"""

import argparse
import csv
import os
import sys
import glob


def parse_quality_summary(filepath):
    """Parse a single CheckV quality_summary.tsv file."""
    records = []
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            records.append(row)
    return records


def evaluate_evidence(record):
    """
    Determine whether a contig triggers Evidence 2 (structural signatures).

    Trigger conditions:
      1. checkv_quality in {Medium-quality, High-quality, Complete}
      2. provirus == "Yes"
      3. viral_genes >= 1

    Returns (triggered: bool, reasons: list[str])
    """
    reasons = []

    # Condition 1: CheckV quality
    quality = record.get('checkv_quality', 'Not-determined')
    quality_trigger = quality in ('Medium-quality', 'High-quality', 'Complete')
    if quality_trigger:
        reasons.append(f'quality:{quality}')

    # Condition 2: Provirus detection (DTR/ITR)
    provirus = record.get('provirus', 'No')
    provirus_trigger = provirus.lower() == 'yes'
    if provirus_trigger:
        reasons.append('provirus:Yes')

    # Condition 3: Viral genes
    try:
        viral_genes = int(record.get('viral_genes', '0'))
    except (ValueError, TypeError):
        viral_genes = 0
    viral_genes_trigger = viral_genes >= 1
    if viral_genes_trigger:
        reasons.append(f'viral_genes:{viral_genes}')

    triggered = quality_trigger or provirus_trigger or viral_genes_trigger
    return triggered, reasons


def find_quality_summaries(checkv_dir):
    """Find all quality_summary.tsv files in CheckV output directory."""
    patterns = [
        os.path.join(checkv_dir, '**', 'quality_summary.tsv'),
        os.path.join(checkv_dir, '*_checkv', 'quality_summary.tsv'),
        os.path.join(checkv_dir, 'quality_summary.tsv'),
    ]

    found = set()
    for pattern in patterns:
        for f in glob.glob(pattern, recursive=True):
            found.add(os.path.abspath(f))

    return sorted(found)


def main():
    parser = argparse.ArgumentParser(
        description='Parse CheckV output for Evidence 2 (structural signatures)'
    )
    parser.add_argument(
        '--checkv-dir',
        required=True,
        help='Path to CheckV output directory'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output evidence TSV path'
    )
    args = parser.parse_args()

    # Find quality_summary.tsv files
    summaries = find_quality_summaries(args.checkv_dir)

    if not summaries:
        print(f"WARNING: No quality_summary.tsv files found in {args.checkv_dir}",
              file=sys.stderr)
        print("  CheckV may not have been run yet.", file=sys.stderr)
        print("  Creating empty evidence file with header only.", file=sys.stderr)

        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        with open(args.output, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow([
                'contig_id', 'checkv_quality', 'completeness',
                'provirus', 'viral_genes', 'host_genes',
                'contig_length', 'trigger_reasons', 'evidence_triggered'
            ])
        print(f"  Empty evidence file written to: {args.output}")
        sys.exit(0)

    print(f"Found {len(summaries)} quality_summary.tsv file(s):")
    for s in summaries:
        print(f"  {s}")

    # Parse all records
    all_records = []
    for summary_file in summaries:
        records = parse_quality_summary(summary_file)
        all_records.extend(records)
        print(f"  Parsed {len(records)} records from {os.path.basename(os.path.dirname(summary_file))}")

    # Deduplicate by contig_id (keep first occurrence)
    seen = set()
    unique_records = []
    for rec in all_records:
        cid = rec.get('contig_id', rec.get('seqname', ''))
        if cid and cid not in seen:
            seen.add(cid)
            unique_records.append(rec)

    print(f"\nTotal unique contigs: {len(unique_records)}")

    # Evaluate evidence for each contig
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    n_triggered = 0
    n_total = 0

    with open(args.output, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            'contig_id', 'checkv_quality', 'completeness',
            'provirus', 'viral_genes', 'host_genes',
            'contig_length', 'trigger_reasons', 'evidence_triggered'
        ])

        for rec in unique_records:
            contig_id = rec.get('contig_id', rec.get('seqname', ''))
            if not contig_id:
                continue

            n_total += 1
            triggered, reasons = evaluate_evidence(rec)
            if triggered:
                n_triggered += 1

            writer.writerow([
                contig_id,
                rec.get('checkv_quality', 'Not-determined'),
                rec.get('completeness', 'NA'),
                rec.get('provirus', 'No'),
                rec.get('viral_genes', '0'),
                rec.get('host_genes', '0'),
                rec.get('contig_length', '0'),
                ';'.join(reasons) if reasons else '',
                'TRUE' if triggered else 'FALSE'
            ])

    # Summary
    print(f"\nEvidence 2 Summary:")
    print(f"  Total contigs:    {n_total}")
    print(f"  Evidence triggered: {n_triggered}")
    print(f"  Not triggered:     {n_total - n_triggered}")
    print(f"\nOutput: {args.output}")

    # Quality distribution
    quality_dist = {}
    for rec in unique_records:
        q = rec.get('checkv_quality', 'Not-determined')
        quality_dist[q] = quality_dist.get(q, 0) + 1
    print("\nCheckV Quality Distribution:")
    for q in sorted(quality_dist.keys()):
        print(f"  {q}: {quality_dist[q]}")


if __name__ == '__main__':
    main()
