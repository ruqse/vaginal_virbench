#!/usr/bin/env python3
"""
build_rca_ground_truth.py — RCA Benchmark Ground Truth from Co-Assembly + Enrichment

Constructs the gold-standard ground truth for the RCA benchmark (Track C) by
combining three independent evidence sources:

  1. Enrichment ratio: Library-size-normalized ratio of RCA/shotgun coverage
     per co-assembled "master contig" (from compute_enrichment.sh Phase 2)
  2. CheckV structural: Viral hallmark genes, completeness, provirus detection
  3. Kraken2 taxonomy: Bacterial/archaeal classification for TN confirmation

Classification logic:
  TP (Gold Standard viral):
    - Enriched (R > 10) AND CheckV viral signal
      (viral_genes > 0 OR quality >= Medium OR provirus = Yes)
  TN (True Negative bacterial):
    - Background (0.5 < R < 2) AND Kraken2 bacterial AND CheckV viral_genes = 0
  Dark Matter:
    - Enriched (R > 10) AND no CheckV signal AND no host genes AND unclassified
  Excluded:
    - Everything else (ambiguous enrichment, conflicting signals)

Output: TSV with columns:
  contig_id, patient, contig_length, enrichment_class, enrichment_R,
  checkv_quality, viral_genes, host_genes, provirus, kraken2_taxon,
  kraken2_domain, label, reason

Usage:
    # Phase 3a: Pool master contigs (run first)
    python build_rca_ground_truth.py pool-contigs \
        --coassembly-dir results/test_real/coassembly \
        --output results/test_real/coassembly/all_master_contigs.fasta

    # Phase 3d: Build ground truth (after CheckV + Kraken2)
    python build_rca_ground_truth.py build \
        --coassembly-dir results/test_real/coassembly \
        --checkv-dir results/test_real/coassembly/checkv \
        --kraken2-output results/test_real/coassembly/kraken2/kraken2_output.tsv \
        --output results/test_real/coassembly/rca_ground_truth.tsv

See the Supplementary Methods for the full ground truth construction pipeline.
"""

import argparse
import csv
import os
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


# =============================================================================
# Constants
# =============================================================================

PATIENTS = [
    "UC028_V2", "UC055_V1", "UC055_V2", "UC062_V2", "UC065_V2",
    "UC074_V2", "UC084_V2", "UC093_V2", "UC093_V3", "UC096_V2",
    "UC115_V2", "UC139_V2", "UC164_V2",
]

# Enrichment thresholds (must match compute_enrichment.sh)
ENRICHED_MIN = 10.0
BACKGROUND_LOW = 0.5
BACKGROUND_HIGH = 2.0

# Bacterial domain keywords for Kraken2 taxonomy
BACTERIAL_DOMAINS = {"Bacteria", "Archaea"}


# =============================================================================
# Pool Contigs (Phase 3a)
# =============================================================================

def pool_contigs(coassembly_dir: Path, output: Path):
    """Concatenate all patient master contigs with patient-prefixed IDs.

    Prefix format: {PATIENT}__ (double underscore, matching existing convention
    from prepare_rca_pooled.sh).
    """
    total_contigs = 0
    patients_found = 0

    with open(output, "w") as out:
        for patient in sorted(PATIENTS):
            fasta = coassembly_dir / f"{patient}_master_contigs.fasta"
            if not fasta.exists():
                print(f"  WARN: Missing {fasta}", file=sys.stderr)
                continue

            n = 0
            seq_lines = []
            current_header = None

            with open(fasta) as f:
                for line in f:
                    if line.startswith(">"):
                        if current_header is not None:
                            out.write(f">{patient}__{current_header}\n")
                            out.write("".join(seq_lines))
                            n += 1
                        current_header = line[1:].strip()
                        seq_lines = []
                    else:
                        seq_lines.append(line)

                # Last record
                if current_header is not None:
                    out.write(f">{patient}__{current_header}\n")
                    out.write("".join(seq_lines))
                    n += 1

            total_contigs += n
            patients_found += 1
            print(f"  {patient}: {n} contigs")

    print(f"\nPooled: {total_contigs} contigs from {patients_found} patients")
    print(f"Output: {output}")
    return total_contigs


# =============================================================================
# Load Enrichment Data (from Phase 2)
# =============================================================================

def load_enrichment(coassembly_dir: Path) -> dict:
    """Load all per-patient enrichment TSVs into unified dict.

    Returns {prefixed_contig_id: {enrichment_class, enrichment_R, contig_length, patient}}
    """
    enrichment = {}
    for patient in PATIENTS:
        tsv = coassembly_dir / f"{patient}_enrichment.tsv"
        if not tsv.exists():
            print(f"  WARN: Missing enrichment: {tsv}", file=sys.stderr)
            continue

        n = 0
        with open(tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                contig_id = row["contig_id"]
                prefixed_id = f"{patient}__{contig_id}"

                # Parse enrichment ratio (may be 'Inf')
                r_str = row["enrichment_ratio_normalized"]
                r_val = float("inf") if r_str == "Inf" else float(r_str)

                enrichment[prefixed_id] = {
                    "enrichment_class": row["classification"],
                    "enrichment_R": r_val,
                    "enrichment_R_raw": row["enrichment_ratio_raw"],
                    "contig_length": int(row["contig_length"]),
                    "shotgun_depth": float(row["shotgun_depth"]),
                    "rca_depth": float(row["rca_depth"]),
                    "patient": patient,
                }
                n += 1

        print(f"  {patient}: {n} contigs with enrichment data")

    return enrichment


# =============================================================================
# Load CheckV Data
# =============================================================================

def load_checkv(checkv_dir: Path) -> dict:
    """Load CheckV quality_summary.tsv.

    Returns {contig_id: {checkv_quality, viral_genes, host_genes, provirus, ...}}
    """
    quality_tsv = checkv_dir / "quality_summary.tsv"
    if not quality_tsv.exists():
        print(f"ERROR: CheckV output not found: {quality_tsv}", file=sys.stderr)
        sys.exit(1)

    checkv = {}
    with open(quality_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cid = row["contig_id"]
            checkv[cid] = {
                "checkv_quality": row.get("checkv_quality", "Not-determined"),
                "viral_genes": int(row.get("viral_genes", 0)),
                "host_genes": int(row.get("host_genes", 0)),
                "provirus": row.get("provirus", "No"),
                "completeness": row.get("completeness", "NA"),
                "miuvig_quality": row.get("miuvig_quality", "NA"),
            }

    print(f"  CheckV loaded: {len(checkv)} contigs")
    return checkv


def has_viral_signal(checkv_row: dict) -> tuple:
    """Determine if contig has CheckV viral signal.

    Returns (bool, reason_string).
    Same logic as run_checkv_pooled.sh Step 3.
    """
    quality = checkv_row.get("checkv_quality", "Not-determined")
    viral_genes = checkv_row.get("viral_genes", 0)
    provirus = checkv_row.get("provirus", "No")

    if quality in ("Medium-quality", "High-quality", "Complete"):
        return True, f"checkv_quality={quality}"
    if viral_genes >= 1:
        return True, f"viral_genes={viral_genes}"
    if provirus == "Yes":
        return True, "provirus=Yes"
    return False, "no_viral_signal"


# =============================================================================
# Load Kraken2 Taxonomy
# =============================================================================

def load_kraken2(kraken2_output: Path) -> dict:
    """Load Kraken2 output TSV (classify format).

    Kraken2 output columns (tab-separated):
      C/U, contig_id, taxon_id (tax_id:taxon_name), length, LCA_mapping
    The --use-names flag puts 'taxon_name (taxid NNNNN)' in column 3.

    Returns {contig_id: {classified, taxon_name, domain}}
    """
    kraken = {}

    with open(kraken2_output) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 3:
                continue

            status = fields[0]  # C or U
            contig_id = fields[1]
            taxon_info = fields[2]  # e.g. "Escherichia coli (taxid 562)"

            # Parse taxon name (remove taxid suffix)
            taxon_name = taxon_info.split(" (taxid")[0].strip() if "(taxid" in taxon_info else taxon_info

            # Determine domain from taxon hierarchy
            # For simplicity, use the taxon name heuristics
            # A more rigorous approach would parse the Kraken2 report for lineage
            classified = status == "C"

            kraken[contig_id] = {
                "classified": classified,
                "taxon_name": taxon_name if classified else "unclassified",
            }

    print(f"  Kraken2 loaded: {len(kraken)} contigs "
          f"({sum(1 for v in kraken.values() if v['classified'])} classified)")
    return kraken


def load_kraken2_report(report_path: Path) -> dict:
    """Parse Kraken2 report.txt to get domain-level classification.

    Report format (tab-separated):
      pct_reads, n_reads_clade, n_reads_taxon, rank_code, taxid, name

    We extract domain (D) assignments and propagate to descendant taxa.
    Returns {taxid: domain_name}
    """
    taxid_to_domain = {}
    current_domain = None

    if not report_path.exists():
        return taxid_to_domain

    with open(report_path) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 6:
                continue
            rank = fields[3].strip()
            taxid = fields[4].strip()
            name = fields[5].strip()

            if rank == "D":
                current_domain = name
                taxid_to_domain[taxid] = name
            elif rank in ("D", "P", "C", "O", "F", "G", "S", "S1", "S2"):
                if current_domain:
                    taxid_to_domain[taxid] = current_domain

    return taxid_to_domain


def classify_domain(kraken2_output: Path, kraken2_report: Path = None) -> dict:
    """Classify each contig into domain using Kraken2 output + report.

    Falls back to name-based heuristics if report not available.
    Returns {contig_id: domain_string}
    """
    # Try to use report for accurate domain mapping
    report_path = kraken2_report
    if report_path is None:
        report_path = kraken2_output.parent / "kraken2_report.txt"

    taxid_domains = load_kraken2_report(report_path) if report_path.exists() else {}

    domains = {}
    with open(kraken2_output) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 3:
                continue

            status = fields[0]
            contig_id = fields[1]
            taxon_info = fields[2]

            if status == "U":
                domains[contig_id] = "unclassified"
                continue

            # Try taxid-based domain lookup
            taxid = None
            if "(taxid" in taxon_info:
                try:
                    taxid = taxon_info.split("(taxid ")[1].rstrip(")")
                except (IndexError, ValueError):
                    pass

            if taxid and taxid in taxid_domains:
                domains[contig_id] = taxid_domains[taxid]
            else:
                # Heuristic: check if taxon name suggests a viral assignment
                taxon_lower = taxon_info.lower()
                if any(kw in taxon_lower for kw in ["virus", "phage", "viridae", "virales"]):
                    domains[contig_id] = "Viruses"
                else:
                    # Default: most contigs in metagenome assemblies are bacterial
                    domains[contig_id] = "Bacteria"

    return domains


# =============================================================================
# Ground Truth Classification (Phase 3d)
# =============================================================================

def build_ground_truth(enrichment: dict, checkv: dict, kraken: dict,
                       domains: dict) -> list:
    """Classify contigs into TP, TN, DarkMatter, or Excluded.

    Returns list of dicts (one per contig).
    """
    results = []
    label_counts = Counter()

    for contig_id in sorted(enrichment.keys()):
        e = enrichment[contig_id]
        cv = checkv.get(contig_id, {
            "checkv_quality": "Not-determined",
            "viral_genes": 0,
            "host_genes": 0,
            "provirus": "No",
            "completeness": "NA",
            "miuvig_quality": "NA",
        })
        kr = kraken.get(contig_id, {"classified": False, "taxon_name": "unclassified"})
        domain = domains.get(contig_id, "unknown")

        is_viral, viral_reason = has_viral_signal(cv)
        enrich_class = e["enrichment_class"]
        enrich_R = e["enrichment_R"]

        # --- Classification logic ---
        label = "excluded"
        reason = ""

        if enrich_class == "enriched":
            if is_viral:
                # TP: Enriched + CheckV viral signal
                label = "TP"
                reason = f"enriched(R={enrich_R:.1f})+{viral_reason}"
            elif cv["viral_genes"] == 0 and cv["host_genes"] == 0 and not kr["classified"]:
                # Dark Matter: Enriched but unknown — potentially novel virus
                label = "dark_matter"
                reason = f"enriched(R={enrich_R:.1f})+no_checkv+unclassified"
            elif cv["host_genes"] > 0 and cv["viral_genes"] == 0:
                # Enriched bacterial — possibly mobile element, exclude
                label = "excluded"
                reason = f"enriched(R={enrich_R:.1f})+host_genes={cv['host_genes']}+no_viral"
            else:
                # Enriched but ambiguous signal
                label = "excluded"
                reason = f"enriched(R={enrich_R:.1f})+ambiguous({viral_reason})"

        elif enrich_class == "background":
            if not is_viral and domain in BACTERIAL_DOMAINS:
                # TN: Background + no viral signal + bacterial taxonomy
                label = "TN"
                reason = f"background(R={enrich_R:.2f})+{viral_reason}+kraken2={kr['taxon_name']}"
            elif is_viral:
                # Background but CheckV says viral — prophage? Exclude.
                label = "excluded"
                reason = f"background(R={enrich_R:.2f})+viral_signal({viral_reason})"
            else:
                # Background but non-bacterial taxonomy
                label = "excluded"
                reason = f"background(R={enrich_R:.2f})+domain={domain}"

        elif enrich_class == "moderate":
            label = "excluded"
            reason = f"moderate_enrichment(R={enrich_R:.2f})"

        elif enrich_class == "depleted":
            label = "excluded"
            reason = f"depleted(R={enrich_R:.2f})"

        label_counts[label] += 1

        results.append({
            "contig_id": contig_id,
            "patient": e["patient"],
            "contig_length": e["contig_length"],
            "enrichment_class": enrich_class,
            "enrichment_R": f"{enrich_R:.4f}" if enrich_R != float("inf") else "Inf",
            "shotgun_depth": e["shotgun_depth"],
            "rca_depth": e["rca_depth"],
            "checkv_quality": cv["checkv_quality"],
            "viral_genes": cv["viral_genes"],
            "host_genes": cv["host_genes"],
            "provirus": cv["provirus"],
            "completeness": cv.get("completeness", "NA"),
            "kraken2_taxon": kr["taxon_name"],
            "kraken2_domain": domain,
            "label": label,
            "reason": reason,
        })

    return results, label_counts


# =============================================================================
# Output
# =============================================================================

def write_ground_truth(results: list, output: Path):
    """Write ground truth TSV."""
    fieldnames = [
        "contig_id", "patient", "contig_length", "enrichment_class",
        "enrichment_R", "shotgun_depth", "rca_depth",
        "checkv_quality", "viral_genes", "host_genes", "provirus",
        "completeness", "kraken2_taxon", "kraken2_domain",
        "label", "reason",
    ]

    with open(output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


def print_summary(results: list, label_counts: Counter):
    """Print classification summary."""
    print(f"\n{'='*60}")
    print("RCA GROUND TRUTH SUMMARY (Co-Assembly + Enrichment Ratio)")
    print(f"{'='*60}\n")

    total = len(results)
    print(f"Total contigs: {total}\n")

    for label in ["TP", "TN", "dark_matter", "excluded"]:
        n = label_counts.get(label, 0)
        pct = 100 * n / total if total > 0 else 0
        print(f"  {label:12s}: {n:6d} ({pct:5.1f}%)")

    # TP breakdown by CheckV quality
    tp_contigs = [r for r in results if r["label"] == "TP"]
    if tp_contigs:
        print(f"\nTP breakdown by CheckV quality:")
        quality_counts = Counter(r["checkv_quality"] for r in tp_contigs)
        for q in ["Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"]:
            if quality_counts.get(q, 0) > 0:
                print(f"  {q}: {quality_counts[q]}")

    # TP breakdown by patient
    patient_counts = Counter(r["patient"] for r in tp_contigs)
    if patient_counts:
        print(f"\nTP per patient:")
        for patient in sorted(PATIENTS):
            n = patient_counts.get(patient, 0)
            if n > 0:
                print(f"  {patient}: {n}")

    # TN breakdown by domain
    tn_contigs = [r for r in results if r["label"] == "TN"]
    if tn_contigs:
        print(f"\nTN breakdown:")
        domain_counts = Counter(r["kraken2_domain"] for r in tn_contigs)
        for domain, count in domain_counts.most_common():
            print(f"  {domain}: {count}")

    # Dark matter
    dm_contigs = [r for r in results if r["label"] == "dark_matter"]
    if dm_contigs:
        print(f"\nDark matter: {len(dm_contigs)} contigs")
        lengths = [r["contig_length"] for r in dm_contigs]
        print(f"  Length range: {min(lengths):,} - {max(lengths):,} bp")
        print(f"  Mean length:  {sum(lengths)/len(lengths):,.0f} bp")

    # Class ratio
    n_tp = label_counts.get("TP", 0)
    n_tn = label_counts.get("TN", 0)
    if n_tp > 0:
        print(f"\nBenchmark class ratio (TP:TN): 1:{n_tn/n_tp:.1f}")
    print()


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="RCA Benchmark Ground Truth — Co-Assembly + Enrichment Ratio",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- pool-contigs subcommand ---
    pool_parser = subparsers.add_parser(
        "pool-contigs",
        help="Pool all patient master contigs into single FASTA (Phase 3a)",
    )
    pool_parser.add_argument(
        "--coassembly-dir", required=True, type=Path,
        help="Directory with {PATIENT}_master_contigs.fasta files",
    )
    pool_parser.add_argument(
        "--output", required=True, type=Path,
        help="Output pooled FASTA (e.g., all_master_contigs.fasta)",
    )

    # --- build subcommand ---
    build_parser = subparsers.add_parser(
        "build",
        help="Build ground truth from enrichment + CheckV + Kraken2 (Phase 3d)",
    )
    build_parser.add_argument(
        "--coassembly-dir", required=True, type=Path,
        help="Directory with {PATIENT}_enrichment.tsv and master contigs",
    )
    build_parser.add_argument(
        "--checkv-dir", required=True, type=Path,
        help="CheckV output directory (must contain quality_summary.tsv)",
    )
    build_parser.add_argument(
        "--kraken2-output", required=True, type=Path,
        help="Kraken2 output TSV (from kraken2 --output)",
    )
    build_parser.add_argument(
        "--kraken2-report", type=Path, default=None,
        help="Kraken2 report TXT (for domain classification; auto-detected if omitted)",
    )
    build_parser.add_argument(
        "--output", required=True, type=Path,
        help="Output ground truth TSV",
    )
    build_parser.add_argument(
        "--enriched-threshold", type=float, default=ENRICHED_MIN,
        help=f"Minimum normalized R for 'enriched' class (default: {ENRICHED_MIN})",
    )

    args = parser.parse_args()

    if args.command == "pool-contigs":
        print(f"[{datetime.now():%H:%M:%S}] Pooling master contigs...")
        os.makedirs(args.output.parent, exist_ok=True)
        n = pool_contigs(args.coassembly_dir, args.output)
        if n == 0:
            print("ERROR: No contigs pooled. Check --coassembly-dir.", file=sys.stderr)
            sys.exit(1)

        # Also copy to spades dir for tool execution
        spades_dir = args.coassembly_dir.parent / "spades"
        if spades_dir.exists():
            import shutil
            dest = spades_dir / "RCA_master_contigs.fasta"
            shutil.copy2(args.output, dest)
            print(f"\nCopied to: {dest}")
            print("  Ready for tool execution:")
            print(f"  sbatch run_secondary_cpu_tools.sh RCA_master")
            print(f"  sbatch run_secondary_gpu_tools.sh RCA_master")

    elif args.command == "build":
        print(f"[{datetime.now():%H:%M:%S}] Building RCA ground truth...\n")

        # Note: enrichment classes (enriched/background/moderate/depleted) are
        # pre-computed by compute_enrichment.sh. The --enriched-threshold flag
        # is stored for reporting but does not reclassify contigs here.
        enriched_threshold = args.enriched_threshold

        # Load all data sources
        print("Loading enrichment data:")
        enrichment = load_enrichment(args.coassembly_dir)
        if not enrichment:
            print("ERROR: No enrichment data loaded.", file=sys.stderr)
            sys.exit(1)

        print("\nLoading CheckV data:")
        checkv = load_checkv(args.checkv_dir)

        print("\nLoading Kraken2 data:")
        kraken = load_kraken2(args.kraken2_output)
        domains = classify_domain(args.kraken2_output, args.kraken2_report)

        # Build ground truth
        print(f"\nClassifying {len(enrichment)} contigs...")
        results, label_counts = build_ground_truth(enrichment, checkv, kraken, domains)

        # Write output
        os.makedirs(args.output.parent, exist_ok=True)
        write_ground_truth(results, args.output)
        print(f"\nGround truth written to: {args.output}")

        # Summary
        print_summary(results, label_counts)

        # Write label-specific contig lists (for filtering)
        for label in ["TP", "TN", "dark_matter"]:
            ids = [r["contig_id"] for r in results if r["label"] == label]
            list_file = args.output.parent / f"rca_{label.lower()}_contig_ids.txt"
            with open(list_file, "w") as f:
                for cid in ids:
                    f.write(f"{cid}\n")
            print(f"  {label} IDs: {list_file} ({len(ids)} contigs)")


if __name__ == "__main__":
    main()
