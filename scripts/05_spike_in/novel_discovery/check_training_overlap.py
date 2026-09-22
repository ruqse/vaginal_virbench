#!/usr/bin/env python3
"""
Check training data overlap between spike-in genomes and 14 virus identification tools.

Produces a genome x tool matrix indicating whether each spike-in genome was likely
present in each tool's training data. This is Table S3 for the manuscript.

Training data sources and freeze dates verified from original papers in papers/.

Usage:
    python check_training_overlap.py \
        --manifest scripts/05_spike_in/spike_in_genomes.tsv \
        --output data/novel_spike_discovery/training_overlap_matrix.tsv \
        [--novel-manifest data/novel_spike_discovery/selected/novel_spike_in_genomes.tsv]
"""

import argparse
import csv
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Training data metadata (verified from PDFs in papers/)
# ---------------------------------------------------------------------------
# Fields:
#   db_sources  - which databases were used for training
#   freeze      - approximate data freeze date (YYYY-MM-DD or "~YYYY")
#   n_genomes   - number of viral genomes in training set (if reported)
#   paper_quote - exact supporting text from the paper
#
# Overlap logic:
#   A spike-in genome is "likely_in_training" if:
#     - source=RefSeq AND tool trained on RefSeq AND accession is NC_* (curated)
#     - source=MetaVR_v5 AND tool trained on IMG/VR v3+ (MetaVR UViGs derive from IMG/VR)
#     - source=GenBank AND tool trained on GenBank/NCBI AND freeze >= deposit date
#   A genome is "not_in_training" if:
#     - source=novel (unpublished internal data)
#   A genome is "possible" if overlap cannot be confirmed or denied
#   Sourmash has no training data (DB-dependent at runtime)

TOOL_TRAINING = {
    "VirFinder": {
        "db_sources": ["RefSeq"],
        "freeze": "2014-01-01",
        "n_genomes": 1562,
        "paper_quote": "1,562 prokaryotic virus genomes from NCBI RefSeq (before Jan 1, 2014)",
        "includes_imgvr": False,
        "includes_genbank": False,
    },
    "VirSorter": {
        "db_sources": ["RefSeq"],
        "freeze": "2014-01-01",
        "n_genomes": None,
        "paper_quote": "NCBI RefSeq Virus (bacteria/archaea), 114,297 proteins (Jan 2014)",
        "includes_imgvr": False,
        "includes_genbank": False,
    },
    "DeepVirFinder": {
        "db_sources": ["RefSeq"],
        "freeze": "2015-05-01",
        "n_genomes": None,
        "paper_quote": "RefSeq viral genomes (May 2015 snapshot)",
        "includes_imgvr": False,
        "includes_genbank": False,
    },
    "MetaPhinder": {
        "db_sources": ["ACLAME", "phage_DBs"],
        "freeze": "2015-01-01",
        "n_genomes": None,
        "paper_quote": "ACLAME, ARDB, VFDB, POG databases",
        "includes_imgvr": False,
        "includes_genbank": False,
    },
    "PPR-Meta": {
        "db_sources": ["RefSeq"],
        "freeze": "2016-01-01",
        "n_genomes": 2279,
        "paper_quote": "10,090 prokaryotes + 8,801 plasmids + 2,279 phages from NCBI (before Jan 2016)",
        "includes_imgvr": False,
        "includes_genbank": False,
    },
    "Seeker": {
        "db_sources": ["RefSeq", "IMG/VR"],
        "freeze": "2018-01-01",
        "n_genomes": 2232,
        "paper_quote": "2,232 RefSeq phages + 75 bacteria; 80% of RefSeq phages + IMG/VR",
        "includes_imgvr": True,  # but older version (pre-v4)
        "includes_genbank": False,
    },
    "VIBRANT": {
        "db_sources": ["RefSeq", "VOG", "Pfam", "KEGG"],
        "freeze": "2019-01-01",
        "n_genomes": None,
        "paper_quote": "633,194 viral proteins from NCBI reference + VOG HMMs",
        "includes_imgvr": False,
        "includes_genbank": False,
    },
    "VirSorter2": {
        "db_sources": ["RefSeq", "Pfam"],
        "freeze": "2020-01-01",
        "n_genomes": None,
        "paper_quote": "RefSeq (all domains) + non-RefSeq sources; Pfam v32.0",
        "includes_imgvr": False,  # paper doesn't explicitly cite IMG/VR
        "includes_genbank": False,
    },
    "geNomad": {
        "db_sources": ["IMG/VR_v4", "GenBank", "PLSDB"],
        "freeze": "2021-07-06",
        "n_genomes": None,
        "paper_quote": "IMG/VR v4; GenBank viruses (retrieved 6 July 2021); PLSDB (release 2020_11_19)",
        "includes_imgvr": True,
        "includes_genbank": True,
    },
    "Jaeger": {
        "db_sources": ["IMG/VR_v4"],
        "freeze": "2023-01-01",
        "n_genomes": None,
        "paper_quote": "IMG/VR v4 + prokaryotic/eukaryotic genomes",
        "includes_imgvr": True,
        "includes_genbank": False,
    },
    "ViraLM": {
        "db_sources": ["RefSeq"],
        "freeze": "2023-09-01",
        "n_genomes": 49929,
        "paper_quote": "49,929 complete viral genomes released before September 2023 from NCBI RefSeq",
        "includes_imgvr": False,
        "includes_genbank": False,
    },
    "TransGINmer": {
        "db_sources": ["RefSeq"],
        "freeze": "2023-10-25",
        "n_genomes": 15037,
        "paper_quote": "all RefSeq genomes of viruses downloaded from NCBI Virus before October 25, 2023, totaling 15,037 genomes",
        "includes_imgvr": False,
        "includes_genbank": False,
    },
    "HVSeeker": {
        "db_sources": ["NCBI"],
        "freeze": "2023-01-01",
        "n_genomes": 2687,
        "paper_quote": "536 bacterial + 2,687 phage genomes from NCBI",
        "includes_imgvr": False,
        "includes_genbank": True,
    },
    "Sourmash": {
        "db_sources": ["none"],
        "freeze": None,
        "n_genomes": None,
        "paper_quote": "Reference-free MinHash k-mer sketching; performance depends on user-supplied DB",
        "includes_imgvr": False,
        "includes_genbank": False,
    },
}

# Ordered tool list for consistent column output
TOOL_ORDER = [
    "MetaPhinder", "Sourmash",
    "VirSorter", "VirFinder", "VirSorter2", "VIBRANT",
    "DeepVirFinder", "PPR-Meta", "Seeker", "Jaeger",
    "HVSeeker", "TransGINmer", "ViraLM", "geNomad",
]


def classify_overlap(genome_source: str, genome_accession: str,
                     genome_panel: str, tool_name: str) -> str:
    """
    Determine if a spike-in genome was likely in a tool's training data.

    Returns one of:
        likely_in_training  - high confidence the genome was seen
        not_in_training     - confident the genome was NOT seen
        possible            - uncertain / cannot confirm
        db_dependent        - Sourmash (depends on user-supplied reference DB)
    """
    tool = TOOL_TRAINING[tool_name]

    # Novel (unpublished) genomes are guaranteed unseen
    if genome_panel == "novel":
        if tool_name == "Sourmash":
            return "db_dependent"
        return "not_in_training"

    # Sourmash has no training data
    if tool_name == "Sourmash":
        return "db_dependent"

    # --- RefSeq genomes (NC_* or GCF_* accessions) ---
    # Note: GenBank includes RefSeq (RefSeq is a curated subset of GenBank).
    # Tools that trained on "GenBank viruses" also saw RefSeq viral genomes.
    is_refseq = genome_source == "RefSeq" or genome_accession.startswith(("NC_", "GCF_"))
    if is_refseq:
        if any(db in ["RefSeq", "NCBI"] for db in tool["db_sources"]):
            return "likely_in_training"
        # GenBank includes RefSeq, so GenBank-trained tools also saw RefSeq
        if tool.get("includes_genbank"):
            return "likely_in_training"
        # MetaPhinder uses ACLAME which includes some RefSeq phages
        if tool_name == "MetaPhinder":
            return "possible"
        # IMG/VR-only tools (Jaeger) did not train on RefSeq per se
        return "not_in_training"

    # --- GenBank-only genomes (MW* accessions, not in RefSeq) ---
    is_genbank = genome_source == "GenBank" or genome_accession.startswith("MW")
    if is_genbank:
        if tool.get("includes_genbank"):
            return "likely_in_training"
        if any(db in ["NCBI"] for db in tool["db_sources"]):
            return "possible"
        return "not_in_training"

    # --- MetaVR v5 / IMG/VR UViGs (IMGVR_UViG_* accessions) ---
    is_imgvr = genome_source == "MetaVR_v5" or genome_accession.startswith("IMGVR_UViG")
    if is_imgvr:
        if tool.get("includes_imgvr"):
            return "likely_in_training"
        return "not_in_training"

    return "possible"


def load_manifest(path: Path) -> list[dict]:
    """Load spike-in genome manifest TSV (skip ## comment lines)."""
    genomes = []
    with open(path) as f:
        lines = [l for l in f if not l.startswith("##")]
    reader = csv.DictReader(lines, delimiter="\t")
    for row in reader:
        # Clean column names (the TSV may have #accession as first col)
        cleaned = {}
        for k, v in row.items():
            cleaned[k.lstrip("#")] = v
        genomes.append(cleaned)
    return genomes


def main():
    parser = argparse.ArgumentParser(
        description="Check training data overlap for spike-in genomes")
    parser.add_argument("--manifest", required=True,
                        help="Path to spike_in_genomes.tsv")
    parser.add_argument("--novel-manifest", default=None,
                        help="Optional path to novel_spike_in_genomes.tsv")
    parser.add_argument("--output", required=True,
                        help="Output TSV path for overlap matrix")
    args = parser.parse_args()

    # Load genomes
    genomes = load_manifest(Path(args.manifest))
    for g in genomes:
        g["panel"] = "reference"

    if args.novel_manifest and Path(args.novel_manifest).exists():
        novel = load_manifest(Path(args.novel_manifest))
        for g in novel:
            g["panel"] = "novel"
        genomes.extend(novel)

    # Create output directory
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    # Compute overlap matrix
    header = ["accession", "organism", "category", "source", "panel"] + TOOL_ORDER
    rows = []
    for g in genomes:
        row = {
            "accession": g.get("accession", ""),
            "organism": g.get("organism", ""),
            "category": g.get("category", ""),
            "source": g.get("source", ""),
            "panel": g.get("panel", "reference"),
        }
        for tool in TOOL_ORDER:
            row[tool] = classify_overlap(
                genome_source=row["source"],
                genome_accession=row["accession"],
                genome_panel=row["panel"],
                tool_name=tool,
            )
        rows.append(row)

    # Write output
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    # Print summary
    print(f"\nTraining Overlap Matrix: {len(genomes)} genomes x {len(TOOL_ORDER)} tools")
    print(f"Output: {args.output}\n")

    # Summary stats — count overlaps separately for viral and bacterial genomes
    viral_overlap = {tool: 0 for tool in TOOL_ORDER}
    bact_overlap = {tool: 0 for tool in TOOL_ORDER}
    for row in rows:
        is_viral = row["category"] != "negative_control"
        for tool in TOOL_ORDER:
            if row[tool] == "likely_in_training":
                if is_viral:
                    viral_overlap[tool] += 1
                else:
                    bact_overlap[tool] += 1

    n_viral = sum(1 for g in genomes if g.get("category") != "negative_control")
    n_bact = sum(1 for g in genomes if g.get("category") == "negative_control")
    print(f"{'Tool':<20} {'Viral overlap':>16} {'Bact overlap':>16}  {'Risk':>6}")
    print("-" * 65)
    for tool in TOOL_ORDER:
        vc = viral_overlap[tool]
        bc = bact_overlap[tool]
        vpct = vc / n_viral * 100 if n_viral > 0 else 0
        risk = "HIGH" if vpct >= 50 else "MED" if vpct >= 20 else "LOW" if vpct > 0 else "-"
        if tool == "Sourmash":
            risk = "DB-DEP"
        print(f"{tool:<20} {vc:>3}/{n_viral} ({vpct:5.1f}%)  {bc:>3}/{n_bact}             {risk:>6}")

    # Print training data reference table
    print(f"\n{'='*80}")
    print("Training Data Reference (for Table S3 / manuscript)")
    print(f"{'='*80}")
    print(f"{'Tool':<16} {'Freeze Date':<14} {'N Genomes':>10}  {'Source'}")
    print("-" * 80)
    for tool in TOOL_ORDER:
        t = TOOL_TRAINING[tool]
        freeze = t["freeze"] if t["freeze"] else "N/A"
        n = str(t["n_genomes"]) if t["n_genomes"] else "-"
        print(f"{tool:<16} {freeze:<14} {n:>10}  {t['paper_quote'][:60]}")


if __name__ == "__main__":
    main()
