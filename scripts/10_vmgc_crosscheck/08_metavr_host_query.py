#!/usr/bin/env python3
"""Count MetaVR host annotations for a selected set of vaginal taxa.

Backs the MetaVR host-taxonomy count in the Discussion and the Supplementary
Methods section "MetaVR host-taxonomy query and its interpretation".

Rule (a selected set, not specified before the export was inspected):
  exact species assignment to Lactobacillus crispatus, L. gasseri, L. iners or
  L. jensenii, or genus assignment to Atopobium, Fannyhessea, Mobiluncus or
  Sneathia. Denominator: records with a populated host_taxonomy field.
Sensitivity rules: the whole Lactobacillus genus plus the four genera, and that
set plus Prevotella, Megasphaera and Dialister. Bifidobacterium species counts
are kept because GTDB places the Gardnerella clade within Bifidobacterium and
many Bifidobacterium assignments lack a species.

Counts host annotations only. It cannot estimate how completely vaginal viral
sequences are covered, or a host-specific deficit in tool sensitivity.

Input : databases/metavr_v5/DownloadUvigMetadata.tsv.gz (local copy dated 2026-02-23)
Output: scripts/10_vmgc_crosscheck/08_metavr_host_counts.json
Usage : python3 scripts/10_vmgc_crosscheck/08_metavr_host_query.py

Adapted without logic changes from a read-only recount written during manuscript review.
"""
import collections
import datetime
import gzip
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "databases/metavr_v5/DownloadUvigMetadata.tsv.gz"
OUT = Path(__file__).with_name("08_metavr_host_counts.json")
SPECIES = {"Lactobacillus crispatus", "Lactobacillus iners", "Lactobacillus gasseri", "Lactobacillus jensenii"}
GENERA = {"Atopobium", "Mobiluncus", "Fannyhessea", "Sneathia"}
EXPANDED = GENERA | {"Lactobacillus", "Prevotella", "Megasphaera", "Dialister"}
MISSING = {"", "\\N", "NA", "None"}


def main():
    counts = collections.Counter()
    genus_counts = collections.Counter()
    species_counts = collections.Counter()
    bifidobacterium_species = collections.Counter()
    host_methods = collections.Counter()
    with gzip.open(SOURCE, "rt") as handle:
        header = next(handle).rstrip("\n").split("\t")
        hi = header.index("host_taxonomy")
        mi = header.index("host_taxonomy_method")
        for line in handle:
            counts["all_records"] += 1
            cells = line.rstrip("\n").split("\t")
            host = cells[hi]
            if host in MISSING:
                continue
            counts["populated_host_taxonomy"] += 1
            host_methods[cells[mi]] += 1
            ranks = dict(rank.split("__", 1) for rank in host.split(";") if "__" in rank)
            genus, species = ranks.get("g", ""), ranks.get("s", "")
            if genus in EXPANDED:
                genus_counts[genus] += 1
            if genus == "Lactobacillus":
                species_counts[species] += 1
            if genus == "Bifidobacterium":
                bifidobacterium_species[species] += 1
            if "gardnerella" in host.lower():
                counts["literal_gardnerella_records"] += 1
            if species in SPECIES:
                counts["four_lactobacillus_species"] += 1
            if genus in GENERA:
                counts["four_other_genera"] += 1
            if species in SPECIES or genus in GENERA:
                counts["narrow_rule"] += 1
            if genus in GENERA | {"Lactobacillus"}:
                counts["five_genus_rule"] += 1
            if genus in EXPANDED:
                counts["eight_genus_rule"] += 1
    counts.setdefault("literal_gardnerella_records", 0)
    result = {
        "source": str(SOURCE.relative_to(ROOT)),
        "source_size_bytes": SOURCE.stat().st_size,
        "source_mtime_utc_not_download_provenance": datetime.datetime.fromtimestamp(
            SOURCE.stat().st_mtime, datetime.timezone.utc).isoformat(),
        "species_rule": sorted(SPECIES),
        "other_genus_rule": sorted(GENERA),
        "counts": dict(counts),
        "selected_genus_counts": dict(sorted(genus_counts.items())),
        "lactobacillus_species_counts": dict(sorted(species_counts.items())),
        "bifidobacterium_species_counts": dict(sorted(bifidobacterium_species.items())),
        "host_taxonomy_methods": dict(sorted(host_methods.items())),
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"counts": result["counts"], "selected_genus_counts": result["selected_genus_counts"]}, indent=2))
    print(f"Full results: {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
