#!/usr/bin/env python3
"""
04_vmgc_coverage_table.py — Part 2 analysis (VMGC benchmark-reference overlap).

For each of the 4,263 VMGC vOTU reps, classify presence in:
  * MetaVR v5 (nucleotide, species-level via the corrected coverage parser):
      present_metavr_species   : ANI >= 95 % over >= 85 % AF of the shorter seq
      present_metavr_relative_only : a hit below the species threshold
      absent_metavr            : no hit (>=90 % perc_identity pre-filter; <90 %
                                 relatives are not captured -- documented limit)
  * RefSeq viral protein (DIAMOND blastx, protein homology present/absent):
      refseq_protein_homolog / no_refseq_protein_homolog  (MIN_ALEN >= 50 aa)

Reported as PROJECT-SPECIFIC absence from MetaVR/RefSeq -- explicitly NOT a
re-derivation of Huang's "85.8 % absent from five databases" (different DB set).

Joins VMGC_virus.info for native family + host (phylum/species) stratification.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vmgc_coverage import best_subject_per_query  # noqa: E402

SPECIES_ANI = 95.0
SPECIES_AF_SHORTER = 85.0
REFSEQ_MIN_ALEN = 50  # aa


def read_fasta_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def load_vmgc_info(path):
    info = {}
    with open(path, newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        fam_c = "Taxonomic_classification_(family_level)"
        php_c = "Host_assignment_(phylum_level)"
        hsp_c = "Host_assignment_(species_level)"
        for row in r:
            info[row["vOTU_ID"]] = {
                "family": row.get(fam_c, "") or "Unclassified",
                "host_phylum": row.get(php_c, "") or "-",
                "host_species": row.get(hsp_c, "") or "-",
            }
    return info


def host_genus(host_species):
    """Representative genus = genus of the first listed species."""
    if not host_species or host_species == "-":
        return "-"
    first = host_species.split(",")[0].strip()
    return first.split()[0] if first else "-"


def is_gardnerella(host_species):
    s = (host_species or "").lower()
    return ("gardnerella" in s) or ("bifidobacterium vaginale" in s)


def best_refseq_protein(refseq_blast):
    """vOTU -> best protein hit dict (alen>=50aa), keyed by query (vOTU_ID).
    DIAMOND outfmt6: qseqid sseqid pident length ... evalue bitscore qlen slen stitle
    """
    best = {}
    if not os.path.isfile(refseq_blast):
        return best
    with open(refseq_blast) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 14:
                continue
            q = f[0]
            pident = float(f[2])
            alen = int(f[3])          # aa
            bitscore = float(f[11])
            stitle = f[14] if len(f) > 14 else ""
            if alen < REFSEQ_MIN_ALEN:
                continue
            if q not in best or bitscore > best[q]["bitscore"]:
                best[q] = {"sseqid": f[1], "pident": pident, "alen": alen,
                           "bitscore": bitscore, "stitle": stitle}
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metavr-blast", required=True)
    ap.add_argument("--refseq-blast", required=True)
    ap.add_argument("--votu-fasta", required=True)
    ap.add_argument("--vmgc-info", required=True)
    ap.add_argument("--tables-dir", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    votu_ids = read_fasta_ids(a.votu_fasta)
    info = load_vmgc_info(a.vmgc_info)
    metavr_best = best_subject_per_query(a.metavr_blast)
    refseq_best = best_refseq_protein(a.refseq_blast)

    os.makedirs(a.tables_dir, exist_ok=True)
    os.makedirs(a.outdir, exist_ok=True)

    rows = []
    metavr_counts = defaultdict(int)
    refseq_counts = defaultdict(int)
    for v in votu_ids:
        m = metavr_best.get(v)
        if m and m["weighted_identity"] >= SPECIES_ANI and m["af_shorter"] >= SPECIES_AF_SHORTER:
            mstatus = "present_metavr_species"
        elif m:
            mstatus = "present_metavr_relative_only"
        else:
            mstatus = "absent_metavr"
        metavr_counts[mstatus] += 1

        rp = refseq_best.get(v)
        rstatus = "refseq_protein_homolog" if rp else "no_refseq_protein_homolog"
        refseq_counts[rstatus] += 1

        meta = info.get(v, {"family": "Unclassified", "host_phylum": "-", "host_species": "-"})
        rows.append({
            "vOTU_ID": v,
            "metavr_status": mstatus,
            "metavr_best_subject": m["sseqid"] if m else "",
            "metavr_weighted_identity": f"{m['weighted_identity']:.2f}" if m else "",
            "metavr_af_shorter": f"{m['af_shorter']:.2f}" if m else "",
            "refseq_status": rstatus,
            "refseq_best_pident": f"{rp['pident']:.1f}" if rp else "",
            "refseq_best_alen_aa": rp["alen"] if rp else "",
            "refseq_best_stitle": rp["stitle"] if rp else "",
            "vmgc_family": meta["family"],
            "vmgc_host_phylum": meta["host_phylum"],
            "vmgc_host_species": meta["host_species"],
            "vmgc_host_genus": host_genus(meta["host_species"]),
            "is_gardnerella": str(is_gardnerella(meta["host_species"])),
        })

    assert len(rows) == 4263, f"vOTU rows {len(rows)} != 4263"
    assert sum(metavr_counts.values()) == 4263

    cols = ["vOTU_ID", "metavr_status", "metavr_best_subject",
            "metavr_weighted_identity", "metavr_af_shorter", "refseq_status",
            "refseq_best_pident", "refseq_best_alen_aa", "refseq_best_stitle",
            "vmgc_family", "vmgc_host_phylum", "vmgc_host_species",
            "vmgc_host_genus", "is_gardnerella"]
    per_votu = os.path.join(a.tables_dir, "table_s24_vmgc_benchmark_overlap.tsv")
    with open(per_votu, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {per_votu} ({len(rows)} vOTUs)")

    # --- contingency summary ---
    # The PRIMARY absence metric is species-level: a vOTU lacks a MetaVR species
    # match unless ANI>=95% over >=85% AF. "no_metavr_species_match" =
    # relative_only + absent (the figure comparable to Huang's species-threshold
    # absence). "absent_metavr" (no >=90% hit at all) is a stricter, deeper-novelty
    # subset reported separately.
    n = len(rows)
    species = metavr_counts["present_metavr_species"]
    relative = metavr_counts["present_metavr_relative_only"]
    absent_metavr = metavr_counts["absent_metavr"]
    no_species = relative + absent_metavr
    no_refseq = refseq_counts["no_refseq_protein_homolog"]
    no_species_and_no_refseq = sum(
        1 for r in rows
        if r["metavr_status"] != "present_metavr_species"
        and r["refseq_status"] == "no_refseq_protein_homolog")
    absent_both = sum(
        1 for r in rows
        if r["metavr_status"] == "absent_metavr"
        and r["refseq_status"] == "no_refseq_protein_homolog")
    summary = os.path.join(a.tables_dir, "table_s24_vmgc_benchmark_overlap_summary.tsv")
    with open(summary, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "n", "pct_of_4263"])
        for k, val in sorted(metavr_counts.items()):
            w.writerow([k, val, f"{100.0*val/n:.1f}"])
        w.writerow(["no_metavr_species_match (relative_only + absent)", no_species,
                    f"{100.0*no_species/n:.1f}"])
        for k, val in sorted(refseq_counts.items()):
            w.writerow([k, val, f"{100.0*val/n:.1f}"])
        w.writerow(["no_metavr_species_AND_no_refseq_protein", no_species_and_no_refseq,
                    f"{100.0*no_species_and_no_refseq/n:.1f}"])
        w.writerow(["absent_metavr(no>=90%hit)_AND_no_refseq_protein", absent_both,
                    f"{100.0*absent_both/n:.1f}"])
    print(f"  wrote {summary}")
    print(f"  MetaVR species match: {species} ({100.0*species/n:.1f}%); "
          f"NO species match: {no_species} ({100.0*no_species/n:.1f}%) "
          f"[of which no detectable >=90% hit at all: {absent_metavr} "
          f"({100.0*absent_metavr/n:.1f}%)]")
    print(f"  no RefSeq protein homolog: {no_refseq} ({100.0*no_refseq/n:.1f}%); "
          f"no MetaVR species AND no RefSeq protein: {no_species_and_no_refseq} "
          f"({100.0*no_species_and_no_refseq/n:.1f}%)")

    # --- stratification by host phylum / family ---
    def stratify(keyfn, label, fname):
        agg = defaultdict(lambda: {"n": 0, "no_species": 0, "absent": 0, "no_refseq": 0})
        for r in rows:
            k = keyfn(r)
            agg[k]["n"] += 1
            if r["metavr_status"] != "present_metavr_species":
                agg[k]["no_species"] += 1
            if r["metavr_status"] == "absent_metavr":
                agg[k]["absent"] += 1
            if r["refseq_status"] == "no_refseq_protein_homolog":
                agg[k]["no_refseq"] += 1
        path = os.path.join(a.outdir, fname)
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow([label, "n_vOTUs", "n_no_metavr_species", "pct_no_metavr_species",
                        "n_absent_no_hit", "pct_absent_no_hit",
                        "n_no_refseq_protein", "pct_no_refseq_protein"])
            for k in sorted(agg, key=lambda x: -agg[x]["n"]):
                d = agg[k]
                w.writerow([k, d["n"], d["no_species"],
                            f"{100.0*d['no_species']/d['n']:.1f}",
                            d["absent"], f"{100.0*d['absent']/d['n']:.1f}",
                            d["no_refseq"], f"{100.0*d['no_refseq']/d['n']:.1f}"])
        print(f"  wrote {path}")

    stratify(lambda r: r["vmgc_host_phylum"], "host_phylum",
             "vmgc_overlap_by_host_phylum.tsv")
    stratify(lambda r: r["vmgc_family"], "viral_family",
             "vmgc_overlap_by_family.tsv")
    stratify(lambda r: ("Gardnerella/B.vaginale" if r["is_gardnerella"] == "True"
                        else r["vmgc_host_genus"]), "host_genus",
             "vmgc_overlap_by_host_genus.tsv")

    # --- consistency cross-check vs manuscript figures ---
    n_classified_family = sum(1 for r in rows if r["vmgc_family"] != "Unclassified")
    n_papilloma = sum(1 for r in rows if r["vmgc_family"] == "Papillomaviridae")
    print("\n  Consistency cross-check (Verification §3b):")
    print(f"    family-classified: {n_classified_family}/4263 "
          f"= {100.0*n_classified_family/n:.1f}% (manuscript: 66%)")
    print(f"    Papillomaviridae vOTUs: {n_papilloma} (manuscript: 61)")
    print("\nPart 2 analysis complete.")


if __name__ == "__main__":
    main()
