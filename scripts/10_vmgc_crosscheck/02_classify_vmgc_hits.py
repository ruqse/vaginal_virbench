#!/usr/bin/env python3
"""
02_classify_vmgc_hits.py — Part 1 classifier (VMGC cross-check).

Labels the benchmark's own novel sequences against the FULL VMGC catalogue,
with DELIBERATELY DIFFERENT treatment by query type:

  * 15-genome ANI-novel panel (full genomes -> species-level claims valid):
      vmgc_species_match : ANI >= 95 % over >= 85 % AF of the shorter sequence
                           (MIUViG/Roux 2019). skani preferred; BLASTn fallback.
      vmgc_relative      : detectable homology below the species threshold
      vmgc_absent        : no BLASTn hit and no skani line

  * 48 RCA dark-matter contigs (short fragments -> homology claims ONLY):
      vmgc_homolog              : VMGC member match >= 95 % identity over
                                  >= 50 % of the CONTIG ("VMGC would recover it")
      vmgc_partial              : detectable but weaker hit
      vmgc_no_detectable_homology : dark even to the full vaginal catalogue
    NB: no species/vOTU-absence claim is made for these short contigs.

Uses the corrected interval-merged coverage parser in vmgc_coverage.py.
"""

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vmgc_coverage import best_subject_per_query, parse_skani  # noqa: E402

SPECIES_ANI = 95.0
SPECIES_AF_SHORTER = 85.0      # full-genome species criterion
FRAG_IDENT = 95.0
FRAG_QCOV = 50.0               # fragment-level "VMGC would recover it"


def read_fasta_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def load_member_map(path):
    m = {}
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            m[row["member_id"]] = row["vOTU_ID"]
    return m


def load_vmgc_info(path):
    """vOTU_ID -> (family, host_phylum, host_species)."""
    info = {}
    with open(path, newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        fam_c = "Taxonomic_classification_(family_level)"
        php_c = "Host_assignment_(phylum_level)"
        hsp_c = "Host_assignment_(species_level)"
        for row in r:
            info[row["vOTU_ID"]] = (
                row.get(fam_c, ""), row.get(php_c, ""), row.get(hsp_c, ""),
            )
    return info


def load_id_mapping(path):
    """short_id -> dict (incl. metavr_best_pident for homology_free flag)."""
    m = {}
    if not path or not os.path.isfile(path):
        return m
    with open(path, newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            m[row["short_id"]] = row
    return m


def resolve_votu(member_id, member_map, info):
    votu = member_map.get(member_id, "")
    fam, php, hsp = info.get(votu, ("", "", ""))
    return votu, fam, php, hsp


def classify_dark(blast_path, fasta_path, member_map, info, outdir):
    best = best_subject_per_query(blast_path)
    ids = read_fasta_ids(fasta_path)
    rows, counts = [], {"vmgc_homolog": 0, "vmgc_partial": 0,
                        "vmgc_no_detectable_homology": 0}
    for cid in ids:
        b = best.get(cid)
        if not b:
            status = "vmgc_no_detectable_homology"
            votu = fam = php = hsp = ""
            wi = qaf = saf = 0.0
            qlen = slen = 0
            sid = ""
        else:
            wi = b["weighted_identity"]
            qaf = b["q_af"]
            saf = b["s_af"]
            qlen, slen = b["qlen"], b["slen"]
            sid = b["sseqid"]
            votu, fam, php, hsp = resolve_votu(sid, member_map, info)
            if wi >= FRAG_IDENT and qaf >= FRAG_QCOV:
                status = "vmgc_homolog"
            else:
                status = "vmgc_partial"
        counts[status] += 1
        rows.append({
            "contig_id": cid, "vmgc_status": status,
            "best_member": sid, "best_vOTU": votu,
            "weighted_identity": f"{wi:.2f}", "contig_cov_pct": f"{qaf:.2f}",
            "subject_cov_pct": f"{saf:.2f}", "contig_len": qlen,
            "vmgc_family": fam, "vmgc_host_phylum": php, "vmgc_host_species": hsp,
        })
    _write(os.path.join(outdir, "dark_matter_vmgc_classification.tsv"), rows,
           ["contig_id", "vmgc_status", "best_member", "best_vOTU",
            "weighted_identity", "contig_cov_pct", "subject_cov_pct",
            "contig_len", "vmgc_family", "vmgc_host_phylum", "vmgc_host_species"])
    _write_counts(os.path.join(outdir, "dark_matter_vmgc_counts.tsv"), counts)
    assert len(rows) == 48, f"dark rows {len(rows)} != 48"
    return counts


def classify_panel(blast_path, fasta_path, skani_path, id_map,
                   member_map, info, outdir):
    best = best_subject_per_query(blast_path)
    skani = {}
    if skani_path and os.path.isfile(skani_path) and os.path.getsize(skani_path) > 0:
        try:
            skani = parse_skani(skani_path)
        except Exception as e:  # noqa: BLE001
            print(f"  [warn] could not parse skani ({e}); BLASTn-only species calls")
    ids = read_fasta_ids(fasta_path)
    rows, counts = [], {"vmgc_species_match": 0, "vmgc_relative": 0,
                        "vmgc_absent": 0}
    for gid in ids:
        b = best.get(gid)
        s = skani.get(gid)
        homology_free = id_map.get(gid, {}).get("metavr_best_pident", "") == "no_hit"

        # species evidence: skani preferred, BLASTn fallback
        skani_species = bool(s) and s["ani"] >= SPECIES_ANI and s["af_shorter"] >= SPECIES_AF_SHORTER
        blast_species = bool(b) and b["weighted_identity"] >= SPECIES_ANI and b["af_shorter"] >= SPECIES_AF_SHORTER

        if skani_species or blast_species:
            status = "vmgc_species_match"
        elif b or s:
            status = "vmgc_relative"
        else:
            status = "vmgc_absent"
        counts[status] += 1

        # matched vOTU: prefer BLASTn member (member_map keyed on those IDs)
        if b:
            votu, fam, php, hsp = resolve_votu(b["sseqid"], member_map, info)
            match_id = b["sseqid"]
        elif s:
            votu, fam, php, hsp = resolve_votu(s["ref"], member_map, info)
            match_id = s["ref"]
        else:
            votu = fam = php = hsp = match_id = ""

        rows.append({
            "genome_id": gid, "homology_free": str(homology_free),
            "vmgc_status": status, "best_match_member": match_id, "best_vOTU": votu,
            "skani_ani": f"{s['ani']:.2f}" if s else "",
            "skani_af_shorter": f"{s['af_shorter']:.2f}" if s else "",
            "blast_weighted_identity": f"{b['weighted_identity']:.2f}" if b else "",
            "blast_af_shorter": f"{b['af_shorter']:.2f}" if b else "",
            "metavr_best_pident": id_map.get(gid, {}).get("metavr_best_pident", ""),
            "vmgc_family": fam, "vmgc_host_phylum": php, "vmgc_host_species": hsp,
        })
    _write(os.path.join(outdir, "novel_panel_vmgc_classification.tsv"), rows,
           ["genome_id", "homology_free", "vmgc_status", "best_match_member",
            "best_vOTU", "skani_ani", "skani_af_shorter",
            "blast_weighted_identity", "blast_af_shorter", "metavr_best_pident",
            "vmgc_family", "vmgc_host_phylum", "vmgc_host_species"])
    _write_counts(os.path.join(outdir, "novel_panel_vmgc_counts.tsv"), counts)
    assert len(rows) == 15, f"panel rows {len(rows)} != 15"
    return counts


def _write(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {path} ({len(rows)} rows)")


def _write_counts(path, counts):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["vmgc_status", "n"])
        for k, v in counts.items():
            w.writerow([k, v])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dark-blast", required=True)
    ap.add_argument("--dark-fasta", required=True)
    ap.add_argument("--panel-blast", required=True)
    ap.add_argument("--panel-fasta", required=True)
    ap.add_argument("--panel-skani", default="")
    ap.add_argument("--id-mapping", default="")
    ap.add_argument("--member-map", required=True)
    ap.add_argument("--vmgc-info", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    member_map = load_member_map(a.member_map)
    info = load_vmgc_info(a.vmgc_info)
    id_map = load_id_mapping(a.id_mapping)

    print("Part 1 — dark-matter contigs (fragment-level):")
    dc = classify_dark(a.dark_blast, a.dark_fasta, member_map, info, a.outdir)
    print(f"  {dc}")
    print("Part 1 — ANI-novel panel (species-level):")
    pc = classify_panel(a.panel_blast, a.panel_fasta, a.panel_skani, id_map,
                        member_map, info, a.outdir)
    print(f"  {pc}")
    print("Part 1 classification complete.")


if __name__ == "__main__":
    main()
