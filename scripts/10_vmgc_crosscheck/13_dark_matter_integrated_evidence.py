#!/usr/bin/env python3
"""
13_dark_matter_integrated_evidence.py — one row per dark-matter contig combining
plasmid-database, virus-database, protein-domain and classifier evidence.

PLSDB / MetaVR (BLASTn, best of megablast and dc-megablast, interval-merged
coverage estimator from vmgc_coverage.py):
  *_near_identical : >= 95 % identity over >= 85 % of the contig (same rule as
                     the MetaVR species-level call in 10_*.py)
  *_partial        : any hit at e <= 1e-5 below that rule
  *_absent         : no hit
Pfam 35.0 (hmmscan --cut_ga on geNomad-predicted proteins):
  plasmid_rep_mob  : Rep_1, Rep_2, Rep_3, RepL, Rep_trans or Mob_Pre domain
  phage_domain     : Zot or any Phage_* domain
geNomad: plasmid call = listed in the plasmid summary (geNomad's own filters).
"""
import csv
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from vmgc_coverage import best_subject_per_query  # noqa: E402

PROJ = Path(__file__).resolve().parents[2]
O = PROJ / "results/vmgc_crosscheck/dark_matter_metavr_refseq"
GEN = PROJ / "results/test_real/full_run/RCA_master/genomad"
REP_MOB = {"Rep_1", "Rep_2", "Rep_3", "RepL", "Rep_trans", "Mob_Pre"}
ID, COV = 95.0, 85.0


def best_of(paths, lens):
    hits = [best_subject_per_query(p) for p in paths]
    out = {}
    for cid, L in lens.items():
        cands = [h[cid] for h in hits if cid in h]
        if not cands:
            continue
        near = [h for h in cands
                if h["weighted_identity"] >= ID and 100 * h["covered_q"] / L >= COV]
        out[cid] = max(near or cands, key=lambda h: (h["covered_q"], h["total_bitscore"]))
        out[cid]["near"] = bool(near)
    return out


def titles(paths):
    t = {}
    for p in paths:
        for line in open(p):
            f = line.rstrip("\n").split("\t")
            if len(f) > 14:
                t[f[1]] = f[14]
    return t


def main():
    per = list(csv.DictReader(open(O / "dark_matter_metavr_refseq_per_contig.tsv"), delimiter="\t"))
    lens = {r["contig_id"]: int(r["length_bp"]) for r in per}
    pl_paths = [O / "dark_vs_plsdb_megablast.tsv", O / "dark_vs_plsdb_dc-megablast.tsv"]
    pl = best_of(pl_paths, lens)
    pt = titles(pl_paths)

    pf = defaultdict(set)
    for line in open(O / "dark_matter_pfam35_family_by_contig.tsv"):
        fam, cid = line.split()
        pf[cid].add(fam)

    plas_call = {r["seq_name"] for r in csv.DictReader(
        open(GEN / "RCA_master_contigs_summary/RCA_master_contigs_plasmid_summary.tsv"), delimiter="\t")}
    agg = {r["seq_name"]: r for r in csv.DictReader(
        open(GEN / "RCA_master_contigs_aggregated_classification/RCA_master_contigs_aggregated_classification.tsv"),
        delimiter="\t")}
    s14 = {r["contig_id"]: r for r in csv.DictReader(
        open(PROJ / "docs/tables/table_s14_dark_matter_per_contig.tsv"), delimiter="\t")}

    rows, tally = [], Counter()
    for r in per:
        cid, L = r["contig_id"], lens[r["contig_id"]]
        h = pl.get(cid)
        pstat = "plsdb_near_identical" if h and h["near"] else "plsdb_partial" if h else "plsdb_absent"
        fams = pf.get(cid, set())
        rep = bool(fams & REP_MOB)
        phage = sorted(f for f in fams if f == "Zot" or f.startswith("Phage_"))
        tally[pstat] += 1
        tally["pfam_rep_mob"] += rep
        tally["pfam_phage_domain"] += bool(phage)
        tally["genomad_plasmid_call"] += cid in plas_call
        rows.append({
            "contig_id": cid, "length_bp": L,
            "enrichment_R": s14.get(cid, {}).get("enrichment_R", ""),
            "plsdb_status": pstat,
            "plsdb_best_subject": h["sseqid"] if h else "",
            "plsdb_best_title": pt.get(h["sseqid"], "") if h else "",
            "plsdb_identity": f"{h['weighted_identity']:.2f}" if h else "",
            "plsdb_contig_cov_pct": f"{100 * h['covered_q'] / L:.1f}" if h else "",
            "metavr_status": r["metavr_status"],
            "metavr_contig_cov_pct": r["metavr_contig_cov_pct"],
            "refseq_status": r["refseq_status"],
            "refseq_title": r["refseq_title"],
            "pfam_families": ";".join(sorted(fams)),
            "pfam_plasmid_rep_mob": int(rep),
            "pfam_phage_domain": ";".join(phage),
            "genomad_plasmid_score": agg.get(cid, {}).get("plasmid_score", ""),
            "genomad_virus_score": agg.get(cid, {}).get("virus_score", ""),
            "genomad_plasmid_call": int(cid in plas_call),
            "n_tools_detecting": s14.get(cid, {}).get("n_tools_detecting", ""),
            "detecting_tools": s14.get(cid, {}).get("detecting_tools", ""),
        })

    out = O / "dark_matter_integrated_evidence.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(dict(tally))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
