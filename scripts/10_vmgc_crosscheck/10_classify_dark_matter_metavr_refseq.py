#!/usr/bin/env python3
"""
10_classify_dark_matter_metavr_refseq.py — per-contig MetaVR / RefSeq status
for the 48 Track C dark-matter contigs (output of 09_*.sh).

MetaVR (nucleotide; best subject per query from the interval-merged estimator
in vmgc_coverage.py, best of megablast and dc-megablast):
  metavr_species   : weighted ANI >= 95 % over >= 85 % of the shorter sequence
                     (Roux et al., 2019 vOTU rule; same as Tables S24/S25)
  metavr_relative  : any hit at e <= 1e-5 below the species rule
  metavr_absent    : no hit at e <= 1e-5 under either task

RefSeq viral protein (DIAMOND blastx --very-sensitive):
  refseq_homolog   : >= 30 % identity over >= 50 aa at e <= 1e-10 (Methods rule)
  refseq_weak      : a hit at e <= 1e-5 failing the homolog rule
  refseq_absent    : no hit at e <= 1e-5
"""

import argparse
import csv
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vmgc_coverage import best_subject_per_query  # noqa: E402

SPECIES_ANI, SPECIES_AF = 95.0, 85.0
RS_ID, RS_ALEN, RS_E = 30.0, 50, 1e-10


def fasta_lengths(path):
    lens, cur = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].split()[0]
                lens[cur] = 0
            else:
                lens[cur] += len(line.strip())
    return lens


def is_species(h):
    return h and h["weighted_identity"] >= SPECIES_ANI and h["af_shorter"] >= SPECIES_AF


def pick_metavr(mb, dc):
    """Prefer a species-level hit; else the hit covering more of the contig."""
    cands = [h for h in (mb, dc) if h]
    if not cands:
        return None, ""
    sp = [h for h in cands if is_species(h)]
    pool = sp or cands
    best = max(pool, key=lambda h: (h["covered_q"], h["total_bitscore"]))
    return best, ("megablast" if best is mb else "dc-megablast")


def refseq_hits(path):
    best_any, best_hom = {}, {}
    if not os.path.isfile(path):
        return best_any, best_hom
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 14:
                continue
            h = {"sseqid": f[1], "pident": float(f[2]), "alen": int(f[3]),
                 "evalue": float(f[10]), "bitscore": float(f[11]),
                 "stitle": f[14] if len(f) > 14 else ""}
            q = f[0]
            if q not in best_any or h["bitscore"] > best_any[q]["bitscore"]:
                best_any[q] = h
            if h["pident"] >= RS_ID and h["alen"] >= RS_ALEN and h["evalue"] <= RS_E:
                if q not in best_hom or h["bitscore"] > best_hom[q]["bitscore"]:
                    best_hom[q] = h
    return best_any, best_hom


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dark-fasta", required=True)
    ap.add_argument("--megablast", required=True)
    ap.add_argument("--dcmegablast", required=True)
    ap.add_argument("--refseq", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    lens = fasta_lengths(a.dark_fasta)
    mb = best_subject_per_query(a.megablast)
    dc = best_subject_per_query(a.dcmegablast)
    rs_any, rs_hom = refseq_hits(a.refseq)

    rows, mc, rc, joint = [], Counter(), Counter(), Counter()
    for cid, L in lens.items():
        m, task = pick_metavr(mb.get(cid), dc.get(cid))
        mstat = ("metavr_species" if is_species(m) else
                 "metavr_relative" if m else "metavr_absent")
        r = rs_hom.get(cid) or rs_any.get(cid)
        rstat = ("refseq_homolog" if cid in rs_hom else
                 "refseq_weak" if r else "refseq_absent")
        mc[mstat] += 1
        rc[rstat] += 1
        joint[(mstat, rstat)] += 1
        rows.append({
            "contig_id": cid, "length_bp": L,
            "metavr_status": mstat, "metavr_task": task,
            "metavr_best_subject": m["sseqid"] if m else "",
            "metavr_weighted_identity": f"{m['weighted_identity']:.2f}" if m else "",
            "metavr_contig_cov_pct": f"{100 * m['covered_q'] / L:.1f}" if m else "",
            "metavr_af_shorter_pct": f"{m['af_shorter']:.1f}" if m else "",
            "metavr_best_evalue": f"{m['best_evalue']:.2g}" if m else "",
            "refseq_status": rstat,
            "refseq_best_subject": r["sseqid"] if r else "",
            "refseq_pident": f"{r['pident']:.1f}" if r else "",
            "refseq_aln_aa": r["alen"] if r else "",
            "refseq_evalue": f"{r['evalue']:.2g}" if r else "",
            "refseq_title": r["stitle"] if r else "",
        })

    os.makedirs(a.outdir, exist_ok=True)
    out = os.path.join(a.outdir, "dark_matter_metavr_refseq_per_contig.tsv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    summ = os.path.join(a.outdir, "dark_matter_metavr_refseq_summary.tsv")
    with open(summ, "w") as fh:
        fh.write("metavr_status\trefseq_status\tn_contigs\n")
        for (ms, rs), n in sorted(joint.items()):
            fh.write(f"{ms}\t{rs}\t{n}\n")
    print(f"n contigs: {len(rows)}")
    print("MetaVR:", dict(mc))
    print("RefSeq:", dict(rc))
    print("Joint:", {f"{k[0]}|{k[1]}": v for k, v in joint.items()})
    print(f"wrote {out}\nwrote {summ}")


if __name__ == "__main__":
    main()
