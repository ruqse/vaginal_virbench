#!/usr/bin/env python3
"""
14_export_gene_map_data.py — source data for the dark-matter gene-map figure.

Writes to results/vmgc_crosscheck/dark_matter_metavr_refseq/:
  gene_map_genes.tsv      one row per geNomad-predicted gene on the 48 contigs,
                          with Pfam 35.0 families (hmmscan --cut_ga) and category
  gene_map_tool_calls.tsv contig x tool positive call (0/1) at default
                          thresholds, current Track C parsers
  gene_map_topology.tsv   geNomad topology for contigs in its plasmid summary
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

PROJ = Path(__file__).resolve().parents[2]
O = PROJ / "results/vmgc_crosscheck/dark_matter_metavr_refseq"
G = PROJ / "results/test_real/full_run/RCA_master/genomad"
sys.path.insert(0, str(PROJ / "scripts/03_assembly"))
import evaluate_rca_benchmark as ev  # noqa: E402

REP_MOB = {"Rep_1", "Rep_2", "Rep_3", "RepL", "Rep_trans", "Mob_Pre"}

dark = [l.strip() for l in open(PROJ / "results/test_real/coassembly/rca_dark_matter_contig_ids.txt") if l.strip()]
dset = set(dark)


def contig_of(gene):
    return gene.rsplit("_", 1)[0]


# Pfam families per gene (best-scoring first)
fams = defaultdict(list)
for line in open(O / "dark_matter_pfam35.domtbl"):
    if line.startswith("#"):
        continue
    f = line.split()
    fam, gene, score = f[0], f[3], float(f[13])
    fams[gene].append((score, fam))

with open(O / "gene_map_genes.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["contig_id", "gene", "start", "end", "strand", "pfam_families", "category"])
    for r in csv.DictReader(open(G / "RCA_master_contigs_annotate/RCA_master_contigs_genes.tsv"), delimiter="\t"):
        c = contig_of(r["gene"])
        if c not in dset:
            continue
        fl = [f for _, f in sorted(set(fams.get(r["gene"], [])), reverse=True)]
        fl = list(dict.fromkeys(fl))
        if set(fl) & REP_MOB:
            cat = "Plasmid replication / mobilisation"
        elif any(f == "Zot" or f.startswith("Phage_") for f in fl):
            cat = "Phage-type domain"
        elif fl:
            cat = "Other Pfam domain"
        else:
            cat = "No Pfam domain"
        w.writerow([c, r["gene"], r["start"], r["end"], r["strand"], ";".join(fl), cat])

scores = ev.discover_tool_outputs(PROJ / "results/test_real/full_run/RCA_master", "RCA_master", gt_ids=dset)
with open(O / "gene_map_tool_calls.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["contig_id", "tool", "call"])
    for tool, sc in scores.items():
        s = ev.fold_scores(sc, tool, dset)
        thr = ev.TOOL_THRESHOLDS.get(tool, 0.5)
        pthr = ev.TOOL_PVALUE_THRESHOLDS.get(tool)
        for c in dark:
            v = s.get(c)
            if v is None:
                call = 0
            elif isinstance(v, tuple):
                call = int(v[0] >= thr and (pthr is None or v[1] < pthr))
            else:
                call = int(float(v) >= thr)
            w.writerow([c, tool, call])
    for tool in set(ev.TOOL_THRESHOLDS) | {"VirSorter", "Sourmash", "VIBRANT", "MetaPhinder"}:
        if tool not in scores:
            for c in dark:
                w.writerow([c, tool, 0])

with open(O / "gene_map_topology.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["contig_id", "topology"])
    for r in csv.DictReader(open(G / "RCA_master_contigs_summary/RCA_master_contigs_plasmid_summary.tsv"), delimiter="\t"):
        if r["seq_name"] in dset:
            w.writerow([r["seq_name"], r["topology"]])
print("done")
