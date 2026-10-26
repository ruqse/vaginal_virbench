#!/usr/bin/env python3
"""Refresh Tables S14/S15 with the current VIBRANT parser and publish Tables S42/S43.

S14/S15 were last written on 2026-03-19, before viral_output_parsers.py (2026-09-21)
began counting only VIBRANT's final accepted phages. Re-scoring the 48 contigs with the
current parser found VIBRANT at 2/48, not 4/48, and every other tool unchanged. This script takes the current calls from the
reassessment outputs, asserts that only VIBRANT moves, and rewrites both the benchmark
sources (which publish_parser_revision_tables.py copies) and the results/tables copies.

Usage: python3 scripts/10_vmgc_crosscheck/16_refresh_dark_matter_tables.py
"""
import csv
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NEW = ROOT / "results/vmgc_crosscheck/dark_matter_metavr_refseq"
BENCH = ROOT / "results/test_real/coassembly/benchmark"
TAB = ROOT / "results/tables"
S14 = ("dark_matter_per_contig.tsv", "table_s14_dark_matter_per_contig.tsv")
S15 = ("dark_matter_detection_rates.tsv", "table_s15_dark_matter_detection_rates.tsv")
DROPPED = {"UC062_V2__NODE_857_length_5246_cov_15.061327",
           "UC065_V2__NODE_337_length_3506_cov_24.741032"}


def read(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(rows, *paths):
    for path in paths:
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\r\n")  # repo tables are CRLF
            w.writeheader()
            w.writerows(rows)


# S15: current per-tool counts are the "DM" column, formatted "41/48 (85%)".
dm = {r["tool"]: int(r["DM"].split("/")[0]) for r in read(NEW / "tool_calls_by_plasmid_group.tsv")}
s15 = read(BENCH / S15[0])
assert {r["tool"] for r in s15} == set(dm), "tool sets differ"
for r in s15:
    old = int(r["n_detected"])
    if r["tool"] == "VIBRANT":
        assert old in (4, 2) and dm["VIBRANT"] == 2, (old, dm["VIBRANT"])  # 2 on a rerun
    else:
        assert old == dm[r["tool"]], (r["tool"], old, dm[r["tool"]])
    r["n_detected"] = str(dm[r["tool"]])
    r["detection_rate"] = str(round(dm[r["tool"]] / int(r["n_dark_matter"]), 4))

# Per-contig current calls (long format, CRLF). dark_matter_integrated_evidence.tsv is
# not used for tools: Pelle's script 13 copied detecting_tools from the stale S14.
cur = {}
for r in read(NEW / "gene_map_tool_calls.tsv"):
    if r["call"].strip() == "1":
        cur.setdefault(r["contig_id"], set()).add(r["tool"])


def fix_tools(r):
    """Drop VIBRANT from the two contigs it no longer calls; every other row must match."""
    old = set(r["detecting_tools"].split(";")) - {""}
    want = old - {"VIBRANT"} if r["contig_id"] in DROPPED else old
    assert want == cur.get(r["contig_id"], set()), (r["contig_id"], old ^ cur.get(r["contig_id"], set()))
    r["detecting_tools"] = ";".join(t for t in r["detecting_tools"].split(";") if t in want)
    r["n_tools_detecting"] = str(len(want))


s14 = read(BENCH / S14[0])
s42 = read(NEW / "dark_matter_integrated_evidence.tsv")
for table in (s14, s42):
    assert len(table) == 48 and {r["contig_id"] for r in table} >= set(cur)
    for r in table:
        fix_tools(r)
    assert sum("VIBRANT" in r["detecting_tools"] for r in table) == dm["VIBRANT"] == 2

write(s15, BENCH / S15[0], TAB / S15[1])
write(s14, BENCH / S14[0], TAB / S14[1])
write(s42, TAB / "table_s42_unclassified_enriched_evidence.tsv")
shutil.copyfile(NEW / "tool_calls_by_plasmid_group.tsv", TAB / "table_s43_tool_calls_by_plasmid_group.tsv")
print("S14/S15 refreshed (VIBRANT 4 -> 2 of 48); S42/S43 written to", TAB)
