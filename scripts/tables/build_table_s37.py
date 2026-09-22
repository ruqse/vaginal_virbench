#!/usr/bin/env python3
"""Build Table S37 -- per-tool recall on the ANI-novel panel, by novelty stratum.

Why this table exists
---------------------
The Q3 result in the manuscript (Results, ANI-novel panel; Figure 4) rested on
numbers that appeared in NO numbered supplementary table: the per-tool
homology-free recall series and the 65/65, 350/427 and 282/427 counts. They were
correct and reproducible, but only from
`data/novel_spike_discovery/full_scale/dark_vs_relatives_stratification.tsv`,
which is not part of the supplement. A reviewer could not check the paper's most
novel claim. Asset audit item B1, 2026-08-11.

S37 appends after S36, so S1-S36 keep their numbers and nothing renumbers.

Two columns are added to the source file, both derived and both documented in the
Table S37 caption:

  stratum_label       'homology-free' / 'homology-detectable' -- the manuscript's
                      wording for the source's 'dark' / 'has_relatives'.
  homology_dependence the axis the Q3 sentence actually claims on. It is NOT the
                      five-approach method taxonomy: VirFinder sits in
                      "Feature-based ML" yet scores pure k-mer composition. The
                      assignment follows the verified_algorithm column of
                      tool_method_taxonomy_verified.tsv, tool by tool.

    python3 scripts/tables/build_table_s37.py
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]         # repo root (scripts/tables/..)
HERE = ROOT / "results" / "tables"                 # published tables
SRC = ROOT / "data/novel_spike_discovery/full_scale/dark_vs_relatives_stratification.tsv"
TAX = HERE / "tool_method_taxonomy_verified.tsv"
OUT = HERE / "table_s37_novel_panel_recall_by_stratum.tsv"

# Consults a reference, marker set, or HMM profile at classification time.
DEPENDENT = {"MetaPhinder", "Sourmash", "VirSorter", "VirSorter2", "VIBRANT", "geNomad"}
# Scores sequence composition alone (k-mer, CNN, LSTM, transformer).
COMPOSITION = {"VirFinder", "DeepVirFinder", "PPR-Meta", "Jaeger",
               "HVSeeker", "Seeker", "ViraLM", "TransGINmer"}

STRATUM = {"dark": "homology-free", "has_relatives": "homology-detectable"}
LENGTH = {"L1500": 1500, "L3000": 3000}

HEADER = ["tool", "method_approach", "homology_dependence", "stratum", "stratum_label",
          "fragment_length_bp", "n_genomes", "n_fragments", "TP", "FN", "recall",
          "recall_delta_vs_homology_detectable", "evaluable"]


def main() -> int:
    rows = list(csv.DictReader(SRC.read_text(encoding="utf-8").splitlines(), delimiter="\t"))
    fam = {r["tool"]: r["method_approach"]
           for r in csv.DictReader(TAX.read_text(encoding="utf-8").splitlines(), delimiter="\t")}

    tools = {r["tool"] for r in rows}
    unassigned = tools - DEPENDENT - COMPOSITION
    if unassigned:
        raise SystemExit(f"tools with no homology_dependence assignment: {sorted(unassigned)}")
    missing_fam = tools - set(fam)
    if missing_fam:
        raise SystemExit(f"tools absent from {TAX.name}: {sorted(missing_fam)}")

    out = []
    for r in rows:
        t = r["tool"]
        length = LENGTH[r["length"]]
        # Jaeger returns no classification below its 2,048 bp window, so its
        # 1,500 bp zeros are a no-op rather than a recall of zero. Flagging this
        # in the table is what keeps the Q3 comparison honest -- without it the
        # 74% floor reads as false.
        evaluable = "no" if (t == "Jaeger" and length < 2048) else "yes"
        out.append([
            t, fam[t],
            "marker/reference-dependent" if t in DEPENDENT else "composition-based",
            r["stratum"], STRATUM[r["stratum"]], length,
            r["n_genomes"], r["n_fragments"], r["TP"], r["FN"], r["recall"],
            r["recall_delta_vs_relatives"], evaluable,
        ])

    # Deterministic order: dependence class, then approach, tool, length, stratum.
    out.sort(key=lambda x: (x[2] != "marker/reference-dependent", x[1], x[0], x[5], x[3]))

    with OUT.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(HEADER)
        w.writerows(out)
    print(f"WROTE {OUT.relative_to(ROOT)}  ({len(out)} rows, {len(tools)} tools)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
