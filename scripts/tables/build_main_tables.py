#!/usr/bin/env python3
"""Build (and --check) the two compact main tables from their cited sources.

  Main Table 1 -> table1_cross_track_ranking.tsv    (curated 5 tools x 4 tracks, MCC)
  Main Table 2 -> table2_tool_recommendations.tsv   (compact: 8 actionable tools)

(The former data-source summary main table was removed -- redundant with Fig. 1,
 which now carries the data sources and per-track sample/ground-truth counts.
 Recover from git history if a consolidated summary is ever needed again.)

Numeric counts/metrics are pulled (and asserted) from NAMED source files so the
main tables stay reproducible -- no hand edits. Usage:
  python build_main_tables.py            # (re)write the two TSVs
  python build_main_tables.py --check    # rebuild in memory, diff vs disk, exit 1 on drift
"""
import sys
import csv
import re
from decimal import Decimal, InvalidOperation, ROUND_HALF_UP
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]         # repo root (scripts/tables/..)
HERE = ROOT / "results" / "tables"                 # published tables
RES = ROOT / "results"
sys.path.insert(0, str(ROOT / "scripts/07_evaluation"))
from metric_display import format_mcc_from_counts

CURATED = ["geNomad", "ViraLM", "VirSorter2", "VIBRANT", "Jaeger"]


def _rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _mcc_map(path, tool_col, mcc_col):
    return {r[tool_col]: r[mcc_col] for r in _rows(path)}


def _fmt(x):
    try:
        return str(Decimal(str(x)).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))
    except (TypeError, ValueError, InvalidOperation):
        return str(x)


# ---------- Main Table 1: cross-track ranking (curated 5 tools) ----------
# Confidence intervals come from the published per-track supplementary tables.
# Points are recomputed from confusion counts and rounded once for display;
# legacy table2_/table3_ files are no longer build inputs.
def build_cross_track_ranking():
    # Read the published per-track point estimates AND intervals. The former
    # builder silently dropped the CI columns already present in Table 1.
    ta = {r["tool"]: r for r in _rows(HERE / "table_s4_bootstrap_confidence_intervals.tsv")
          if r["length_bp"] == "1500" and r["metric"] == "MCC"}
    ta_counts = {r["tool"]: r for r in _rows(HERE / "table_s34_track_a_overall_metrics_full.tsv")}
    tc = {r["tool"]: r for r in _rows(HERE / "table_s35_track_c_overall_metrics_full.tsv")}
    me = {r["tool"]: r for r in _rows(HERE / "table_s12_secondary_benchmark_overall.tsv")}
    mi = {r["tool"]: r for r in _rows(HERE / "table_s27_mitch_ranking_ci.tsv")}
    mi_counts = {r["tool"]: r for r in _rows(RES / "test_real/viralm_investigation/mitch_scaleup30/overall_metrics.tsv")}
    header = ["Tool", "Track A MCC (1,500 bp)", "Track A MCC (1,500 bp) 95% CI",
              "Track C MCC", "Track C MCC 95% CI", "Multi-evidence MCC",
              "Multi-evidence MCC 95% CI", "MiTCH MCC", "MiTCH MCC 95% CI"]
    def interval(row, lo, hi):
        return f"{_fmt(row[lo])} to {_fmt(row[hi])}"
    rows = []
    for t in CURATED:
        # Jaeger cannot run below its 2,048 bp window, so Track A at 1,500 bp is undefined.
        ta_cell = "N/A (min 2,048 bp)" if t == "Jaeger" else format_mcc_from_counts(ta_counts[t])
        ta_ci = "N/A" if t == "Jaeger" else interval(ta[t], "ci_lo", "ci_hi")
        rows.append([t, ta_cell, ta_ci,
                     format_mcc_from_counts(tc[t]), interval(tc[t], "MCC_CI_low", "MCC_CI_high"),
                     format_mcc_from_counts(me[t]), interval(me[t], "MCC_CI_low", "MCC_CI_high"),
                     format_mcc_from_counts(mi_counts[t]), interval(mi[t], "CI_low", "CI_high")])
    return [header] + rows


# ---------- Main Table 2: compact recommendations (8 actionable tools) ----------
# Compact, actionable view: Tier | Tool | Best for | Key evidence | Requirements.
# Per-tool MCCs are in Table 1; the full 12-column matrix for all 14 tools (both
# MCCs, 3-way eukaryotic recall, runtimes, strengths/caveats) lives in
# table2_tool_recommendations_full.tsv -> supplementary workbook sheet Table2_full.
# Key-evidence phrases are faithful one-line condensations of the "Key Strength"
# column of that full table (this benchmark's own results, not external claims).
KEY_EVIDENCE = {
    "geNomad": "Highest Track A MCC at 1,500 bp; primary-assembly MCC difference from VIBRANT unresolved",
    "ViraLM": "100% recall in all three eukaryote-infecting virus categories at 1,500 bp; 65/65 homology-free fragments",
    "Jaeger": "Recovered 41/48 candidate dark-matter contigs; fastest where it runs (>= 2,048 bp)",
    "VirSorter2": "Highest Track C MCC (0.60); Track A MCC 0.501 at 3,000 bp",
    "PPR-Meta": "Ranked 4th on Track C; three-class phage/chromosome/plasmid model",
    "TransGINmer": "81% herpesvirus recall; beats feature-based ML on short contigs",
    "VIBRANT": "Highest MiTCH MCC; high precision at lower recall on primary assemblies",
    "DeepVirFinder": "58% dark-matter detection; high Lactobacillus-phage recall",
}


# Peak resident memory is a deployment constraint on a par with VRAM, but the
# Requirements cell used to declare it only for GPU tools -- so geNomad, the
# tool with by far the largest footprint (~18.5 GB, 310x the smallest in the
# panel), read as a bare "CPU". Surface RSS for CPU tools once it is large
# enough to dictate node choice.
RSS_DECLARE_MB = 2048


def _peak_rss_mb(tool):
    """Track A 1,500 bp peak RSS for one tool, from Table S16 (None if absent)."""
    path = HERE / "table_s16_resource_usage_L1500.tsv"
    with open(path, newline="") as fh:
        body = (ln for ln in fh if not ln.startswith("#"))
        for r in csv.DictReader(body, delimiter="\t"):
            if r["tool"] == tool and r["condition"] == "TrackA_L1500":
                return float(r["peak_rss_mb"])
    return None


def _requirements(gpu, min_len, tool=None):
    """Collapse 'GPU Required' + 'Min Contig Length' into one compact cell."""
    if gpu.strip().startswith("No"):
        base = "CPU"
        rss = _peak_rss_mb(tool) if tool else None
        if rss is not None and rss >= RSS_DECLARE_MB:
            base = f"CPU ({rss / 1024:.1f} GB RAM)"
    else:
        vram = gpu[gpu.find("(") + 1: gpu.find(")")].replace("VRAM", "").strip()
        base = f"GPU ({vram})"
    ml = (min_len or "").strip()
    if ml and ml not in ("None", "-"):
        # thousands separator on the numeric part: "2048 bp" -> "2,048 bp"
        ml = re.sub(r"\d{4,}", lambda m: f"{int(m.group()):,}", ml)
        return f"{base}; min {ml}"
    return base


def build_tool_recommendations():
    full = _rows(HERE / "table2_tool_recommendations_full.tsv")
    keep = ("Recommended", "Complementary")
    header = ["Tier", "Tool", "Best for", "Key evidence", "Requirements"]
    rows = []
    prophage = {r["tool"]: r for r in _rows(HERE / "table_s13_q2_prophage_detection.tsv")}
    for r in sorted(full, key=lambda row: {"Recommended": 0, "Complementary": 1}.get(row["Tier"], 2)):
        if r["Tier"] not in keep:
            continue
        tool = r["Tool"]
        assert tool in KEY_EVIDENCE, f"missing key-evidence phrase for {tool!r}"
        evidence = KEY_EVIDENCE[tool]
        if tool == "VIBRANT":
            p = prophage[tool]
            recall_pct = 100 * int(p["prophage_tp"]) / int(p["prophage_total"])
            evidence += f"; {recall_pct:.1f}% prophage recall"
        rows.append([r["Tier"], tool, r["Research Goal"],
                     evidence,
                     _requirements(r["GPU Required"], r["Min Contig Length"], tool)])
    assert len(rows) == 8, f"expected 8 actionable rows, got {len(rows)}"
    return [header] + rows


OUTPUTS = {
    "table1_cross_track_ranking.tsv": build_cross_track_ranking,
    "table2_tool_recommendations.tsv": build_tool_recommendations,
}


def to_tsv(table):
    return "".join("\t".join(row) + "\n" for row in table)


def main():
    check = "--check" in sys.argv[1:]
    drift = []
    for name, fn in OUTPUTS.items():
        text = to_tsv(fn())
        path = HERE / name
        if check:
            cur = path.read_text(encoding="utf-8") if path.exists() else ""
            status = "ok   " if cur == text else "DRIFT"
            if cur != text:
                drift.append(name)
            print(f"{status} {name}")
        else:
            path.write_text(text, encoding="utf-8")
            print(f"WROTE {name}  ({len(fn()) - 1} data rows)")
    if check and drift:
        print(f"\nFAIL: {len(drift)} table(s) differ from sources: {', '.join(drift)}")
        sys.exit(1)
    if check:
        print("PASS: both main tables reproduce from their cited sources.")


if __name__ == "__main__":
    main()
