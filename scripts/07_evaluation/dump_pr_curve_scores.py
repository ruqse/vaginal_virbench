#!/usr/bin/env python3
"""
dump_pr_curve_scores.py — export the per-fragment score matrix behind AUPRC.

Why this exists
---------------
`benchmark_fragments.py` writes a scalar AUPRC. This script re-parses the raw
tool outputs with the same parsers and exports full-precision scores for
Figure 2B. The reported flag distinguishes observed scores from imputed zeros,
so curves can end at the last available score without extrapolation. Dots are
loaded separately from the published metrics and include any decision filters.

It imports `benchmark_fragments` rather than reimplementing any parsing, so the
scores here are identical by construction to those used for the published
AUPRC values.

Usage
-----
  python3 scripts/07_evaluation/dump_pr_curve_scores.py \
      --manifest data/spike_in/fragments/fragment_manifest.tsv \
      --results-root results/spike_in_benchmark \
      --lengths L1500 \
      --output results/spike_in_benchmark/pr_curve_scores.tsv.gz
"""

import argparse
import gzip
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import benchmark_fragments as bf  # noqa: E402
from deterministic_gzip import DeterministicGzipWriter  # noqa: E402


def dump_length(results_root: Path, length_bin: str, truth: dict, out_rows: list) -> None:
    """Append (length, tool, fragment, label, score, pvalue) rows for one length bin."""
    results_dir = results_root / f"{length_bin}_fragments"
    if not results_dir.is_dir():
        print(f"  [{length_bin}] SKIP - no directory {results_dir}")
        return

    target_len = int(length_bin[1:])
    sub_truth = {fid: info for fid, info in truth.items()
                 if info["fragment_length"] == target_len}
    if not sub_truth:
        print(f"  [{length_bin}] SKIP - manifest has no fragments of this length")
        return

    for tool, parser in bf.TOOL_PARSERS.items():
        try:
            if tool in {"VirSorter", "VirSorter2", "VIBRANT", "geNomad"}:
                scores = parser(results_dir, manifest_ids=set(sub_truth))
            else:
                scores = parser(results_dir)
        except Exception as exc:                      # noqa: BLE001
            print(f"  [{length_bin}] {tool:14} parser error: {exc}")
            continue
        if not scores:
            # e.g. Jaeger at L500-L1500: requires fragments >= 2,048 bp, so it
            # produced no output at all. Emit nothing; the figure marks it n/a.
            print(f"  [{length_bin}] {tool:14} no output")
            continue

        n_written = 0
        # Match evaluate_tool()'s deterministic input order. Full float precision
        # below preserves distinct raw scores instead of creating extra ties.
        for fid in sorted(sub_truth):
            info = sub_truth[fid]
            val = scores.get(fid)
            if val is None:
                # Tool returned no row for this fragment -> scored as negative,
                # matching evaluate_tool()'s treatment of missing predictions.
                # `reported = 0` marks it as IMPUTED, not measured: tools that
                # pre-filter their output (geNomad --min-score, VirSorter2,
                # VIBRANT) leave thousands of fragments unranked, and a PR curve
                # drawn through them shows tie-order noise rather than tool
                # behaviour. Keeping the flag lets the figure plot the available
                # scores and omit unsupported continuations.
                score, pval, reported = 0.0, "", "0"
            elif isinstance(val, tuple):
                score, pval, reported = float(val[0]), val[1], "1"
            elif isinstance(val, bool):
                score, pval, reported = (1.0 if val else 0.0), "", "1"
            else:
                score, pval, reported = float(val), "", "1"
            out_rows.append((
                length_bin, tool, fid,
                "1" if info["label"] == "viral" else "0",
                f"{score:.17g}", str(pval), reported,
            ))
            n_written += 1
        print(f"  [{length_bin}] {tool:14} {n_written} fragments")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--results-root", required=True)
    ap.add_argument("--lengths", nargs="*", default=["L1500"])
    ap.add_argument("--output", required=True)
    ap.add_argument("--metrics-table", type=Path,
                    default=HERE.parent.parent / "results/tables/table_s7_track_a_overall_metrics_by_length.tsv",
                    help="Current published metrics used for the AUPRC self-check (default: Table S7)")
    args = ap.parse_args()

    truth = bf.load_manifest(args.manifest)
    print(f"Manifest: {len(truth)} fragments")

    rows: list = []
    for lb in args.lengths:
        dump_length(Path(args.results_root), lb, truth, rows)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    opener = DeterministicGzipWriter if out.suffix == ".gz" else (lambda path: open(path, "wt"))
    with opener(out) as fh:
        fh.write("length_bin\ttool\tfragment_id\tlabel\tscore\tpvalue\treported\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"\nWrote {len(rows)} rows -> {out}")

    # Self-check: recompute AUPRC from the dumped scores and compare against the
    # current Table S7. Historical per-run benchmark_results.tsv files can
    # contain an older AUPRC implementation and must not define figure values.
    import csv
    import numpy as np
    with args.metrics_table.open() as fh:
        published = {(r["length_bin"], r["tool"]): r["AUPRC"]
                     for r in csv.DictReader(fh, delimiter="\t")}
    failures = []
    print("\nAUPRC self-check (dumped vs published):")
    for lb in args.lengths:
        for tool in bf.TOOL_PARSERS:
            sub = [(int(r[3]), float(r[4])) for r in rows if r[0] == lb and r[1] == tool]
            if not sub:
                continue
            y = np.array([s[0] for s in sub])
            s = np.array([s[1] for s in sub])
            if len(np.unique(s)) <= 1:
                continue
            got = round(bf._auprc(y, s), 4)
            want = published.get((lb, tool), "NA")
            flag = "ok " if want not in ("NA", "") and abs(got - float(want)) < 5e-5 else "DIFF"
            if flag == "DIFF":
                failures.append(f"{lb}/{tool}")
            print(f"  {flag} {lb} {tool:14} dumped={got:.4f} published={want}")
    if failures:
        raise SystemExit("AUPRC mismatch: " + ", ".join(failures))


if __name__ == "__main__":
    main()
