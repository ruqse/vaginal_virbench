#!/usr/bin/env python3
"""Validate the curated four-per-arm Track B backgrounds, or print their runs.

PRJNA1170175 is retained for Track A, but its pooled enriched viromes are not
eligible bacterial-community backgrounds for the primary Track B CST contrast.
Saved outputs from the former ten-background analysis are historical inputs;
consumers must use this manifest, never discover backgrounds from result dirs.
"""
from __future__ import annotations

import argparse
import csv
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "data/expansion_cohorts/track_b_backgrounds.tsv"
EXCLUDED_COHORTS = frozenset({"PRJNA1170175"})
EXPECTED_TOOLS = (
    "DeepVirFinder", "VirSorter", "VirSorter2", "VirFinder", "PPR-Meta",
    "VIBRANT", "Seeker", "MetaPhinder", "Sourmash", "geNomad", "HVSeeker",
    "Jaeger", "TransGINmer", "ViraLM",
)
# Completion artifacts for the current eight-background Track B input set.
# Empty/header-only results can represent a valid no-hit result. Existence,
# together with exit_code=0, distinguishes those from missing/failed outputs.
OUTPUT_PATTERNS = {
    "DeepVirFinder": ("deepvirfinder/*_dvfpred.txt",),
    "VirSorter": ("virsorter/work/VIRSorter_global-phage-signal.csv",),
    "VirSorter2": ("virsorter2/final-viral-score.tsv",),
    "VirFinder": ("virfinder/virfinder_results.tsv", "virfinder/results.txt"),
    "PPR-Meta": ("pprmeta/*.csv",),
    "VIBRANT": ("vibrant/**/VIBRANT_machine_*.tsv", "vibrant/**/VIBRANT_phages_*.fasta"),
    "Seeker": ("seeker/seeker_results.tsv",),
    "MetaPhinder": ("metaphinder/output.txt",),
    "Sourmash": ("sourmash/sourmash_results.csv",),
    "geNomad": ("genomad/**/*_virus_summary.tsv",),
    "HVSeeker": ("hvseeker/hvseeker_results.tsv",),
    "Jaeger": ("jaeger/**/*_default_jaeger.tsv", "jaeger/jaeger_results.tsv"),
    "TransGINmer": ("transginmer/transginmer_results.tsv",),
    "ViraLM": ("viralm/viralm_results.tsv",),
}
EXPECTED_ARMS = {
    "SRR27287961": "CST-I",
    "SRR27287984": "CST-I",
    "SRR35936354": "CST-I",
    "SRR27287985": "CST-I",
    "SRR27287964": "CST-IV-B",
    "SRR27287966": "CST-IV-B",
    "SRR35936363": "CST-IV-B",
    "SRR34514383": "CST-IV-B",
}


def validate_backgrounds(rows):
    if any(row["cohort"] in EXCLUDED_COHORTS for row in rows):
        raise ValueError("PRJNA1170175 enriched viromes are excluded from Track B")
    run_arms = {row["run"]: row["stratum"] for row in rows}
    if len(rows) != 8 or run_arms != EXPECTED_ARMS:
        raise ValueError("Track B requires the curated eight backgrounds with their original CST arms")
    if Counter(row["stratum"] for row in rows) != {"CST-I": 4, "CST-IV-B": 4}:
        raise ValueError("Track B requires four CST-I and four CST-IV-B backgrounds")
    return rows


def load_backgrounds(path=MANIFEST):
    with Path(path).open(newline="") as handle:
        return validate_backgrounds(list(csv.DictReader(handle, delimiter="\t")))


def restore_metric_precision(result):
    """Derive unrounded metrics from the shared evaluator's exact confusion counts.

The shared evaluator rounds its reported metrics to four decimals. Aggregates
and bootstrap endpoints should instead use the underlying integer counts.
"""
    tp, fp, tn, fn = (result[name] for name in ("TP", "FP", "TN", "FN"))
    denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return dict(result,
                MCC=(tp * tn - fp * fn) / denominator if denominator else 0.0,
                precision=tp / (tp + fp) if tp + fp else 0.0,
                recall=tp / (tp + fn) if tp + fn else 0.0)


def load_execution_records(backgrounds, root=ROOT):
    """Load complete terminal records for exactly the manifest-selected runs."""
    records = []
    for background in backgrounds:
        run = background["run"]
        run_records = []
        run_dir = Path(root) / "results/expansion/track_b/benchmark" / f"trackB_{run}_cov10"
        for name in ("resource_usage.tsv", "resource_usage_gpu.tsv"):
            with (run_dir / name).open(newline="") as handle:
                run_records.extend(dict(row, background=run) for row in csv.DictReader(handle, delimiter="\t"))
        if len(run_records) != 14 or {row["tool"] for row in run_records} != set(EXPECTED_TOOLS):
            raise ValueError(f"{run}: incomplete execution records for 14 tools")
        records.extend(run_records)
    return records


def validate_successful_executions(records, root=ROOT):
    """Reject failed executions and missing artifacts before publication scoring.

The CPU/GPU wrappers record skipped=yes, exit_code=0 both for cached outputs
and for unavailable dependencies. A skipped record is accepted only when its
expected result artifact exists. No positive viral call is required.
"""
    failures = [f"{r['background']} {r['tool']} exit_code={r.get('exit_code', 'missing')}"
                for r in records if r.get("exit_code") != "0"]
    if failures:
        raise ValueError("Unsuccessful Track B executions: " + "; ".join(failures))
    for row in records:
        run, tool = row["background"], row["tool"]
        run_dir = Path(root) / "results/expansion/track_b/benchmark" / f"trackB_{run}_cov10"
        if not any(path.is_file() for pattern in OUTPUT_PATTERNS[tool]
                   for path in run_dir.glob(pattern)):
            raise ValueError(f"{run} {tool}: missing completed result artifact "
                             f"(skipped={row.get('skipped', 'unknown')})")
    return records


def require_successful_executions(backgrounds, root=ROOT):
    return validate_successful_executions(load_execution_records(backgrounds, root), root)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--check-completion", action="store_true",
                        help="Require successful tool records and output artifacts for the eight runs")
    args = parser.parse_args()
    backgrounds = load_backgrounds(args.manifest)
    if args.check_completion:
        records = require_successful_executions(backgrounds)
        print(f"All {len(records)} selected run/tool records and output artifacts passed.")
    else:
        for background in backgrounds:
            print(background["run"])
