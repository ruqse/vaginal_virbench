#!/usr/bin/env python3
"""Render Fig. S15: what the retained multi-evidence sets contain.

Counts are recomputed from the 13 per-sample ``ground_truth_with_kraken2.tsv``
files so the figure stays tied to the benchmark source data.

Run:
  python scripts/07_evaluation/viralm_investigation/06_plot_gt_construction.py
  python .../06_plot_gt_construction.py --with-flow   # legacy 3-panel render

Outputs:
  results/figures/fig_s15.png                              (canonical)
  results/figures/fig_multievidence_gt_construction.{png,pdf}
  results/figures/fig_s15.pdf                              (canonical vector)
  results/figures/fig_s14_source_data.tsv

Scope note (2026-09-10): the hand-authored tier schematic is now Figure S14
(formerly main Figure 2), kept at ``results/figures/fig_s14.png``. The retained-set
composition formerly numbered S14 is now Figure S15. This script writes S15
without overwriting the authored S14. ``draw_flow`` remains reachable via
``--with-flow`` to reproduce the tier counts behind S14. Both figures share
``fig_s14_source_data.tsv``, which records the tier and composition counts.

The negative-grounding labels are disambiguated so that "host" is never used
for two different things: renaming the human category to "Kraken2 human" (not
"human (host)") leaves "host" to refer only to CheckV's cellular "host genes"
field, whose exact wording is kept; the Kraken2 non-human call is written
"Kraken2 bacterial/archaeal".
"""
from __future__ import annotations

import csv
import os
import tempfile
from collections import Counter
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


ROOT = Path(__file__).resolve().parents[3]
GTROOT = ROOT / "results/test_real/ground_truth"
OUTDIR = ROOT / "results/figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

SAMPLES = [
    "UC028_V2",
    "UC055_V1",
    "UC055_V2",
    "UC062_V2",
    "UC065_V2",
    "UC074_V2",
    "UC084_V2",
    "UC093_V2",
    "UC093_V3",
    "UC096_V2",
    "UC115_V2",
    "UC139_V2",
    "UC164_V2",
]
MIN_LEN = 1500

HOMOLOGY_LINES = {"E1a", "E1b", "E2"}
ORTHOGONAL_LINES = {"E4", "E5"}

COLORS = {
    "homology": "#4C72B0",
    "orthogonal": "#DD8452",
    "positive": "#55A868",
    "negative": "#5B7FA6",
    "excluded": "#C44E52",
    "input": "#333333",
    "rule": "#E9E9E9",
    "muted": "#6B6B6B",
    "greybox": "#7F7F7F",
    "bracket": "#9A9A9A",
}


def read_counts() -> dict[str, object]:
    """Recompute all counts used in the figure from per-sample source tables."""
    tier = Counter()
    negative = Counter()
    positive_orthogonal = 0
    positive_homology_only = 0
    tier1 = 0
    tier2 = 0
    n_input = 0

    for sample in SAMPLES:
        path = GTROOT / sample / "ground_truth_with_kraken2.tsv"
        if not path.exists():
            raise FileNotFoundError(f"missing ground-truth source table: {path}")

        with path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if int(row["length"]) < MIN_LEN:
                    continue

                n_input += 1
                tier_id = int(row["tier"])
                tier[tier_id] += 1

                raw_evidence = row.get("evidence_lines", "-") or "-"
                evidence = set(raw_evidence.split(",")) if raw_evidence != "-" else set()
                notes = row.get("notes", "") or ""

                if tier_id in (1, 2):
                    tier1 += tier_id == 1
                    tier2 += tier_id == 2

                    if evidence & ORTHOGONAL_LINES:
                        positive_orthogonal += 1
                    elif evidence and evidence.issubset(HOMOLOGY_LINES):
                        positive_homology_only += 1
                elif tier_id == 0:
                    # Internal bucket keys are stable identifiers; display labels
                    # (draw_composition_bars) are disambiguated separately.
                    if "Homo sapiens" in notes or "9606" in notes:
                        negative["Kraken2 human (host)"] += 1
                    elif notes.startswith("kraken2:"):
                        negative["Kraken2 bacterial"] += 1
                    elif notes.startswith("checkv"):
                        negative["CheckV host genes, 0 viral genes"] += 1

    n_positive = tier1 + tier2
    n_negative = tier[0]
    n_excluded = tier[3] + tier[-1]
    n_benchmark = n_positive + n_negative

    return {
        "tier": tier,
        "negative": negative,
        "positive_orthogonal": positive_orthogonal,
        "positive_homology_only": positive_homology_only,
        "positive_other": n_positive - positive_orthogonal - positive_homology_only,
        "tier1": tier1,
        "tier2": tier2,
        "n_input": n_input,
        "n_positive": n_positive,
        "n_negative": n_negative,
        "n_excluded": n_excluded,
        "n_benchmark": n_benchmark,
        "ratio": n_negative / max(n_positive, 1),
    }


def rounded_box(
    ax: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    text: str,
    facecolor: str,
    *,
    edgecolor: str | None = None,
    textcolor: str = "white",
    fontsize: float = 9.0,
    weight: str = "bold",
    alpha: float = 1.0,
    hatch: str | None = None,
    rounding: float = 1.5,
) -> FancyBboxPatch:
    patch = FancyBboxPatch(
        (x - width / 2, y - height / 2),
        width,
        height,
        boxstyle=f"round,pad=0.38,rounding_size={rounding}",
        linewidth=1.2,
        facecolor=facecolor,
        edgecolor=edgecolor or facecolor,
        alpha=alpha,
        hatch=hatch,
        zorder=2,
    )
    ax.add_patch(patch)
    ax.text(
        x,
        y,
        text,
        ha="center",
        va="center",
        color=textcolor,
        fontsize=fontsize,
        fontweight=weight,
        linespacing=1.15,
        zorder=3,
    )
    return patch


def arrow(
    ax: plt.Axes,
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    *,
    color: str = "#555555",
    linewidth: float = 1.6,
    radius: float = 0.0,
    mutation_scale: float = 13.0,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            (x1, y1),
            (x2, y2),
            arrowstyle="-|>",
            mutation_scale=mutation_scale,
            linewidth=linewidth,
            color=color,
            connectionstyle=f"arc3,rad={radius}",
            zorder=1,
        )
    )


def bracket(
    ax: plt.Axes,
    x0: float,
    x1: float,
    y: float,
    label: str,
    *,
    color: str = "#9A9A9A",
    fontsize: float = 10.5,
    drop: float = 0.9,
) -> None:
    """Draw a labelled bracket spanning [x0, x1] with end-ticks dropping toward
    the boxes below and the label centred above."""
    ax.plot(
        [x0, x0, x1, x1],
        [y - drop, y, y, y - drop],
        color=color,
        linewidth=1.2,
        zorder=1,
    )
    ax.text(
        (x0 + x1) / 2,
        y + 1.1,
        label,
        ha="center",
        va="bottom",
        fontsize=fontsize,
        color=color,
        fontweight="bold",
    )


def pct(count: int, denom: int) -> float:
    return 100 * count / denom if denom else 0.0


def draw_flow(ax: plt.Axes, counts: dict[str, object]) -> None:
    tier: Counter = counts["tier"]  # type: ignore[assignment]

    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    # ---- title -----------------------------------------------------------
    ax.text(
        50,
        99.2,
        "Construction of the multi-evidence ground truth",
        ha="center",
        va="top",
        fontsize=18,
        fontweight="bold",
    )
    ax.text(
        50,
        95.4,
        "13 UChoose vaginal shotgun metagenomes; contig counts at the >=1,500 bp benchmark window",
        ha="center",
        va="top",
        fontsize=11,
        color=COLORS["muted"],
    )

    # ================= Panel A -- Evidence inputs =========================
    ax.text(2, 90.6, "A", ha="left", va="center", fontsize=15, fontweight="bold")
    ax.text(6.5, 90.6, "Evidence inputs", ha="left", va="center", fontsize=13, fontweight="bold")
    ax.text(
        6.5,
        87.6,
        "Viral evidence groups are counted per contig; non-viral (human/bacterial) evidence grounds negatives.",
        ha="left",
        va="center",
        fontsize=9.5,
        color=COLORS["muted"],
    )

    # brackets over the two evidence classes
    bracket(ax, 3.0, 72.5, 84.4, "viral evidence groups", color=COLORS["bracket"], fontsize=10.5)
    bracket(ax, 78.0, 97.0, 84.4, "negative grounding", color=COLORS["bracket"], fontsize=10.5)

    # viral evidence boxes E1-E5 (colour still encodes reference-dependent vs
    # reference-free); a single non-viral grounding box sits apart on the right.
    viral_boxes = [
        (10.0, "E1 homology\nDIAMOND + MetaVR", COLORS["homology"], 13.6),
        (24.5, "E2 CheckV\nviral genes", COLORS["homology"], 12.2),
        (38.0, "E3 Phigaro\nprophage", COLORS["homology"], 12.2),
        (52.0, "E4 CRISPR\nhost spacer", COLORS["orthogonal"], 12.2),
        (66.0, "E5 RCA phi29\nvirome enrichment", COLORS["orthogonal"], 13.2),
    ]
    for x, label, color, width in viral_boxes:
        rounded_box(ax, x, 78.5, width, 6.8, label, color, fontsize=8.6, rounding=1.25)

    rounded_box(
        ax,
        87.5,
        78.5,
        18.6,
        7.6,
        "Non-viral evidence\nKraken2 human or bacterial/archaeal,\nor CheckV host genes",
        COLORS["greybox"],
        fontsize=7.6,
        rounding=1.25,
    )

    # ================= Panel B -- Decision funnel =========================
    ax.text(2, 70.6, "B", ha="left", va="center", fontsize=15, fontweight="bold")
    ax.text(6.5, 70.6, "Decision funnel", ha="left", va="center", fontsize=13, fontweight="bold")

    rounded_box(
        ax,
        50,
        65.0,
        44,
        5.6,
        f"All assembled contigs >= 1,500 bp    n = {counts['n_input']:,}",
        COLORS["input"],
        fontsize=12.0,
    )
    arrow(ax, 50, 62.0, 50, 59.9, color=COLORS["input"], linewidth=1.5)

    rounded_box(
        ax,
        50,
        56.7,
        66,
        5.4,
        "Count independent viral evidence groups per contig   {E1, E2, E3, E4, E5}",
        COLORS["rule"],
        textcolor="#222222",
        fontsize=11,
    )

    tier_specs = [
        (10.0, ">= 3 groups\nTier 1\nn = {0:,}".format(tier[1]), COLORS["positive"], 1.0, None),
        (28.0, "2 groups\nTier 2\nn = {0:,}".format(tier[2]), COLORS["positive"], 1.0, None),
        (48.0, "0 groups +\nnon-viral evidence\nTier 0\nn = {0:,}".format(tier[0]), COLORS["negative"], 1.0, None),
        (68.0, "1 group\nTier 3\nexcluded\nn = {0:,}".format(tier[3]), COLORS["excluded"], 0.58, "////"),
        (88.0, "0 groups +\nno evidence\nTier -1\nn = {0:,}".format(tier[-1]), COLORS["excluded"], 0.58, "////"),
    ]
    for x, label, color, alpha, hatch in tier_specs:
        rounded_box(ax, x, 46.0, 16.0, 8.6, label, color, fontsize=8.5, alpha=alpha, hatch=hatch)
        arrow(ax, 50, 53.9, x, 50.5, color="#8C8C8C", linewidth=1.1,
              radius=(0.10 if x < 50 else -0.10) if abs(x - 50) > 2 else 0.0)

    # tier -> outcome
    arrow(ax, 10.0, 41.6, 18.0, 33.4, color=COLORS["positive"], linewidth=1.7, radius=-0.07)
    arrow(ax, 28.0, 41.6, 20.5, 33.4, color=COLORS["positive"], linewidth=1.7, radius=0.07)
    arrow(ax, 48.0, 41.6, 48.0, 33.4, color=COLORS["negative"], linewidth=1.7)
    arrow(ax, 68.0, 41.6, 79.0, 33.4, color=COLORS["excluded"], linewidth=1.7, radius=-0.07)
    arrow(ax, 88.0, 41.6, 81.0, 33.4, color=COLORS["excluded"], linewidth=1.7, radius=0.07)

    rounded_box(
        ax,
        19.0,
        29.8,
        30.0,
        6.6,
        f"Retain as viral positives\nTier 1 + Tier 2    n = {counts['n_positive']:,}",
        COLORS["positive"],
        fontsize=10.2,
    )
    rounded_box(
        ax,
        48.5,
        29.8,
        26.0,
        6.6,
        f"Retain as grounded negatives\nTier 0    n = {counts['n_negative']:,}",
        COLORS["negative"],
        fontsize=10.2,
    )
    rounded_box(
        ax,
        80.0,
        29.8,
        30.0,
        6.6,
        "Exclude ambiguous / unsupported\nn = {0:,} ({1:.0f}%)".format(
            counts["n_excluded"], pct(int(counts["n_excluded"]), int(counts["n_input"]))
        ),
        COLORS["excluded"],
        fontsize=10.2,
        alpha=0.66,
        hatch="////",
    )

    # outcome -> benchmark
    arrow(ax, 19.0, 26.5, 33.0, 19.4, color=COLORS["positive"], linewidth=1.8, radius=0.08)
    arrow(ax, 48.5, 26.5, 40.0, 19.4, color=COLORS["negative"], linewidth=1.8, radius=-0.08)
    rounded_box(
        ax,
        35.0,
        15.6,
        48.0,
        6.2,
        "Benchmark set\n{0:,} contigs    viral : non-viral = 1 : {1:.0f}".format(
            counts["n_benchmark"], counts["ratio"]
        ),
        COLORS["input"],
        fontsize=12.0,
    )


def draw_stacked_bar(
    ax: plt.Axes,
    segments: list[tuple[str, int, str]],
    total: int,
    *,
    title: str,
    small_label: str | None = None,
) -> None:
    left = 0
    ax.set_xlim(0, total)
    ax.set_ylim(-0.55, 0.55)
    ax.axis("off")
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold", pad=6)

    for label, value, color in segments:
        ax.barh(0, value, left=left, height=0.48, color=color, edgecolor="white", linewidth=1.1)
        if value / total >= 0.065:
            ax.text(
                left + value / 2,
                0,
                f"{value:,}\n{pct(value, total):.0f}%",
                ha="center",
                va="center",
                color="white",
                fontsize=9.2,
                fontweight="bold",
                linespacing=1.0,
            )
        left += value

    if small_label:
        ax.text(total * 0.99, 0.33, small_label, ha="right", va="center", fontsize=8.5, color="#555555")

    legend_handles = [Rectangle((0, 0), 1, 1, color=color) for label, value, color in segments]
    ax.legend(
        legend_handles,
        [label for label, value, color in segments],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=len(segments),
        fontsize=8.1,
        frameon=False,
        handlelength=1.1,
        columnspacing=1.25,
    )


def draw_composition_bars(
    fig: plt.Figure,
    counts: dict[str, object],
    *,
    standalone: bool = False,
) -> None:
    """Draw the two retained-set composition bars.

    ``standalone=True`` is the Fig. S15 layout (the bars are the whole figure,
    tagged A and B). ``standalone=False`` is the legacy ``--with-flow`` layout,
    where they sit under the decision funnel as a single panel C.
    """
    if standalone:
        fig.text(0.035, 0.90, "A", ha="left", va="center", fontsize=15, fontweight="bold")
        fig.text(0.535, 0.90, "B", ha="left", va="center", fontsize=15, fontweight="bold")
        ax_left = fig.add_axes([0.06, 0.30, 0.40, 0.30])
        ax_right = fig.add_axes([0.56, 0.30, 0.40, 0.30])
    else:
        fig.text(0.045, 0.235, "C", ha="left", va="center", fontsize=15, fontweight="bold")
        fig.text(0.075, 0.235, "What the retained sets contain", ha="left", va="center",
                 fontsize=13, fontweight="bold")
        ax_left = fig.add_axes([0.06, 0.075, 0.40, 0.13])
        ax_right = fig.add_axes([0.56, 0.075, 0.40, 0.13])

    n_positive = int(counts["n_positive"])
    positive_segments = [
        ("orthogonal corroboration (E4/E5)", int(counts["positive_orthogonal"]), COLORS["orthogonal"]),
        ("homology-only (E1/E2)", int(counts["positive_homology_only"]), COLORS["homology"]),
        ("other evidence combos", int(counts["positive_other"]), "#B7B7B7"),
    ]
    draw_stacked_bar(
        ax_left,
        positive_segments,
        n_positive,
        title=f"Positive set composition  (n = {n_positive:,})",
        small_label=f"other: {counts['positive_other']:,} ({pct(int(counts['positive_other']), n_positive):.0f}%)",
    )

    negative: Counter = counts["negative"]  # type: ignore[assignment]
    n_negative = int(counts["n_negative"])
    # Display labels disambiguate "host": renaming the human category to
    # "Kraken2 human" leaves "host" to mean only CheckV's cellular host-genes
    # field, so CheckV's exact term "host genes" is kept as-is.
    negative_segments = [
        ("Kraken2 human", negative.get("Kraken2 human (host)", 0), "#8C6BB1"),
        ("Kraken2 bacterial/archaeal", negative.get("Kraken2 bacterial", 0), "#3690C0"),
        ("CheckV host genes, 0 viral", negative.get("CheckV host genes, 0 viral genes", 0), COLORS["negative"]),
    ]
    draw_stacked_bar(
        ax_right,
        negative_segments,
        n_negative,
        title=f"Negative set grounding  (n = {n_negative:,})",
    )


def write_source_data(counts: dict[str, object]) -> None:
    tier: Counter = counts["tier"]  # type: ignore[assignment]
    negative: Counter = counts["negative"]  # type: ignore[assignment]
    n_input = int(counts["n_input"])
    n_positive = int(counts["n_positive"])
    n_negative = int(counts["n_negative"])

    rows = [
        ("input", "assembled_contigs_ge_1500bp", n_input, ""),
        ("tier", "tier_1_ge_3_groups", tier[1], pct(tier[1], n_input)),
        ("tier", "tier_2_eq_2_groups", tier[2], pct(tier[2], n_input)),
        ("tier", "tier_0_nonviral_evidence", tier[0], pct(tier[0], n_input)),
        ("tier", "tier_3_eq_1_group_excluded", tier[3], pct(tier[3], n_input)),
        ("tier", "tier_minus_1_ambiguous_excluded", tier[-1], pct(tier[-1], n_input)),
        ("outcome", "positive_tier_1_plus_2", n_positive, pct(n_positive, n_input)),
        ("outcome", "negative_tier_0", n_negative, pct(n_negative, n_input)),
        ("outcome", "excluded_tier_3_plus_minus_1", counts["n_excluded"], pct(int(counts["n_excluded"]), n_input)),
        ("outcome", "benchmark_positive_plus_negative", counts["n_benchmark"], pct(int(counts["n_benchmark"]), n_input)),
        (
            "positive_composition",
            "orthogonal_corroboration_E4_or_E5",
            counts["positive_orthogonal"],
            pct(int(counts["positive_orthogonal"]), n_positive),
        ),
        (
            "positive_composition",
            "homology_only_E1_or_E2",
            counts["positive_homology_only"],
            pct(int(counts["positive_homology_only"]), n_positive),
        ),
        (
            "positive_composition",
            "other_evidence_combos",
            counts["positive_other"],
            pct(int(counts["positive_other"]), n_positive),
        ),
        (
            "negative_grounding",
            "Kraken2_human",
            negative.get("Kraken2 human (host)", 0),
            pct(negative.get("Kraken2 human (host)", 0), n_negative),
        ),
        (
            "negative_grounding",
            "Kraken2_bacterial_archaeal",
            negative.get("Kraken2 bacterial", 0),
            pct(negative.get("Kraken2 bacterial", 0), n_negative),
        ),
        (
            "negative_grounding",
            "CheckV_host_genes_0_viral",
            negative.get("CheckV host genes, 0 viral genes", 0),
            pct(negative.get("CheckV host genes, 0 viral genes", 0), n_negative),
        ),
    ]

    with (OUTDIR / "fig_s14_source_data.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["section", "item", "count", "percent_within_section"])
        for section, item, count, percent in rows:
            writer.writerow(
                [
                    section,
                    item,
                    int(count),
                    "" if percent == "" else f"{float(percent):.2f}",
                ]
            )


def build_figure(counts: dict[str, object], *, with_flow: bool = False) -> plt.Figure:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
            "font.size": 10,
            "axes.titleweight": "bold",
            "savefig.facecolor": "white",
        }
    )

    if with_flow:
        # Legacy 3-panel render: evidence inputs + decision funnel + composition.
        # The evidence inputs and tier funnel are shown in Figure S14; kept here to regenerate and
        # verify the counts drawn in that hand-authored schematic.
        fig = plt.figure(figsize=(16, 10), facecolor="white")
        flow_ax = fig.add_axes([0.03, 0.275, 0.94, 0.70])
        draw_flow(flow_ax, counts)
        draw_composition_bars(fig, counts)
        return fig

    fig = plt.figure(figsize=(14, 4.2), facecolor="white")
    draw_composition_bars(fig, counts, standalone=True)
    return fig


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--with-flow",
        action="store_true",
        help="Render the legacy 3-panel figure (evidence inputs + decision funnel "
             "+ composition) instead of the composition-only Fig. S15. Use this to "
             "regenerate the counts behind the hand-made Figure S14.",
    )
    args = parser.parse_args()

    counts = read_counts()
    write_source_data(counts)

    fig = build_figure(counts, with_flow=args.with_flow)
    if args.with_flow:
        # Never write canonical figure names from the legacy render — S15 is the
        # composition-only figure, and the flow now belongs to Figure S14.
        outputs = [
            OUTDIR / "fig_multievidence_gt_construction.png",
            OUTDIR / "fig_multievidence_gt_construction.pdf",
        ]
    else:
        outputs = [
            OUTDIR / "fig_s15.png",                          # canonical raster
            OUTDIR / "fig_s15.pdf",                          # canonical vector (submission)
        ]
    for path in outputs:
        fig.savefig(path, dpi=400)
    plt.close(fig)

    print(
        "counts:",
        {
            "tier": dict(counts["tier"]),  # type: ignore[arg-type]
            "positive_orthogonal": counts["positive_orthogonal"],
            "positive_homology_only": counts["positive_homology_only"],
            "negative": dict(counts["negative"]),  # type: ignore[arg-type]
        },
    )
    for path in outputs:
        print("wrote", path)
    print("wrote", OUTDIR / "fig_s14_source_data.tsv")


if __name__ == "__main__":
    main()
