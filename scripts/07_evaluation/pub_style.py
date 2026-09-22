#!/usr/bin/env python3
"""
Shared publication style configuration for all benchmark figures.

Sets matplotlib rcParams for journal-ready figures (300 DPI, Arial/Helvetica,
consistent font sizes, clean spines).

Usage:
    import pub_style
    pub_style.apply()
    # ... then create figures as normal
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

# ---------------------------------------------------------------------------
# Publication rcParams
# ---------------------------------------------------------------------------

def apply():
    """Apply publication-ready matplotlib rcParams globally."""
    plt.rcParams.update({
        # Font
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 10,
        # Axes
        "axes.labelsize": 12,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        # Ticks
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        # Legend
        "legend.fontsize": 9,
        "legend.title_fontsize": 10,
        "legend.frameon": False,
        "legend.borderpad": 0.3,
        # Figure
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.1,
        # Lines
        "lines.linewidth": 1.5,
        "lines.markersize": 6,
        # Patches
        "patch.linewidth": 0.5,
    })


# ---------------------------------------------------------------------------
# Tool metadata (single source of truth)
# ---------------------------------------------------------------------------

TOOL_FAMILY = {
    "DeepVirFinder": "Deep learning",
    "Seeker":        "Deep learning",
    "PPR-Meta":      "Deep learning",
    "HVSeeker":      "Deep learning",
    "Jaeger":        "Deep learning",
    "TransGINmer":   "Attention/transformer",
    "ViraLM":        "Attention/transformer",
    "VirSorter":     "Feature-based ML",
    "VirSorter2":    "Feature-based ML",
    "VIBRANT":       "Feature-based ML",
    "VirFinder":     "Feature-based ML",
    "geNomad":       "DNN + markers",
    "MetaPhinder":   "Reference-based",
    "Sourmash":      "Reference-based",
}

FAMILY_PALETTE = {
    "Deep learning":          "#1f77b4",
    "Attention/transformer":  "#9467bd",
    "Feature-based ML":       "#2ca02c",
    "DNN + markers":          "#ff7f0e",
    "Reference-based":        "#d62728",
}

TOOL_ORDER = [
    "geNomad", "ViraLM", "TransGINmer", "VirSorter2", "Jaeger",
    "HVSeeker", "PPR-Meta", "VIBRANT", "DeepVirFinder", "Seeker",
    "VirSorter", "MetaPhinder", "VirFinder", "Sourmash",
]

# Distinct per-tool colors (colorblind-friendlier palette)
_tab20 = plt.cm.tab20(np.linspace(0, 1, 20))
TOOL_COLORS = {t: mcolors.to_hex(_tab20[i]) for i, t in enumerate(TOOL_ORDER)}

TOOL_SCOPE = {
    "DeepVirFinder": "all_virus",
    "VirSorter":     "phage_only",
    "VirSorter2":    "all_virus",
    "VirFinder":     "phage_only",
    "PPR-Meta":      "phage_and_plasmid",
    "VIBRANT":       "all_virus",
    "Seeker":        "phage_only",
    "MetaPhinder":   "phage_only",
    "Sourmash":      "db_dependent",
    "geNomad":       "all_virus",
    "HVSeeker":      "all_virus",
    "Jaeger":        "phage_only",
    "TransGINmer":   "all_virus",
    "ViraLM":        "all_virus",
}

SCOPE_LABEL_COLORS = {
    "phage_only":        "#e74c3c",
    "all_virus":         "#2980b9",
    "phage_and_plasmid": "#8e44ad",
    "db_dependent":      "#7f8c8d",
}

# Length bins
ALL_LENGTHS = ["L500", "L1000", "L1500", "L3000", "L5000", "L10000"]

LENGTH_COLORS = {
    "L500":   "#c6dbef",
    "L1000":  "#9ecae1",
    "L1500":  "#6baed6",
    "L3000":  "#3182bd",
    "L5000":  "#08519c",
    "L10000": "#08306b",
}

# Category metadata
PHAGE_CATS = ["lactobacillus_phage", "gardnerella_phage",
              "megasphaera_phage", "fannyhessea_phage", "sneathia_phage"]
EUKARYOTIC_CATS = ["hpv", "herpesvirus", "anellovirus"]

CAT_DISPLAY = {
    "lactobacillus_phage": "Lactobacillus\nphage",
    "gardnerella_phage":   "Gardnerella\nphage",
    "megasphaera_phage":   "Megasphaera\nphage",
    "fannyhessea_phage":   "Fannyhessea\nphage",
    "sneathia_phage":      "Sneathia\nphage",
    "hpv":                 "HPV",
    "herpesvirus":         "Herpesvirus",
    "anellovirus":         "Anellovirus",
}

# CST display names
CST_LABELS = {
    "UC115_V2": "CST-I",
    "UC093_V3": "CST-IV",
}

CST_COLORS = {
    "UC115_V2": "#3182bd",
    "UC093_V3": "#e6550d",
}


def add_panel_label(ax, label, x=-0.08, y=1.05, fontsize=14):
    """Add a bold panel label (A, B, C, ...) to an axes."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold", va="bottom", ha="right")
