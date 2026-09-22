#!/usr/bin/env python3
"""
track_c_gmm_sensitivity.py — Gaussian-mixture-model alternative to the
hand-picked R thresholds used to define Track C TP / TN / dark-matter labels.

The Track C primary benchmark uses conservative round-number cutoffs aligned
with the bimodal log10(R) distribution (Methods § "VLP-enrichment benchmark";
Limitations paragraph). A 3 x 3 grid of alternative round-number cutoffs
(track_c_threshold_sensitivity.py) shows that rankings are robust across
nine threshold combinations. This script tests the ranking against a
*methodologically distinct* alternative: model-based valley-finding via a
2-component Gaussian mixture model on log10(R), with per-contig posterior
probabilities defining the TP / TN candidate sets.

Pipeline:
  1. Pool per-sample enrichment files; drop R = 0 and non-finite values.
  2. Fit a 2-component GMM on log10(R).
  3. Identify the high-R component (mean) and the low-R component.
  4. Define posterior cutoffs:
        TP candidate: P(high-R | R) > P_HIGH  (default 0.95)
        TN candidate: P(high-R | R) < P_LOW   (default 0.05)
        excluded:     P_LOW <= P(high-R | R) <= P_HIGH
  5. Apply the same orthogonal CheckV / Kraken2 conjunctive rules used by the
     primary benchmark (TP must also have CheckV viral signal; TN must also
     have Kraken2 bacterial domain and CheckV viral_genes = 0).
  6. Score all 14 tools against the GMM-derived TP/TN set.
  7. Compare the resulting ranking to the primary Track C ranking.

Outputs (under --output-dir):
  - track_c_gmm_fit.tsv           (component means/variances/weights;
                                   implied R thresholds at P_LOW and P_HIGH)
  - track_c_gmm_counts.tsv        (TP / TN / dark / excluded under GMM rules)
  - track_c_gmm_per_tool.tsv      (per-tool MCC, recall, precision, F1)
  - track_c_gmm_vs_default.tsv    (side-by-side ranking comparison;
                                   Spearman rho)
  - track_c_gmm_relabel.tsv       (per-contig: default label vs GMM label;
                                   only contigs that change category)

Re-uses helpers from track_c_threshold_sensitivity.py + build_rca_ground_truth.py.
"""

import argparse
import csv
import importlib.util
import sys
from collections import Counter
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Module imports (mirrors track_c_threshold_sensitivity.py)
# ---------------------------------------------------------------------------

_SCRIPTS = Path(__file__).resolve().parent.parent
_BUILD_GT = _SCRIPTS / "03_assembly" / "build_rca_ground_truth.py"
_EVAL_RCA = _SCRIPTS / "03_assembly" / "evaluate_rca_benchmark.py"


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_bgt = _load_module("build_rca_ground_truth", _BUILD_GT)
_eval = _load_module("evaluate_rca_benchmark", _EVAL_RCA)

PATIENTS = _bgt.PATIENTS
BACTERIAL_DOMAINS = _bgt.BACTERIAL_DOMAINS
has_viral_signal = _bgt.has_viral_signal
load_enrichment = _bgt.load_enrichment
load_checkv = _bgt.load_checkv
load_kraken2 = _bgt.load_kraken2
classify_domain = _bgt.classify_domain

TOOL_THRESHOLDS = _eval.TOOL_THRESHOLDS
evaluate_tool = _eval.evaluate_tool
TOOLS = list(TOOL_THRESHOLDS.keys())

# import sibling sensitivity script for the tool-prediction loader
_SENS = _SCRIPTS / "07_evaluation" / "track_c_threshold_sensitivity.py"
_sens = _load_module("track_c_threshold_sensitivity", _SENS)
load_all_tool_predictions = _sens.load_all_tool_predictions

# Default Track C labels (from the primary benchmark)
DEFAULT_TP_R = 10.0
DEFAULT_TN_BAND = (0.5, 2.0)


# ---------------------------------------------------------------------------
# GMM fitting
# ---------------------------------------------------------------------------

def fit_gmm_log10R(enrichment, n_components=3, random_state=0):
    """Fit an n-component GMM on log10(enrichment_R).

    Empirically (Track C, 8,164 contigs >= 2,500 bp with R > 0 and finite):
      k=2: BIC = 21875.6
      k=3: BIC = 21078.9   <-- best fit, manuscript's "bimodal" is actually trimodal:
                                depleted (R~0.01), baseline (R~0.6), enriched (R~129)
      k=4: BIC = 21120.1
      k=5: BIC = 21141.0
    The 3-component model maps cleanly to the manuscript's intended categories.

    Contigs with R <= 0 or non-finite (R == 0 or shotgun_depth == 0 -> R = inf)
    are excluded from the fit. Inf-R contigs are handled separately downstream.
    """
    from sklearn.mixture import GaussianMixture

    contig_ids = []
    r_values = []
    for cid in sorted(enrichment.keys()):
        r = enrichment[cid].get("enrichment_R", 0.0)
        contig_ids.append(cid)
        r_values.append(r)
    r_values = np.asarray(r_values, dtype=float)
    contig_ids = np.asarray(contig_ids)

    # Mask: positive, finite R (excludes R = 0 and R = inf)
    mask = np.isfinite(r_values) & (r_values > 0)
    log10R_fit = np.log10(r_values[mask]).reshape(-1, 1)

    gmm = GaussianMixture(n_components=n_components, covariance_type="full",
                          n_init=10, random_state=random_state, max_iter=500)
    gmm.fit(log10R_fit)

    # Sort components by mean (low -> high)
    means = gmm.means_.ravel()
    order = np.argsort(means)
    high_idx = int(order[-1])     # highest-mean = enriched component
    low_idx = int(order[0])       # lowest-mean = depleted component
    middle_idx = int(order[1]) if n_components >= 3 else None

    summary = {
        "n_components": n_components,
        "n_contigs_fit": int(mask.sum()),
        "n_contigs_R0": int((r_values == 0).sum()),
        "n_contigs_Rinf": int(np.isinf(r_values).sum()),
        "BIC": float(gmm.bic(log10R_fit)),
        "converged": bool(gmm.converged_),
        "log_likelihood": float(gmm.score(log10R_fit) * len(log10R_fit)),
        "low_mean_log10R": float(means[low_idx]),
        "low_std_log10R": float(np.sqrt(gmm.covariances_[low_idx].ravel()[0])),
        "low_weight": float(gmm.weights_[low_idx]),
        "low_R_mean": float(10 ** means[low_idx]),
        "high_mean_log10R": float(means[high_idx]),
        "high_std_log10R": float(np.sqrt(gmm.covariances_[high_idx].ravel()[0])),
        "high_weight": float(gmm.weights_[high_idx]),
        "high_R_mean": float(10 ** means[high_idx]),
    }
    if middle_idx is not None:
        summary.update({
            "middle_mean_log10R": float(means[middle_idx]),
            "middle_std_log10R": float(np.sqrt(gmm.covariances_[middle_idx].ravel()[0])),
            "middle_weight": float(gmm.weights_[middle_idx]),
            "middle_R_mean": float(10 ** means[middle_idx]),
        })

    return gmm, r_values, contig_ids, mask, high_idx, low_idx, middle_idx, summary


def posterior_high_component(gmm, log10R_values, high_idx):
    """Return P(high-R component | R) for each value of log10R."""
    proba = gmm.predict_proba(log10R_values.reshape(-1, 1))
    return proba[:, high_idx]


def find_R_at_posterior(gmm, high_idx, p_target, search_range=(-3, 4)):
    """Find the R value (linear scale) at which P(high | R) = p_target.

    Uses 1D bisection on log10(R) space.
    """
    lo, hi = search_range
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        p = gmm.predict_proba(np.array([[mid]]))[0, high_idx]
        if p < p_target:
            lo = mid
        else:
            hi = mid
    return float(10 ** mid)


# ---------------------------------------------------------------------------
# GMM-based ground truth construction
# ---------------------------------------------------------------------------

def build_gt_from_gmm(enrichment, checkv, kraken, domains, gmm,
                      high_idx, low_idx, middle_idx,
                      p_high=0.95, p_baseline=0.95):
    """Re-classify contigs using 3-component GMM posterior cutoffs.

    The 3 components correspond to (lowest -> highest mean):
      low    = depleted     (R << 1)  -- not used as TN by itself
      middle = bacterial baseline (R ~ 1) -- TN candidate
      high   = enriched     (R >> 1)  -- TP / dark candidate

    Decision rule:
      - TP candidate:        P(high | R) > p_high   AND CheckV viral signal
      - TN candidate:        P(middle | R) > p_baseline AND not viral
                                                       AND domain in BACTERIAL_DOMAINS
      - dark-matter candidate: P(high | R) > p_high  AND no viral signal
                                                       AND no host genes
                                                       AND no Kraken2 classification
      - excluded: everything else (including R = 0, R = inf without conjunctive
                                   evidence, low-mode contigs, and the GMM
                                   "I don't know" zones)

    R = inf contigs (zero shotgun depth) are treated as P(high) = 1.0.
    R = 0 contigs are treated as P(low) = 1.0 -> excluded.
    """
    gt = {}
    counts = Counter()
    relabel_records = []

    for contig_id in sorted(enrichment.keys()):
        e = enrichment[contig_id]
        r = e["enrichment_R"]
        cv = checkv.get(contig_id, {
            "checkv_quality": "Not-determined",
            "viral_genes": 0,
            "host_genes": 0,
            "provirus": "No",
        })
        kr = kraken.get(contig_id, {"classified": False, "taxon_name": "unclassified"})
        domain = domains.get(contig_id, "unknown")
        is_viral, _ = has_viral_signal(cv)

        # Posterior probabilities under the GMM
        if r == 0:
            p_hi = 0.0
            p_md = 0.0
        elif np.isinf(r):
            p_hi = 1.0    # zero shotgun depth -> definitely enriched mode
            p_md = 0.0
        elif np.isfinite(r) and r > 0:
            proba = gmm.predict_proba(np.array([[np.log10(r)]]))[0]
            p_hi = float(proba[high_idx])
            p_md = float(proba[middle_idx]) if middle_idx is not None else 0.0
        else:
            p_hi = float("nan")
            p_md = float("nan")

        # Classification
        if p_hi > p_high:
            if is_viral:
                gt[contig_id] = {
                    "label": "viral",
                    "contig_length": e["contig_length"],
                    "checkv_quality": cv.get("checkv_quality", "NA"),
                    "viral_genes": cv.get("viral_genes", 0),
                    "host_genes": cv.get("host_genes", 0),
                    "provirus": cv.get("provirus", "No"),
                }
                counts["TP"] += 1
                gmm_label = "TP"
            elif (cv.get("viral_genes", 0) == 0
                  and cv.get("host_genes", 0) == 0
                  and not kr["classified"]):
                counts["dark_matter"] += 1
                gmm_label = "dark_matter"
            else:
                counts["excluded"] += 1
                gmm_label = "excluded"
        elif p_md > p_baseline:
            if not is_viral and domain in BACTERIAL_DOMAINS:
                gt[contig_id] = {
                    "label": "negative",
                    "contig_length": e["contig_length"],
                    "checkv_quality": cv.get("checkv_quality", "NA"),
                    "viral_genes": cv.get("viral_genes", 0),
                    "host_genes": cv.get("host_genes", 0),
                    "provirus": cv.get("provirus", "No"),
                }
                counts["TN"] += 1
                gmm_label = "TN"
            else:
                counts["excluded"] += 1
                gmm_label = "excluded"
        else:
            counts["excluded"] += 1
            gmm_label = "excluded"

        relabel_records.append({
            "contig_id": contig_id,
            "enrichment_R": r if np.isfinite(r) else float("inf"),
            "log10R": np.log10(r) if (np.isfinite(r) and r > 0) else float("nan"),
            "P_high": p_hi,
            "P_middle": p_md,
            "gmm_label": gmm_label,
        })

    return gt, counts, relabel_records


# ---------------------------------------------------------------------------
# Default-label loader (for ranking comparison)
# ---------------------------------------------------------------------------

def load_default_labels(rca_gt_path: Path):
    """Read the primary-benchmark labels from rca_ground_truth.tsv (label col)."""
    out = {}
    with open(rca_gt_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            out[row["contig_id"]] = row["label"]
    return out


# ---------------------------------------------------------------------------
# Spearman rank correlation (no scipy dependency)
# ---------------------------------------------------------------------------

def spearman_rho(x, y):
    """Spearman rank correlation between two equal-length lists."""
    def ranks(a):
        order = np.argsort(a)
        r = np.empty_like(order, dtype=float)
        r[order] = np.arange(len(a))
        # ties: average ranks (use scipy if available; else acceptable approximation)
        return r
    rx, ry = ranks(np.asarray(x)), ranks(np.asarray(y))
    return float(np.corrcoef(rx, ry)[0, 1])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--coassembly-dir", required=True, type=Path)
    parser.add_argument("--checkv-dir", required=True, type=Path)
    parser.add_argument("--kraken2-output", required=True, type=Path)
    parser.add_argument("--kraken2-report", type=Path, default=None)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--sample-id", default="RCA_master")
    parser.add_argument("--rca-ground-truth", required=True, type=Path,
                        help="Path to the primary-benchmark rca_ground_truth.tsv")
    parser.add_argument("--default-mcc", type=Path, default=None,
                        help="Path to track_c_sensitivity_mcc.tsv (or table3); "
                             "if provided, used to extract default-label MCC ranking")
    parser.add_argument("--n-components", type=int, default=3,
                        help="Number of GMM components (default 3 = depleted/baseline/enriched). "
                             "BIC strongly favours k=3 over k=2 on this dataset.")
    parser.add_argument("--p-high", type=float, default=0.95,
                        help="P(enriched component) cutoff for TP candidates")
    parser.add_argument("--p-baseline", type=float, default=0.95,
                        help="P(baseline component) cutoff for TN candidates")
    parser.add_argument("--output-dir", "-o", required=True, type=Path)
    parser.add_argument("--random-state", type=int, default=0)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Track C — GMM-based threshold sensitivity")
    print("=" * 70)

    # --- Load data once -----------------------------------------------------
    print("\nLoading enrichment / CheckV / Kraken2 ...")
    enrichment = load_enrichment(args.coassembly_dir)
    checkv = load_checkv(args.checkv_dir)
    kraken = load_kraken2(args.kraken2_output)
    domains = classify_domain(args.kraken2_output, args.kraken2_report)
    print(f"  enrichment: {len(enrichment)}, checkv: {len(checkv)}, "
          f"kraken2: {len(kraken)}, domains: {len(domains)}")

    # --- Fit GMM ------------------------------------------------------------
    print(f"\nFitting {args.n_components}-component GMM on log10(R) ...")
    gmm, r_values, contig_ids, mask, high_idx, low_idx, middle_idx, fit = \
        fit_gmm_log10R(enrichment, n_components=args.n_components,
                       random_state=args.random_state)
    print(f"  n contigs fit (R > 0 & finite):  {fit['n_contigs_fit']}")
    print(f"  n contigs R == 0  (excluded):    {fit['n_contigs_R0']}")
    print(f"  n contigs R == inf (treated TP): {fit['n_contigs_Rinf']}")
    print(f"  low (depleted)  : log10R={fit['low_mean_log10R']:+7.3f}  "
          f"(R={fit['low_R_mean']:.3f}), std={fit['low_std_log10R']:.3f}, "
          f"w={fit['low_weight']:.3f}")
    if middle_idx is not None:
        print(f"  middle (baseline): log10R={fit['middle_mean_log10R']:+7.3f}  "
              f"(R={fit['middle_R_mean']:.3f}), std={fit['middle_std_log10R']:.3f}, "
              f"w={fit['middle_weight']:.3f}")
    print(f"  high (enriched) : log10R={fit['high_mean_log10R']:+7.3f}  "
          f"(R={fit['high_R_mean']:.3f}), std={fit['high_std_log10R']:.3f}, "
          f"w={fit['high_weight']:.3f}")
    print(f"  BIC={fit['BIC']:.1f}, converged={fit['converged']}")

    # --- Implied R thresholds at posterior cutoffs --------------------------
    R_at_phigh = find_R_at_posterior(gmm, high_idx, args.p_high)
    print(f"\n  Implied R at P(high) = {args.p_high}: R = {R_at_phigh:.3f}")
    print(f"  (compare to default Track C: TP at R > {DEFAULT_TP_R}, "
          f"TN band {DEFAULT_TN_BAND})")

    # --- Build GMM ground truth --------------------------------------------
    print("\nBuilding GMM-based ground truth ...")
    gt_gmm, counts_gmm, relabel = build_gt_from_gmm(
        enrichment, checkv, kraken, domains, gmm,
        high_idx, low_idx, middle_idx,
        p_high=args.p_high, p_baseline=args.p_baseline,
    )
    print(f"  TP:      {counts_gmm['TP']}")
    print(f"  TN:      {counts_gmm['TN']}")
    print(f"  dark:    {counts_gmm['dark_matter']}")
    print(f"  excluded:{counts_gmm['excluded']}")

    # --- Compare to default labels ------------------------------------------
    default_labels = load_default_labels(args.rca_ground_truth)
    n_changed = 0
    flip_records = []
    for rec in relabel:
        cid = rec["contig_id"]
        default_label = default_labels.get(cid, "excluded")
        # Map default labels to {TP, TN, dark_matter, excluded}
        default_mapped = {"TP": "TP", "TN": "TN",
                          "dark_matter": "dark_matter",
                          "excluded": "excluded"}.get(default_label, "excluded")
        rec["default_label"] = default_mapped
        if rec["gmm_label"] != default_mapped:
            n_changed += 1
            flip_records.append(rec)
    print(f"\n  contigs changing category vs default: {n_changed}")

    # --- Load tool predictions and score against GMM ground truth ----------
    print(f"\nLoading tool predictions from {args.results_dir} ...")
    tool_preds = load_all_tool_predictions(args.results_dir, args.sample_id,
                                           gt_ids=set(enrichment))

    print("\nScoring 14 tools against GMM-derived TP/TN set ...")
    per_tool_gmm = {}
    for tool in TOOLS:
        threshold = TOOL_THRESHOLDS[tool]
        m = evaluate_tool(tool, tool_preds[tool], gt_gmm, threshold)
        per_tool_gmm[tool] = m
        print(f"  {tool:14s}  MCC={m['MCC']:.4f}  "
              f"P={m['precision']:.3f}  R={m['recall']:.3f}  F1={m['F1']:.3f}")

    # --- Default-ranking comparison ----------------------------------------
    # Score the same tools against the *default* ground truth for a clean comparison.
    print("\nScoring 14 tools against default-label TP/TN set (recomputed) ...")
    default_gt = {}
    for cid, lbl in default_labels.items():
        if lbl in ("TP", "TN"):
            e = enrichment.get(cid, {})
            cv = checkv.get(cid, {})
            default_gt[cid] = {
                "label": "viral" if lbl == "TP" else "negative",
                "contig_length": e.get("contig_length", 0),
                "checkv_quality": cv.get("checkv_quality", "NA"),
                "viral_genes": cv.get("viral_genes", 0),
                "host_genes": cv.get("host_genes", 0),
                "provirus": cv.get("provirus", "No"),
            }
    per_tool_default = {}
    for tool in TOOLS:
        threshold = TOOL_THRESHOLDS[tool]
        per_tool_default[tool] = evaluate_tool(tool, tool_preds[tool],
                                               default_gt, threshold)

    # Spearman rho between GMM and default rankings
    mcc_gmm = [per_tool_gmm[t]["MCC"] for t in TOOLS]
    mcc_def = [per_tool_default[t]["MCC"] for t in TOOLS]
    rho = spearman_rho(mcc_gmm, mcc_def)
    print(f"\n  Spearman rho (GMM-MCC vs default-MCC, n=14): {rho:.4f}")

    # Top-5 set comparison
    top5_gmm = set(sorted(TOOLS, key=lambda t: -per_tool_gmm[t]["MCC"])[:5])
    top5_def = set(sorted(TOOLS, key=lambda t: -per_tool_default[t]["MCC"])[:5])
    print(f"  top-5 (GMM):     {sorted(top5_gmm)}")
    print(f"  top-5 (default): {sorted(top5_def)}")
    print(f"  intersection:    {sorted(top5_gmm & top5_def)} "
          f"({len(top5_gmm & top5_def)}/5)")

    # --- Write outputs ------------------------------------------------------
    def write_tsv(rows, path, fields):
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            for row in rows:
                w.writerow(row)
        print(f"  wrote {path}")

    # GMM fit summary
    fit_row = {**fit,
               "p_high": args.p_high, "p_baseline": args.p_baseline,
               "R_at_p_high": R_at_phigh,
               "default_TP_R": DEFAULT_TP_R,
               "default_TN_low": DEFAULT_TN_BAND[0],
               "default_TN_high": DEFAULT_TN_BAND[1]}
    write_tsv(
        [fit_row],
        args.output_dir / "track_c_gmm_fit.tsv",
        list(fit_row.keys()),
    )

    # GMM counts
    write_tsv(
        [{"TP": counts_gmm["TP"], "TN": counts_gmm["TN"],
          "dark_matter": counts_gmm["dark_matter"],
          "excluded": counts_gmm["excluded"],
          "n_changed_vs_default": n_changed}],
        args.output_dir / "track_c_gmm_counts.tsv",
        ["TP", "TN", "dark_matter", "excluded", "n_changed_vs_default"],
    )

    # Per-tool metrics (GMM)
    rows = []
    for tool in TOOLS:
        m = per_tool_gmm[tool]
        rows.append({"tool": tool,
                     "MCC": m["MCC"], "precision": m["precision"],
                     "recall": m["recall"], "F1": m["F1"],
                     "TP": m.get("TP", "NA"), "FP": m.get("FP", "NA"),
                     "TN": m.get("TN", "NA"), "FN": m.get("FN", "NA")})
    write_tsv(rows, args.output_dir / "track_c_gmm_per_tool.tsv",
              ["tool", "MCC", "precision", "recall", "F1",
               "TP", "FP", "TN", "FN"])

    # Side-by-side ranking comparison
    rows = []
    rank_gmm = sorted(TOOLS, key=lambda t: -per_tool_gmm[t]["MCC"])
    rank_def = sorted(TOOLS, key=lambda t: -per_tool_default[t]["MCC"])
    for i, t in enumerate(TOOLS):
        rows.append({
            "tool": t,
            "MCC_default": per_tool_default[t]["MCC"],
            "MCC_gmm": per_tool_gmm[t]["MCC"],
            "rank_default": rank_def.index(t) + 1,
            "rank_gmm": rank_gmm.index(t) + 1,
            "delta_MCC": per_tool_gmm[t]["MCC"] - per_tool_default[t]["MCC"],
            "spearman_rho_overall": rho,
        })
    rows.sort(key=lambda r: r["rank_default"])
    write_tsv(rows, args.output_dir / "track_c_gmm_vs_default.tsv",
              ["tool", "MCC_default", "MCC_gmm",
               "rank_default", "rank_gmm", "delta_MCC", "spearman_rho_overall"])

    # Per-contig relabel records (only changed)
    if flip_records:
        write_tsv(
            flip_records,
            args.output_dir / "track_c_gmm_relabel.tsv",
            ["contig_id", "enrichment_R", "log10R", "P_high", "P_middle",
             "default_label", "gmm_label"],
        )
    else:
        with open(args.output_dir / "track_c_gmm_relabel.tsv", "w") as f:
            f.write("# No contigs changed category between default and GMM labels.\n")
            f.write("contig_id\tenrichment_R\tlog10R\tP_high\tP_middle\tdefault_label\tgmm_label\n")
        print(f"  (no contigs flipped)")

    print("\nDone.")


if __name__ == "__main__":
    main()
