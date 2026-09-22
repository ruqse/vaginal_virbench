#!/usr/bin/env python3
"""Recompute the cascade tables that depend on the multi-evidence ranking,
using the LOCKED 13-sample pooled benchmark (secondary_benchmark_pooled13).

Outputs (to results/test_real/secondary_benchmark_pooled13/cascade/):
  - q2_prophage_detection.tsv      (-> Table S13)
  - q4_group_tests.tsv             (-> Table S13 caption; marker-vs-sequence tests)
  - cross_track_panelA.tsv         (-> Table S22 Panel A, multi-evidence swapped)
  - cross_track_panelB.tsv         (-> Table S22 Panel B, 14 tools, BH-adjusted)
  - cross_track_panelC.tsv         (-> Table S22 Panel C, top-5 by Track A L1500)
  - SUMMARY.txt                    (the headline §7/§8 numbers, for the manuscript)

Pure recompute from existing artifacts; no tool execution.

Use --only-cross-track to refresh Table S22 and its cascade inputs after a
Track B scoring change, without recomputing the Q4 or evidence-line analyses.
"""
from __future__ import annotations
import argparse, csv, gzip, sys
from itertools import combinations
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr, rankdata

ROOT = Path(__file__).resolve().parents[2]
POOL = ROOT / "results/test_real/secondary_benchmark_pooled13"
OUT = POOL / "cascade"
OUT.mkdir(parents=True, exist_ok=True)

MARKER = {"geNomad", "VirSorter2", "VirSorter", "VIBRANT"}
# Q2 marker-vs-sequence contrast excludes Jaeger & Seeker (length incompat in the
# original design) and Sourmash (degenerate: emits ~2 contigs).
Q2_EXCLUDE = {"Jaeger", "Seeker", "Sourmash"}


def exact_rank_permutation(marker, sequence):
    """Enumerate group assignments using midranks, including tied observations.

    Two-sided p is twice the smaller inclusive tail, capped at one. Unlike
    scipy's Mann-Whitney method='exact', this null distribution retains ties.
    """
    n = len(marker)
    ranks = rankdata(list(marker) + list(sequence), method="average")
    offset = n * (n + 1) / 2
    observed = float(ranks[:n].sum() - offset)
    null = np.asarray([sum(ranks[list(idx)]) - offset
                       for idx in combinations(range(len(ranks)), n)])
    p_less = float(np.mean(null <= observed))
    p_greater = float(np.mean(null >= observed))
    return observed, min(1.0, 2 * min(p_less, p_greater)), p_greater, p_less

# ---------------------------------------------------------------- Q2 / S13
def recompute_q2():
    # category per contig
    cat = {}
    with open(POOL / "pooled_benchmark_ready_ground_truth.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cat[r["contig_id"]] = r["category"]
    n_pro = sum(1 for c in cat.values() if c == "prophage")
    n_free = sum(1 for c in cat.values() if c == "free_phage")
    # per-tool counts of viral calls on prophage / free-phage
    pro_tp = {}; free_pos = {}; types = {}
    with gzip.open(POOL / "pooled_tool_predictions.tsv.gz", "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            t = r["tool"]; c = cat.get(r["contig_id"])
            types.setdefault(t, "marker" if t in MARKER else "sequence")
            yp = r["y_pred"] == "1"
            if c == "prophage":
                pro_tp[t] = pro_tp.get(t, 0) + (1 if yp else 0)
            elif c == "free_phage":
                free_pos[t] = free_pos.get(t, 0) + (1 if yp else 0)
    rows = []
    for t in sorted(pro_tp, key=lambda x: -(pro_tp[x] / max(n_pro, 1))):
        rows.append({
            "tool": t, "type": types[t],
            "prophage_tp": pro_tp[t], "prophage_total": n_pro,
            "prophage_recall": round(pro_tp[t] / n_pro, 4),
            "freephage_pos": free_pos.get(t, 0), "freephage_total": n_free,
            "freephage_rate": round(free_pos.get(t, 0) / max(n_free, 1), 4),
            "included_in_group_mean": "no" if t in Q2_EXCLUDE else "yes",
        })
    with open(OUT / "q2_prophage_detection.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    # group means (exclude Jaeger/Seeker/Sourmash; see included_in_group_mean)
    def gmean(grp, key):
        vals = [r[key] for r in rows
                if r["type"] == grp and r["included_in_group_mean"] == "yes"]
        return (sum(vals) / len(vals), len(vals)) if vals else (0.0, 0)
    mk_pro, n_mk = gmean("marker", "prophage_recall")
    sq_pro, n_sq = gmean("sequence", "prophage_recall")
    mk_fp, _ = gmean("marker", "freephage_rate")
    sq_fp, _ = gmean("sequence", "freephage_rate")
    # Q4 group contrasts. Unit of analysis is the tool, so the group sizes are
    # 4 and 7 and the exact (permutation) null is used rather than the normal
    # approximation. Two-sided: the free-phage direction was not pre-specified.
    def gvals(grp, key):
        return [r[key] for r in rows
                if r["type"] == grp and r["included_in_group_mean"] == "yes"]
    tests = []
    for key, label in (("prophage_recall", "Tier 1 prophage recall"),
                       ("freephage_rate", "free-phage positive rate")):
        mv, sv = gvals("marker", key), gvals("sequence", key)
        u, p, p_gt, p_lt = exact_rank_permutation(mv, sv)
        tests.append({
            "contrast": label, "metric": key,
            "n_marker": len(mv), "n_sequence": len(sv),
            "marker_mean": round(sum(mv) / len(mv), 4),
            "marker_min": min(mv), "marker_max": max(mv),
            "sequence_mean": round(sum(sv) / len(sv), 4),
            "sequence_min": min(sv), "sequence_max": max(sv),
            "U_marker": round(float(u), 1), "p_two_sided": round(float(p), 4),
            "p_one_sided_marker_greater": round(float(p_gt), 4),
            "p_one_sided_marker_less": round(float(p_lt), 4),
            "test": "exact rank permutation with ties; unit = tool; two-sided = twice smaller inclusive tail",
        })
    with open(OUT / "q4_group_tests.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(tests[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(tests)
    return dict(n_pro=n_pro, n_free=n_free, mk_pro=mk_pro, sq_pro=sq_pro,
                mk_fp=mk_fp, sq_fp=sq_fp, n_mk=n_mk, n_sq=n_sq, tests=tests)

# ---------------------------------------------------------------- Cross-track / S22
S22_TEMPLATE = ROOT / "results/tables/templates/table_s22_cross_track_spearman.template.tsv"

def recompute_cross_track():
    # new multi-evidence MCC per tool
    me = {}
    with open(POOL / "overall_metrics.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            me[r["tool"]] = float(r["MCC"])
    # Preserve the other published contexts, but refresh Track B from its
    # source table. Copying its old S22 column would retain stale scoring.
    contexts = ["Track A (L1500 MCC)", "Track A (L3000 MCC)",
                "Track B (10x coverage, CST-I MCC)",
                "Multi-evidence benchmark (MCC)", "Track C (MCC)"]
    # Tool order comes from the S22 template; every number is recomputed from
    # current sources below (S22 is never read back as its own input).
    panelA = {}
    with open(S22_TEMPLATE) as fh:
        in_a = False
        for line in fh:
            if line.startswith("tool\tTrack A (L1500"):
                in_a = True; continue
            if in_a:
                if not line.strip() or line.startswith("#"):
                    break
                panelA[line.rstrip("\n").split("\t")[0]] = [None] * 5
    if set(panelA) - set(me):
        raise ValueError("Multi-evidence metrics must contain one current MCC per S22 tool")
    with open(ROOT / "results/tables/table_s8_track_b_metrics.tsv") as fh:
        track_b_rows = [r for r in csv.DictReader(fh, delimiter="\t")
                        if r["background"] == "UC115_V2" and float(r["coverage"]) == 10]
    track_b = {r["tool"]: float(r["MCC"]) for r in track_b_rows}
    if len(track_b_rows) != len(panelA) or set(track_b) != set(panelA):
        raise ValueError("Table S8 must contain one UC115_V2 10x MCC per S22 tool")
    with open(ROOT / "results/test_real/coassembly/benchmark/overall_metrics.tsv") as fh:
        track_c = {r["tool"]: float(r["MCC"]) for r in csv.DictReader(fh, delimiter="\t")}
    if set(track_c) != set(panelA):
        raise ValueError("Track C must contain one current MCC per S22 tool")
    with open(ROOT / "results/tables/table_s7_track_a_overall_metrics_by_length.tsv") as fh:
        track_a = {(r["tool"], int(r["length_bp"])): float(r["MCC"])
                   for r in csv.DictReader(fh, delimiter="\t")}
    # All assembly contexts use current sources, not frozen S22 values.
    for tool in panelA:
        panelA[tool][0] = track_a[tool, 1500]
        panelA[tool][1] = track_a[tool, 3000]
        panelA[tool][2] = track_b[tool]
        panelA[tool][3] = me[tool]
        panelA[tool][4] = track_c[tool]
    tools = list(panelA.keys())
    mat = np.array([panelA[t] for t in tools])  # tools x 5 contexts

    def pairwise(idx_tools):
        sub = np.array([panelA[t] for t in idx_tools])
        rows = []
        for i in range(5):
            for j in range(i + 1, 5):
                rho, p = spearmanr(sub[:, i], sub[:, j])
                rows.append([contexts[i], contexts[j], len(idx_tools),
                             round(rho, 3), p])
        # BH adjust
        ps = [r[4] for r in rows]
        order = np.argsort(ps)
        m = len(ps); adj = [0.0] * m
        prev = 1.0
        for rank, oi in enumerate(reversed(order)):
            k = m - rank
            val = min(prev, ps[oi] * m / k)
            adj[oi] = val; prev = val
        for r, a in zip(rows, adj):
            r.append(round(r[4], 6)); r[4] = round(r[4], 6); r.append(round(a, 4))
        return rows

    top5 = ["geNomad", "ViraLM", "VirSorter2", "TransGINmer", "VirSorter"]
    B = pairwise(tools)
    C = pairwise(top5)

    # write Panel A
    with open(OUT / "cross_track_panelA.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["tool"] + contexts)
        for t in tools:
            w.writerow([t] + [round(x, 4) for x in panelA[t]])
    # write Panel B/C
    for name, rows in [("cross_track_panelB.tsv", B), ("cross_track_panelC.tsv", C)]:
        with open(OUT / name, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["context_1", "context_2", "n_tools", "spearman_rho",
                        "spearman_p", "spearman_p_dup", "spearman_p_BH"])
            w.writerows(rows)
    # summary stats
    def stats(rows):
        rhos = sorted(r[3] for r in rows)
        sig = sum(1 for r in rows if r[6] < 0.05)
        return min(rhos), max(rhos), float(np.median(rhos)), sig, len(rows)
    return dict(B=stats(B), C=stats(C),
                me_rank=sorted(me.items(), key=lambda kv: -kv[1]))


def write_cross_track_table():
    """Publish only S22 from the refreshed cascade panels."""
    table = ROOT / "results/tables/table_s22_cross_track_spearman.tsv"
    preamble = S22_TEMPLATE.read_text().split("## Panel A:", 1)[0]
    titles = {
        "A": "per-tool MCC by context",
        "B": "pairwise Spearman rank correlations",
        "C": "restricted to top-5 Track A L1500 tools "
             "(geNomad, ViraLM, VirSorter2, TransGINmer, VirSorter)",
    }
    with table.open("w", newline="") as fh:
        fh.write(preamble)
        for panel, title in titles.items():
            if panel != "A":
                fh.write("\n")
            fh.write(f"## Panel {panel}: {title}\n")
            with (OUT / f"cross_track_panel{panel}.tsv").open() as source:
                reader = csv.DictReader(source, delimiter="\t")
                fields = [field for field in reader.fieldnames if field != "spearman_p_dup"]
                writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                                        extrasaction="ignore")
                writer.writeheader()
                writer.writerows(reader)


def refresh_cross_track_summary(ct):
    """Update the two cross-track lines without changing Q4/S11 summaries."""
    bmin, bmax, bmed, bsig, bn = ct["B"]
    cmin, cmax, cmed, csig, cn = ct["C"]
    updated = {
        "Panel B (": f"Panel B (14 tools, 10 pairs): rho {bmin:.2f}-{bmax:.2f}  "
                     f"median {bmed:.2f}  {bsig}/{bn} significant (q<0.05)",
        "Panel C (": f"Panel C (top-5 Track A): rho {cmin:.2f}-{cmax:.2f}  "
                     f"median {cmed:.2f}  {csig}/{cn} significant",
    }
    summary = OUT / "SUMMARY.txt"
    if summary.exists():
        lines = summary.read_text().splitlines()
        for prefix, replacement in updated.items():
            matches = [i for i, line in enumerate(lines) if line.startswith(prefix)]
            if len(matches) != 1:
                raise ValueError(f"Expected one {prefix} line in {summary}")
            lines[matches[0]] = replacement
        summary.write_text("\n".join(lines) + "\n")
    print("\n".join(updated.values()))

# ---------------------------------------------------------------- S11 tier sensitivity
def recompute_s11():
    variants = [("primary", POOL),
                ("no_E1b", POOL.parent / "secondary_benchmark_pooled13_no_e1b"),
                ("merged_E1E2", POOL.parent / "secondary_benchmark_pooled13_merged_e1e2"),
                ("no_E5", POOL.parent / "secondary_benchmark_pooled13_no_e5")]
    data = {}
    for name, d in variants:
        f = d / "overall_metrics.tsv"
        if not f.exists():
            print(f"[S11] WARN missing {f}"); continue
        mccs = {}
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                mccs[r["tool"]] = float(r["MCC"])
        ranked = sorted(mccs, key=lambda t: -mccs[t])
        data[name] = (mccs, {t: i + 1 for i, t in enumerate(ranked)})
    if "primary" not in data:
        return None
    tools = sorted(data["primary"][0], key=lambda t: -data["primary"][0][t])
    rows = []
    for t in tools:
        row = {"tool": t}
        for name, _ in variants:
            if name in data:
                row[f"MCC_{name}"] = round(data[name][0][t], 4)
                row[f"rank_{name}"] = data[name][1][t]
        rows.append(row)
    with open(OUT / "s11_multi_evidence_sensitivity.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    # geNomad rank-1 stability
    g = {name: data[name][1].get("geNomad") for name, _ in variants if name in data}
    return g

# ---------------------------------------------------------------- main
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only-cross-track", action="store_true",
                        help="Refresh S22 and cross-track cascade outputs only")
    args = parser.parse_args()
    if args.only_cross_track:
        ct = recompute_cross_track()
        write_cross_track_table()
        refresh_cross_track_summary(ct)
        print("Updated Table S22 and cross-track cascade panels")
        return
    q2 = recompute_q2()
    ct = recompute_cross_track()
    s11 = recompute_s11()
    with open(OUT / "SUMMARY.txt", "w") as fh:
        def pr(*a):
            s = " ".join(str(x) for x in a); print(s); fh.write(s + "\n")
        pr("=== Q2 / Table S13 (13-sample pooled, >=1500bp) ===")
        pr(f"prophage contigs: {q2['n_pro']}   free-phage contigs: {q2['n_free']}")
        pr(f"prophage recall  marker={q2['mk_pro']*100:.1f}% (n={q2['n_mk']})  "
           f"sequence={q2['sq_pro']*100:.1f}% (n={q2['n_sq']})  "
           f"[excl Jaeger/Seeker/Sourmash]")
        pr(f"free-phage rate  marker={q2['mk_fp']*100:.1f}%  sequence={q2['sq_fp']*100:.1f}%  "
           f"(fold {q2['sq_fp']/max(q2['mk_fp'],1e-9):.1f}x)")
        for t in q2["tests"]:
            pr(f"  {t['contrast']:26s} marker {t['marker_mean']*100:.1f}% "
               f"[{t['marker_min']*100:.1f}-{t['marker_max']*100:.1f}]  "
               f"sequence {t['sequence_mean']*100:.1f}% "
               f"[{t['sequence_min']*100:.1f}-{t['sequence_max']*100:.1f}]  "
               f"U={t['U_marker']:.1f}  p={t['p_two_sided']:.3f} (exact, two-sided)")
        pr("")
        pr("=== Cross-track / Table S22 + §8 ===")
        bmin, bmax, bmed, bsig, bn = ct["B"]
        cmin, cmax, cmed, csig, cn = ct["C"]
        pr(f"Panel B (14 tools, 10 pairs): rho {bmin:.2f}-{bmax:.2f}  median {bmed:.2f}  "
           f"{bsig}/{bn} significant (q<0.05)")
        pr(f"Panel C (top-5 Track A): rho {cmin:.2f}-{cmax:.2f}  median {cmed:.2f}  "
           f"{csig}/{cn} significant")
        pr("")
        pr("=== S11 tier-sensitivity: geNomad rank by variant ===")
        pr(f"  {s11}" if s11 else "  (variant runs not found)")
        pr("")
        pr("=== Multi-evidence ranking (new) ===")
        for i, (t, m) in enumerate(ct["me_rank"], 1):
            pr(f"  {i:2d}. {t:14s} {m:.4f}")
    print("\nWrote cascade outputs to", OUT)


if __name__ == "__main__":
    main()
