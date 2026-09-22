#!/usr/bin/env python3
"""Assemble all supplementary tables into a single multi-sheet workbook.

Post-renumber scheme (S1-S36 contiguous). Retired ranking tables are folded:
former S19 -> S18_rank sheet, former S21 -> S19_rank sheet. The five demoted
main tables now live as S33-S36 (+ Table2_full for recommendations).
- Comment-aware: strips leading '#' lines.
- Splits the multi-panel S22 (cross-track) file into S22_A/B/C sheets.
- Preserves the S23 GMM metadata block as an S23_metadata sheet.
- Keeps original TSVs untouched (this only reads them).
"""
import re
from pathlib import Path
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment

HERE = Path(__file__).resolve().parents[2] / "results" / "tables"   # published tables
OUT = HERE / "supplementary_tables.xlsx"

HEAD_FONT = Font(name="Arial", size=10, bold=True, color="FFFFFF")
BODY_FONT = Font(name="Arial", size=10)
HEAD_FILL = PatternFill("solid", fgColor="44546A")
WRAP = Alignment(vertical="top", wrap_text=True)

num_int = re.compile(r"-?\d+$")
num_flt = re.compile(r"(-?\d*\.\d+([eE][+-]?\d+)?|-?\d+[eE][+-]?\d+)$")

def conv(s):
    if s is None:
        return None
    t = s.strip()
    if t == "":
        return ""
    if num_int.fullmatch(t):
        try: return int(t)
        except ValueError: return t
    if num_flt.fullmatch(t):
        try: return float(t)
        except ValueError: return t
    return t

def read_tsv(path):
    """Rows of a plain TSV, with leading '#' comment lines stripped.

    Comment stripping is not cosmetic. add_sheet() styles row 1 as the header
    and freezes the pane beneath it, so a comment line surviving into the grid
    both mis-styles the sheet and pushes the real header into the body. Table
    S16 shipped that way until 2026-08-11; its comment block is reproduced in
    the Table S16 caption, which is what a reader actually sees.
    """
    rows = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        if line.strip() == "" and not rows:
            continue
        rows.append(line.rstrip("\r").split("\t"))
    return rows

def parse_panels(path):
    panels = {}
    cur = None
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("## Panel"):
            m = re.search(r"Panel ([A-D])", line)
            cur = m.group(1)
            panels[cur] = {"title": line[2:].strip(" #"), "rows": []}
            continue
        if line.startswith("#") or line.strip() == "" or cur is None:
            continue
        panels[cur]["rows"].append(line.split("\t"))
    return panels

def parse_meta_table(path):
    meta, table = [], []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            meta.append(re.sub(r"^#\s?", "", line))
        elif line.strip() != "":
            table.append(line.split("\t"))
    return meta, table

# Terminology normalisation for shipped headers. The five-way method taxonomy
# of Table S21 was called a "family" until 2026-09-02; a collaborator flagged
# the clash with taxonomic family (Tables S24, S25 carry viral_family /
# vmgc_family in the same workbook), so every shipped sheet says
# method_approach. Every results/tables TSV was renamed at the same time, so this
# map is normally a no-op; it is the safety net for a sheet re-copied from a
# results/ file, where the pipeline still writes the internal tool_family.
# Only these exact tokens are touched; viral_family and vmgc_family are
# deliberately absent from the map.
APPROACH_HEADERS = {"family": "method_approach",
                    "method_family": "method_approach",
                    "tool_family": "method_approach"}

def normalise_header(row):
    return [APPROACH_HEADERS.get(c, c) if isinstance(c, str) else c for c in row]

# Declared per-sheet transformations of the source TSV. These are the only ways a
# sheet may differ from its source file; the SI QC gate reads this table from the
# AST (never by import, since importing this module rewrites the workbook).
SHEET_TRANSFORMS = {
    # Reader-facing wording: the novel panel is split into homology groups.
    "S37": {"rename": {"stratum": "group", "stratum_label": "group_label"}},
    # Internal verification provenance (local literature-cache paths) is not part
    # of the published taxonomy; the remaining columns are unchanged.
    "S37_taxonomy": {"drop": ["in_docs_references_md", "verification_source"]},
}

def apply_transforms(name, rows):
    spec = SHEET_TRANSFORMS.get(name)
    if not spec or not rows:
        return rows
    # Columns are dropped when present: the public copy of the source TSV already
    # omits them, so the same builder yields the same sheet from either copy.
    drop = set(spec.get("drop", ()))
    keep = [i for i, c in enumerate(rows[0]) if c not in drop]
    rows = [[r[i] for i in keep if i < len(r)] for r in rows]
    ren = spec.get("rename", {})
    return [[ren.get(c, c) for c in rows[0]]] + rows[1:]

def add_sheet(wb, name, rows, header=True):
    rows = apply_transforms(name, rows)
    ws = wb.create_sheet(title=name[:31])
    if header and rows:
        rows = [normalise_header(rows[0])] + list(rows[1:])
    widths = {}
    for r, row in enumerate(rows, 1):
        for c, val in enumerate(row, 1):
            cell = ws.cell(row=r, column=c, value=conv(val) if isinstance(val, str) else val)
            cell.font = HEAD_FONT if (header and r == 1) else BODY_FONT
            if header and r == 1:
                cell.fill = HEAD_FILL
                cell.alignment = WRAP
            widths[c] = min(max(widths.get(c, 8), len(str(val)) + 2), 60)
    for c, w in widths.items():
        ws.column_dimensions[ws.cell(row=1, column=c).column_letter].width = w
    if header and rows:
        ws.freeze_panes = "A2"
    return ws

# ordered (sheet_name, source_file, description)
PLAN = [
    ("Table2_full", "table2_tool_recommendations_full.tsv", "Main Table 2 (tool recommendations), full: all 14 tools incl. runtimes, key strengths, caveats"),
    ("S1", "table_s1_track_a_by_category_L1500.tsv", "Track A per-category counts and per-tool recall at 1,500 bp"),
    ("S2", "table_s2_track_c_by_cst.tsv", "Track C performance by CST (with constituent VISTA mgCSTs)"),
    ("S3", "table_s3_track_c_by_virus_category.tsv", "Track C true-positive composition by virus category"),
    ("S4", "table_s4_bootstrap_confidence_intervals.tsv", "Bootstrap 95% CIs for every metric x length"),
    ("S5", "table_s5_cochrans_q_test.tsv", "Cochran's Q test per Track A length"),
    ("S6", "table_s6_mcnemar_significant_pairs.tsv", "Pairwise McNemar (BH-corrected) tool pairs"),
    ("S7", "table_s7_track_a_overall_metrics_by_length.tsv", "Track A per-tool metrics across six lengths"),
    ("S8", "table_s8_track_b_metrics.tsv", "Track B metrics (background x coverage x tool)"),
    ("S9", "table_s9b_track_b_multibackground_excl.tsv", "Track B four-vs-four multi-background CST-I/CST-IV-B contrast (MCC, precision and FPR differences with background-block bootstrap CIs; per-arm sample counts, precision/recall/FPR means, relative precision drop, post hoc leave-one-background-out precision range)"),
    ("S10", "table_s10_track_b_recall_by_coverage.tsv", "Track B recall by coverage (>=1,500 bp)"),
    ("S11", "table_s11_multi_evidence_sensitivity.tsv", "Multi-evidence tier robustness (3 checks)"),
    ("S12", "table_s12_secondary_benchmark_overall.tsv", "Per-tool metrics on the 13,030-contig multi-evidence benchmark"),
    ("S13", "table_s13_q2_prophage_detection.tsv", "Q2 prophage recall + free-phage FPR"),
    ("S14", "table_s14_dark_matter_per_contig.tsv", "Dark-matter per-contig evidence + BLASTx"),
    ("S15", "table_s15_dark_matter_detection_rates.tsv", "Per-tool dark-matter detection rate"),
    ("S16", "table_s16_resource_usage_L1500.tsv", "Runtime, peak RAM/VRAM"),
    ("S17", "table_s17_genomad_dark_matter_threshold_sweep.tsv", "geNomad virus-score sweep on 48 dark-matter contigs"),
    ("S18_MCC", "table_s18_track_c_R_sensitivity_mcc.tsv", "Track C MCC across 3x3 R grid"),
    ("S18_rank", "table_s18_track_c_R_sensitivity_ranking.tsv", "Track C per-cell ranking (former S19)"),
    ("S19_recall", "table_s19_rca_sensitivity_grid.tsv", "RCA-validated recall across 4x4 Evidence-5 mapping grid"),
    ("S19_rank", "table_s19_rca_ranking_stability.tsv", "Evidence-5 grid recall-based ranking (former S21)"),
    ("S20", "table_s20_track_a_by_source_genome_L1500.tsv", "Per-source-genome Track A recall/precision at L1500"),
    ("S21", "table_s21_tool_panel_metadata.tsv", "Tool versions, container digests, training DBs, thresholds, scope"),
    # S22 panels handled specially (after S21)
    ("S23", None, "GMM relabelling: per-tool MCC default vs GMM"),
    ("S23_metadata", None, "GMM fit parameters, BIC, set sizes, ranking comparison"),
    ("S24_per_vOTU", "table_s24_vmgc_benchmark_overlap.tsv", "Per-vOTU MetaVR/RefSeq overlap (4,263 vOTUs)"),
    ("S24_summary", "table_s24_vmgc_benchmark_overlap_summary.tsv", "VMGC overlap summary contingency"),
    ("S24_by_phylum", "table_s24_vmgc_overlap_by_host_phylum.tsv", "VMGC overlap by host phylum"),
    ("S24_by_genus", "table_s24_vmgc_overlap_by_host_genus.tsv", "VMGC overlap by host genus"),
    ("S24_by_family", "table_s24_vmgc_overlap_by_family.tsv", "VMGC overlap by viral family"),
    ("S25_dark_matter", "table_s25_dark_matter_vmgc.tsv", "VMGC cross-check of 48 dark-matter contigs"),
    ("S25_novel_panel", "table_s25_novel_panel_vmgc.tsv", "VMGC cross-check of 15 ANI-novel genomes"),
    ("S26_summary", "table_s26_e1b_parser_audit.tsv", "E1b parser audit per-sample counts"),
    ("S26_tier_changes", "table_s26_e1b_parser_audit_tier_changes.tsv", "E1b per-contig tier changes (large)"),
    ("S27", "table_s27_mitch_ranking_ci.tsv", "MITCH external validation: per-tool pooled MCC + sample-block CIs vs discovery"),
    ("S27_pairwise", "table_s27_pairwise_mcc.tsv", "Paired MCC differences among the five decision-relevant tools: primary assemblies and Track C (contig resampling), MiTCH (sample resampling); 2,000 draws"),
    ("S28", "table_s28_mitch_per_cst.tsv", "MITCH: top tools per traditional CST"),
    ("S29", "table_s29_mitch_per_mgcst.tsv", "MITCH: top tools per metagenomic CST (mgCST)"),
    ("S30", "table_s30_mitch_per_sample_mcc.tsv", "MITCH: per-sample MCC matrix (30 samples x 14 tools)"),
    ("S31", "table_s31_mitch_diversity_contrast.tsv", "MITCH: high- vs low-diversity delta-MCC + sample-block CIs"),
    ("S32", "table_s32_mitch_per_sample_cst.tsv", "MITCH: per-sample CST/mgCST assignments"),
    ("S33", "table_s33_spike_in_panel.tsv", "Spike-in reference panel: 14 viral + 6 bacterial genomes (former Main Table 1)"),
    ("S34", "table_s34_track_a_overall_metrics_full.tsv", "Track A per-tool metrics at 1,500 bp, all 14 tools, full confusion counts + metrics (former Main Table 2)"),
    ("S35", "table_s35_track_c_overall_metrics_full.tsv", "Track C per-tool full confusion counts + metrics, all 14 tools (former Main Table 3)"),
    ("S36", "table_s36_evidence_lines.tsv", "Evidence lines (E1a-E5) plus negative grounding: signals, thresholds, and references for the multi-evidence benchmark (former Main Table 5)"),
    ("S37", "table_s37_novel_panel_recall_by_stratum.tsv", "ANI-novel panel: per-tool recall by novelty group and fragment length (source data for Figure 4)"),
    # The Table S37 caption tells readers the homology_dependence column is
    # assigned from this file's verified_algorithm column. It therefore has to
    # ship with the supplement; before 2026-08-11 it did not, so the caption
    # pointed at a file no reader could obtain.
    ("S37_taxonomy", "tool_method_taxonomy_verified.tsv", "Per-tool verified method taxonomy: methodological approach, verified algorithm, verification signal, and primary reference (source of the homology_dependence assignment in Table S37)"),
    ("S38", "table_s38_negative_grounding_specificity.tsv", "Per-tool false-positive rate on the multi-evidence negatives split by grounding category (human, bacterial/archaeal, CheckV host genes), with a Track C human replication"),
    ("S39", "table_s39_mitch_contrast_estimands.tsv", "MITCH: CST-IV minus CST-I/III/V delta-MCC under pooled, mean per-sample and prevalence-standardised estimands, with sample-block CIs"),
    ("S40a", "table_s40a_mitch_per_sample_diversity.tsv", "MITCH per sample: bacterial Shannon diversity, richness, Lactobacillus fraction, labelled viral contigs by category, viral fraction and per-tool MCC"),
    ("S40b", "table_s40b_mitch_diversity_correlations.tsv", "MITCH: Spearman correlations of bacterial diversity with labelled viral contigs, viral fraction and per-sample MCC (all samples; CST-IV-B only)"),
    ("S41", "table_s41_threshold_calibration.tsv", "Track A threshold calibration at 1,500 and 3,000 bp: default rules versus optimistic in-sample MCC optima"),
]

wb = Workbook()
wb.remove(wb.active)
idx_rows = [["Sheet", "Contents"]]

s22 = parse_panels(HERE / "table_s22_cross_track_spearman.tsv")
s23_meta, s23_table = parse_meta_table(HERE / "table_s23_track_c_gmm_sensitivity.tsv")

def emit(name, desc):
    idx_rows.append([name, desc])

for name, src, desc in PLAN:
    if name == "S21":
        add_sheet(wb, name, read_tsv(HERE / src)); emit(name, desc)
        for L in ("A", "B", "C", "D"):
            p = s22[L]
            add_sheet(wb, f"S22_{L}", p["rows"]); emit(f"S22_{L}", p["title"])
        continue
    if name == "S23":
        add_sheet(wb, "S23", s23_table); emit("S23", desc); continue
    if name == "S23_metadata":
        add_sheet(wb, "S23_metadata", [["GMM relabelling metadata (Table S23)"]] + [[m] for m in s23_meta], header=True)
        emit("S23_metadata", desc); continue
    add_sheet(wb, name, read_tsv(HERE / src)); emit(name, desc)

index = wb.create_sheet("Index", 0)
for r, row in enumerate(idx_rows, 1):
    for c, val in enumerate(row, 1):
        cell = index.cell(row=r, column=c, value=val)
        cell.font = HEAD_FONT if r == 1 else BODY_FONT
        if r == 1:
            cell.fill = HEAD_FILL
index.column_dimensions["A"].width = 18
index.column_dimensions["B"].width = 75
index.freeze_panes = "A2"

wb.save(OUT)
print(f"WROTE {OUT}  ({len(wb.sheetnames)} sheets)")
print("Sheets:", ", ".join(wb.sheetnames))
