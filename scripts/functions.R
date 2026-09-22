#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   functions.R  - shared style and helpers for scripts/figures.R   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Sourced by scripts/figures.R after it has attached the packages and defined
# ROOT, TABLES_DIR and FIGURES_DIR. Contents:
#   1. Publication style: tool metadata, palettes, themes, type scale
#      (mirrors scripts/07_evaluation/pub_style.py constants).
#   2. Save helper used by every figure.
#   3. Data loaders and one panel builder that more than one figure uses.
# Nothing figure-specific lives here.


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   1. Publication style                                            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Tool metadata (single source of truth) ---
# TOOL_ORDER is a FIXED DISPLAY ORDER, inherited from the first figure commit.
# It is not a ranking: it reproduces no metric ordering (checked against MCC,
# AUPRC and recall at all six Track A lengths) and does not group by
# TOOL_FAMILY. Its only job is to keep the tool axis identical across figures
# so a reader can track one tool between them. Do not describe it in a legend
# or comment as ordered by performance.
TOOL_ORDER <- c(
  "geNomad", "ViraLM", "TransGINmer", "VirSorter2", "Jaeger",
  "HVSeeker", "PPR-Meta", "VIBRANT", "DeepVirFinder", "Seeker",
  "VirSorter", "MetaPhinder", "VirFinder", "Sourmash"
)

TOOL_FAMILY <- c(
  DeepVirFinder = "Deep learning",
  Seeker        = "Deep learning",
  `PPR-Meta`    = "Deep learning",
  HVSeeker      = "Deep learning",
  Jaeger        = "Deep learning",
  TransGINmer   = "Attention/transformer",
  ViraLM        = "Attention/transformer",
  VirSorter     = "Feature-based ML",
  VirSorter2    = "Feature-based ML",
  VIBRANT       = "Feature-based ML",
  VirFinder     = "Feature-based ML",
  geNomad       = "DNN + markers",
  MetaPhinder   = "Reference-based",
  Sourmash      = "Reference-based"
)

FAMILY_PALETTE <- c(
  "Deep learning"         = "#1f77b4",
  "Attention/transformer" = "#9467bd",
  "Feature-based ML"      = "#2ca02c",
  "DNN + markers"         = "#ff7f0e",
  "Reference-based"       = "#d62728"
)

# --- Marker/reference dependence vs sequence composition ---
#
# A second, coarser cut across the same 14 tools, orthogonal to TOOL_FAMILY: does
# classifying a contig require matching it against an external sequence resource
# at run time (whole-genome references, MinHash signatures, or protein/HMM marker
# profiles), or is the score computed from the contig's own sequence composition
# alone? Assignments are read off the verified per-tool algorithm column of
# tool_method_taxonomy_verified.tsv in TABLES_DIR (each row there traced to the
# primary tool paper), not from TOOL_FAMILY, which splits tools by architecture
# and therefore puts the k-mer classifier VirFinder next to the HMM-driven
# VIBRANT. Used by fig5, where this is the axis the result is about.
# The same assignment ships as the homology_dependence column of Table S37;
# qc_numeric_claims.py::D-fig5-hybrid-flag holds the two in agreement.
TOOL_APPROACH <- c(
  MetaPhinder   = "marker_ref",   # BLASTn/ANI vs whole-genome phage references
  Sourmash      = "marker_ref",   # MinHash containment vs a signature database
  VirSorter     = "marker_ref",   # viral protein/HMM hits plus genome context
  VirSorter2    = "marker_ref",   # one RF over HMM hits + gene-organisation stats
  VIBRANT       = "marker_ref",   # NN over KEGG/Pfam/VOG annotation profiles
  geNomad       = "marker_ref",   # IGLOO sequence branch plus marker-gene branch
  VirFinder     = "composition",  # logistic regression over k-mer frequencies
  DeepVirFinder = "composition",  # CNN over one-hot encoded nucleotide sequence
  `PPR-Meta`    = "composition",  # bi-path CNN over base and codon encodings
  Jaeger        = "composition",  # dilated CNN over six-frame translations
  HVSeeker      = "composition",  # bidirectional LSTM over nucleotide sequence
  Seeker        = "composition",  # LSTM over nucleotide sequence
  ViraLM        = "composition",  # DNABERT-2 transformer over raw sequence
  TransGINmer   = "composition"   # self-attention over k-mer/codon tokens + GIN
)

# geNomad is the only tool in the marker_ref group with a SEPARATE alignment-free
# classifier: its IGLOO branch scores the nucleotide sequence on its own, and that
# score is then aggregated with a marker-branch score (Camargo et al., Fig. 1, which
# calls the two "an alignment-free classifier (sequence branch)" and "a gene-based
# classifier (marker branch)"; verified against PMC11324519 on 2026-09-15).
#
# An earlier version of this comment added "It still cannot classify without the
# marker database." That is contradicted by the same source, which states the
# marker and sequence branches CAN be run independently, with good classification
# performance from the sequence branch alone. What is true of THIS benchmark is
# narrower and is what the captions now say: geNomad was run as `genomad
# end-to-end` (modules/local/virus_id/genomad/main.nf), so both branches were
# active. Do not restate the stronger claim.
#
# Flagged with a dagger by tool_approach_label() so the figure does not imply it is
# purely reference-based.
#
# The flag is deliberately NOT "also scores non-homology features". That criterion
# is met by three of the six, and flagging only some of them was the bug this
# replaced:
#   VirSorter2  ONE random forest per viral group over a single mixed vector; of 27
#               features (Guo et al. 2021, Table 1) only 7 are HMM-derived, the
#               other 20 gene-organisation and GC statistics. No separate branch.
#   VirSorter   2 of its 6 metrics (short-gene enrichment, strand-switch depletion)
#               are genome structure; the abstract claims both a "reference-dependent
#               and reference-independent manner". No separate branch.
#   VIBRANT, MetaPhinder, Sourmash are genuinely homology-only.
# Separate-branch architecture is the sharper cut and the one Fig. 5 is about: it
# picks out the tool you would predict should cope with novel sequence.
#
# Avoid the word "composition" for VirSorter2: Guo et al. use it for the k-mer camp
# (VirFinder, DeepVirFinder) and place VirSorter2 in the "gene content and genomic
# structural features" camp opposite it.
TOOL_DUAL_BRANCH <- c("geNomad")

APPROACH_LEVELS <- c("marker_ref", "composition")
APPROACH_LABELS <- c(
  marker_ref  = "Marker/reference-dependent",
  composition = "Composition-based (alignment-free)"
)

# Append a dagger to the dual-branch tool(s); used for axis labels wherever tools
# are shown split by TOOL_APPROACH.
tool_approach_label <- function(tools) {
  tools <- as.character(tools)
  ifelse(tools %in% TOOL_DUAL_BRANCH, paste0(tools, "\u2020"), tools)
}

# "Marker/reference-dependent\n(n = 6)" etc., for facet strips. `tools` is the
# tool subset actually plotted, so the counts cannot drift from the panel. The
# count goes on its own line because the marker/reference strip sits over only
# 6 of 14 columns, and the one-line form overran that panel and was clipped.
approach_strip_labels <- function(tools) {
  n <- table(factor(TOOL_APPROACH[as.character(tools)],
                    levels = APPROACH_LEVELS))
  setNames(sprintf("%s\n(n = %d)", APPROACH_LABELS[APPROACH_LEVELS],
                   as.integer(n[APPROACH_LEVELS])),
           APPROACH_LEVELS)
}

# Per-tool colours (tab20 equivalents)
.tab20 <- c(
  "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c",
  "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
  "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f",
  "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"
)
TOOL_COLORS <- setNames(.tab20[seq_along(TOOL_ORDER)], TOOL_ORDER)

# --- Declared target scope -------------------------------------------------
#
# What virus types does the tool CLAIM to detect? A third cut across the same 14
# tools, orthogonal to both TOOL_FAMILY (architecture) and TOOL_APPROACH
# (homology dependence). This is the axis Fig. 3B is about: a tool that targets
# prokaryotic viruses only and scores 0.00 on herpesvirus is out of scope, not
# failing, and the figure must not let those two read alike.
#
# Values are read off the scope_label column of
# table_s21_tool_panel_metadata.tsv in TABLES_DIR, which is the verified,
# published source (Table S21). Do not edit them here from recollection: the
# code and the table are compared cell-by-cell by scripts/qc_tool_scope.py, and
# they disagreed for three tools before 2026-09-15 --
#   HVSeeker      was "all_virus", Table S21 and manuscript both say phage-only
#   VIBRANT       was "all_virus", Table S21 says prokaryotic viruses
#   DeepVirFinder was "all_virus", Table S21 says prokaryotic-virus-focused
# -- which mis-coloured the tool labels in Fig. S7 (fig_s7 in figures.R) and would
# have put all three in the wrong block here.
TOOL_SCOPE <- c(
  DeepVirFinder = "prokaryotic_virus",   # "prokaryotic-virus-focused"
  VirSorter     = "phage_only",
  VirSorter2    = "all_virus",
  VirFinder     = "phage_only",
  `PPR-Meta`    = "phage_and_plasmid",
  VIBRANT       = "prokaryotic_virus",   # "prokaryotic viruses"
  Seeker        = "phage_only",
  MetaPhinder   = "phage_only",
  Sourmash      = "db_dependent",        # "user-provided signature database"
  geNomad       = "all_virus",
  HVSeeker      = "phage_only",
  Jaeger        = "phage_only",
  TransGINmer   = "all_virus",
  ViraLM        = "all_virus"
)

# Coarse grouping used for the FACET BLOCKS of Fig. 3B. Fig. 3B answers Q2 --
# "Does a tool's declared target scope predict eukaryote-infecting virus
# detection?" -- so it is grouped by SCOPE, not by TOOL_APPROACH (Fig. 4 / Q3)
# and not by TOOL_FAMILY. Grouping it by either of those would show the panel
# against the axis of a different question than its own subsection asks.
#
# The six scope_label values in Table S21 are too fine to facet on (four hold a
# single tool), so they collapse to three blocks on one question: can the tool be
# expected to call a eukaryote-infecting virus at all?
#   all_virus    claims eukaryote-infecting viruses too            (4 tools)
#   prokaryotic  phage-only + prokaryotic-virus + phage+plasmid    (9 tools)
#   db_dependent scope is whatever signature database was supplied (1: Sourmash)
# Sourmash keeps its own block rather than being folded into either: its Table
# S21 database is "user-provided signature database", so its scope is a property
# of the run, not of the tool.
#
# NOTE the blocks are expected to OVERLAP on the eukaryote-infecting rows. That
# overlap is the Q2 result ("Published scope alone does not describe
# eukaryote-infecting virus recall"), not a defect in the stratification. Do not
# "fix" it by regrouping.
TARGET_LEVELS <- c("all_virus", "prokaryotic", "db_dependent")
TARGET_LABELS <- c(
  all_virus    = "All-virus targets",
  prokaryotic  = "Prokaryotic-virus targets",
  # Wrapped short because this block is ONE tile wide (Sourmash only).
  db_dependent = "DB-\ndep."
)
TOOL_TARGET <- setNames(
  ifelse(TOOL_SCOPE == "all_virus",    "all_virus",
  ifelse(TOOL_SCOPE == "db_dependent", "db_dependent", "prokaryotic")),
  names(TOOL_SCOPE)
)

# Strip labels carrying each block's n, counted from the tools actually plotted
# so the counts cannot drift from the panel (same contract as
# approach_strip_labels()). No parentheses around the n: "(n = 1)" is two
# characters wider than "n = 1" and overran the one-tile DB-dep. strip.
target_strip_labels <- function(tools) {
  n <- table(factor(TOOL_TARGET[as.character(tools)], levels = TARGET_LEVELS))
  setNames(sprintf("%s\nn = %d", TARGET_LABELS[TARGET_LEVELS],
                   as.integer(n[TARGET_LEVELS])),
           TARGET_LEVELS)
}

# One entry per LEVEL OF TOOL_SCOPE -- keep the two in step. Used by
# Fig. S7 (fig_s7 in figures.R) to colour tool labels; a missing level silently
# yields NA and paints that label black, so adding a scope value here is not
# optional. scripts/qc_tool_scope.py checks the pairing.
SCOPE_LABEL_COLORS <- c(
  phage_only        = "#e74c3c",
  prokaryotic_virus = "#e67e22",
  all_virus         = "#2980b9",
  phage_and_plasmid = "#8e44ad",
  db_dependent      = "#7f8c8d"
)

# --- Length bins ---
ALL_LENGTHS <- c("L500", "L1000", "L1500", "L3000", "L5000", "L10000")

LENGTH_COLORS <- c(
  L500   = "#c6dbef",
  L1000  = "#9ecae1",
  L1500  = "#6baed6",
  L3000  = "#3182bd",
  L5000  = "#08519c",
  L10000 = "#08306b"
)

# --- Category metadata ---
PHAGE_CATS <- c("lactobacillus_phage", "gardnerella_phage",
                "megasphaera_phage", "fannyhessea_phage", "sneathia_phage")
EUKARYOTIC_CATS <- c("hpv", "herpesvirus", "anellovirus")

# The five PHAGE_CATS are phage categories named by HOST genus, so the bare
# genus name ("Lactobacillus") reads on an axis as recall on the bacterium
# itself rather than on its phages -- exactly the wrong claim for a viral-
# detection benchmark. Carry " phage" in the display string so no axis or
# caption has to supply it. The three EUKARYOTIC_CATS are virus names already
# and are left as they are.
CAT_DISPLAY <- c(
  lactobacillus_phage = "Lactobacillus phage",
  gardnerella_phage   = "Gardnerella phage",
  megasphaera_phage   = "Megasphaera phage",
  fannyhessea_phage   = "Fannyhessea phage",
  sneathia_phage      = "Sneathia phage",
  hpv                 = "HPV",
  herpesvirus         = "Herpesvirus",
  anellovirus         = "Anellovirus"
)

# --- CST metadata ---
SAMPLE_CST <- c(
  UC028_V2 = "CST-IV", UC055_V1 = "CST-III", UC055_V2 = "CST-III",
  UC062_V2 = "CST-I",  UC065_V2 = "CST-IV",  UC074_V2 = "CST-III",
  UC084_V2 = "CST-I",  UC093_V2 = "CST-IV",  UC093_V3 = "CST-IV",
  UC096_V2 = "CST-I",  UC115_V2 = "CST-I",   UC139_V2 = "CST-I",
  UC164_V2 = "CST-III"
)

CST_COLORS <- c("CST-I" = "#3182bd", "CST-III" = "#31a354",
                "CST-IV" = "#e6550d", "CST-IV-B" = "#e6550d")

CST_LABELS <- c(UC115_V2 = "CST-I", UC093_V3 = "CST-IV-B")

TIER_COLORS <- c("1" = "#d62728", "2" = "#ff7f0e", "3" = "#fdd49e", "0" = "#bdbdbd")

# --- Publication theme ---
theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      text             = element_text(family = "sans"),
      axis.title       = element_text(size = 12),
      axis.text        = element_text(size = 10),
      plot.title       = element_text(size = 12, face = "bold", hjust = 0),
      legend.title     = element_text(size = 10),
      legend.text      = element_text(size = 9),
      legend.background = element_blank(),
      # Facet strips as filled label chips rather than theme_classic's hard
      # black-outlined white box. theme_classic draws strip.background with a
      # rel(2) black border, which at print size reads as a second panel frame
      # stacked on the axis line and competes with the data. A soft fill with no
      # border keeps the strip a label. Applied here, not per figure, so every
      # faceted panel in the pack matches: fig3A (MCC/AUPRC), fig4A
      # (AUPRC/MCC/recall), fig5 (the two tool blocks), fig7C, fig9B, fig_s6.
      # NB %+replace% swaps the whole element, so margin must be given here or
      # the strip text sits flush against the chip edge.
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text       = element_text(size = 10, face = "bold",
                                      margin = margin(3.5, 3.5, 3.5, 3.5)),
      panel.grid       = element_blank()
    )
}

# --- Print-size type scale (>= 8 pt at the final journal width) ---
#
# A figure's *printed* font size is the declared size times
# (print width / rendered canvas width). The figure driver historically rendered on
# oversized canvases (e.g. fig4 at 14.4 in) and let the journal shrink them, so
# a declared 10 pt printed at ~4.9 pt in a 180 mm column. Figures that set
# `scale = 1` and a canvas equal to the print width instead get
# declared pt == printed pt, and can use the constants below directly.
#
# PT_MIN is the floor: no text in such a figure may be declared smaller.
PT_MIN        <- 8    # axis text, legend text, in-tile labels
PT_AXIS_TITLE <- 9
PT_STRIP      <- 9
PT_TAG        <- 10

# geom_text()/geom_label() take `size` in mm, not points; ggplot2::.pt converts.
GEOM_TEXT_MIN <- PT_MIN / .pt   # ~2.81 mm == 8 pt

# Applied with patchwork's `&` after the panels are assembled, so it overrides
# whatever theme_pub() set. Element-specific settings inside a builder (e.g.
# axis.text.x angle/size) survive, because they are distinct theme elements that
# only inherit *unset* properties from their parent.
theme_print <- function() {
  theme(
    axis.title   = element_text(size = PT_AXIS_TITLE),
    axis.text    = element_text(size = PT_MIN),
    strip.text   = element_text(size = PT_STRIP, face = "bold"),
    legend.title = element_text(size = PT_MIN),
    legend.text  = element_text(size = PT_MIN),
    plot.tag     = element_text(size = PT_TAG, face = "bold")
  )
}

# Colour y-axis tick labels by methodological approach
colour_axis_by_family <- function(p, axis = "y") {
  g <- ggplot_build(p)
  if (axis == "y") {
    labs <- ggplot_build(p)$layout$panel_params[[1]]$y$get_labels()
  } else {
    labs <- ggplot_build(p)$layout$panel_params[[1]]$x$get_labels()
  }
  cols <- FAMILY_PALETTE[TOOL_FAMILY[labs]]
  cols[is.na(cols)] <- "black"
  if (axis == "y") {
    p + theme(axis.text.y = element_text(colour = cols))
  } else {
    p + theme(axis.text.x = element_text(colour = cols))
  }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   2. Save helper                                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Global element scaling for saved figures. scale < 1 enlarges all text, data
# points, lines, and annotations relative to the canvas (lower = larger), while
# keeping each figure's output dimensions unchanged. Tune for legibility.
ELEMENT_SCALE <- 0.8

# Every figure is written as a 600 dpi PNG with a white background. ggsave
# renders on a (w x scale) by (h x scale) inch device, so a figure can opt out
# of the global ELEMENT_SCALE by passing its own `scale` and sizing its canvas
# to the intended print width instead. Returns the output path.
save_fig <- function(p, out_path, w, h, scale = ELEMENT_SCALE) {
  ggsave(out_path, p, width = w, height = h, dpi = 600,
         scale = scale, limitsize = FALSE, bg = "white")
  out_path
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   3. Shared data loaders and builders                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Track A (fragment benchmark) -------------------------------------------

# Table S7 per-length overall metrics. Used by fig2 and fig_s2.
.load_fig2_trackA_overall <- function(root) {
  # Table S7 contains the current estimates underlying Tables S4 and S34.
  # Historical per-run benchmark_results.tsv files contain older AUPRC values.
  # Use a figure-specific name: the Figure 3 module defines its own loader.
  read_tsv(file.path(TABLES_DIR,
                     "table_s7_track_a_overall_metrics_by_length.tsv"),
           show_col_types = FALSE) %>%
    mutate(tool = factor(tool, levels = TOOL_ORDER),
           length_bin = factor(length_bin, levels = ALL_LENGTHS))
}

# Per-length benchmark_results.tsv, stacked. Used by fig3, fig4 and fig_s17.
.load_trackA_overall <- function(root) {
  frames <- list()
  for (len in ALL_LENGTHS) {
    f <- file.path(root, "results/spike_in_benchmark",
                   paste0(len, "_fragments"), "benchmark_results.tsv")
    if (!file.exists(f)) next
    frames[[len]] <- read_tsv(f, show_col_types = FALSE) %>%
      mutate(length_bin = len)
  }
  bind_rows(frames) %>%
    mutate(tool = factor(tool, levels = TOOL_ORDER),
           length_bp = as.integer(sub("L", "", length_bin)),
           length_bin = factor(length_bin, levels = ALL_LENGTHS))
}

# Per-length benchmark_results_by_category.tsv, stacked. Used by fig4, fig_s17.
.load_trackA_bycategory <- function(root) {
  frames <- list()
  for (len in ALL_LENGTHS) {
    f <- file.path(root, "results/spike_in_benchmark",
                   paste0(len, "_fragments"), "benchmark_results_by_category.tsv")
    if (!file.exists(f)) next
    frames[[len]] <- read_tsv(f, show_col_types = FALSE) %>%
      mutate(length_bin = len)
  }
  bind_rows(frames) %>%
    mutate(tool = factor(tool, levels = TOOL_ORDER),
           length_bin = factor(length_bin, levels = ALL_LENGTHS))
}

# --- Track B (assembly benchmark) -------------------------------------------

# Spike-in co-assembly ground truth, both 10x backgrounds. Used by fig6, fig_s18.
.load_trackB_assemblies <- function(root) {
  assembly_dir <- file.path(root, "data/spike_in/assemblies")
  coverages <- c("0.1", "0.5", "1", "5", "10", "50")
  backgrounds <- c("UC115_V2", "UC093_V3")
  gt_list <- list()
  for (bg in backgrounds) {
    for (cov in coverages) {
      gt_file <- file.path(assembly_dir, bg,
                           paste0("assembly_cov", cov), "ground_truth.tsv")
      if (!file.exists(gt_file)) next
      gt_list[[paste(bg, cov)]] <- read_tsv(gt_file, show_col_types = FALSE) %>%
        mutate(background = bg, coverage = as.numeric(cov),
               CST = CST_LABELS[bg])
    }
  }
  bind_rows(gt_list)
}

# --- Multi-background CST contrast (membership from the active manifest) ---
# Replaces the original single-pair panel. The main figure and Fig. S19 show
# the nine highest CST-I mean MCC tools in the current Table S9; Fig. S12 shows
# all 14. Select again after rescoring instead of retaining an obsolete list.

.ARM_COLORS <- c("CST-I" = unname(CST_COLORS["CST-I"]),
                 "CST-IV-B" = unname(CST_COLORS["CST-IV-B"]))

.load_trackB_manifest <- function(root) {
  f <- file.path(root, "data/expansion_cohorts/track_b_backgrounds.tsv")
  if (!file.exists(f)) stop("Missing track_b_backgrounds.tsv")
  manifest <- read_tsv(f, show_col_types = FALSE) %>%
    transmute(background = run, stratum)
  if (anyDuplicated(manifest$background) ||
      !setequal(unique(manifest$stratum), names(.ARM_COLORS))) {
    stop("Track B manifest must contain unique backgrounds in both CST arms")
  }
  manifest
}

.trackB_sample_sizes <- function(root) {
  manifest <- .load_trackB_manifest(root)
  sprintf("CST-I: n = %d; CST-IV-B: n = %d backgrounds",
          sum(manifest$stratum == "CST-I"),
          sum(manifest$stratum == "CST-IV-B"))
}

.load_trackB_deltas <- function(root) {
  f <- file.path(TABLES_DIR, "table_s9b_track_b_multibackground_excl.tsv")
  if (!file.exists(f)) stop("Missing table_s9b_track_b_multibackground_excl.tsv")
  manifest <- .load_trackB_manifest(root)
  d <- read_tsv(f, show_col_types = FALSE)
  for (arm in names(.ARM_COLORS)) {
    n_column <- if (arm == "CST-I") "n_CST_I" else "n_CST_IVB"
    if (n_column %in% names(d) &&
        any(is.na(d[[n_column]]) | d[[n_column]] != sum(manifest$stratum == arm))) {
      stop("Track B delta table is stale relative to the active background manifest")
    }
  }
  d %>%
    mutate(resolved = ifelse(as.character(delta_excludes_0) %in% c("TRUE", "True"),
                             "resolved", "not resolved"))
}

.trackB_main_tools <- function(root) {
  .load_trackB_deltas(root) %>%
    arrange(desc(MCC_CST_I_mean), tool) %>%
    slice_head(n = 9) %>%
    pull(tool)
}

# Fig. 6B (and, with all 14 tools, Fig. S12): per-tool delta (CST-IV-B minus
# CST-I) with bootstrap 95% CI.
build_fig5d_delta_ci <- function(root = ROOT,
                                 tools = .trackB_main_tools(root)) {
  d <- .load_trackB_deltas(root) %>%
    filter(tool %in% tools) %>%
    mutate(tool = factor(tool, levels = rev(tools)))
  ggplot(d, aes(x = delta_IVB_minus_I, y = tool)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = delta_CIl, xmax = delta_CIh), height = 0.25,
                   linewidth = 0.5) +
    geom_point(aes(shape = resolved), size = 2.6) +
    # Display only: the internal levels stay "resolved"/"not resolved" (set from
    # delta_excludes_0 in .load_trackB_deltas), but "resolved" does not say what
    # was resolved. Spell out the actual criterion on the key.
    scale_shape_manual(values = c("resolved" = 19, "not resolved" = 1),
                       labels = c("resolved"     = "CI excludes 0",
                                  "not resolved" = "CI includes 0"),
                       name = NULL) +
    labs(x = "ΔMCC (CST-IV-B − CST-I)", y = NULL,
         subtitle = .trackB_sample_sizes(root)) +
    theme_pub()
}

# --- Track C / multi-evidence ground truth ----------------------------------

# Used by fig7, fig_s20 and fig_s21.
.load_gt <- function(root) {
  # Single source of truth: the LOCKED 13-sample pooled benchmark-ready ground
  # truth (>= 1500 bp, tier 0/1/2; Gate 0), shared with the pooled multi-evidence
  # benchmark so the figure denominators match Table 3 / Table S12 exactly.
  gt_file <- file.path(root,
    "results/test_real/secondary_benchmark_pooled13/pooled_benchmark_ready_ground_truth.tsv")
  read_tsv(gt_file, show_col_types = FALSE) %>%
    mutate(CST = SAMPLE_CST[sample])
}

# Used by fig7 and fig_s21.
.load_rca <- function(root) {
  read_tsv(file.path(root, "results/test_real/coassembly/rca_ground_truth.tsv"),
           show_col_types = FALSE)
}

# --- MiTCH external validation ----------------------------------------------

# Used by fig8, fig_s13 and fig_s22.
.scaleup_dir <- function(root) {
  file.path(root, "results/test_real/viralm_investigation/mitch_scaleup30")
}

.load_scaleup <- function(root, name) {
  f <- file.path(.scaleup_dir(root), name)
  if (!file.exists(f)) stop("Missing ", name)
  read_tsv(f, show_col_types = FALSE)
}
