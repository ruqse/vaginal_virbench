#!/usr/bin/env Rscript
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   figures.R  - all R-drawn manuscript figures                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Writes one PNG per figure (<id>.png) into FIGURES_DIR, or into --out-dir.
# Every figure below has its own boxed header and its own function, so the
# objects it builds cannot leak into another figure. Shared style, palettes,
# loaders and the save helper live in scripts/functions.R.
#
# Usage (from any directory):
#   Rscript scripts/figures.R                          # every figure -> FIGURES_DIR
#   Rscript scripts/figures.R --only fig_s1            # one figure
#   Rscript scripts/figures.R --only fig2,fig_s9 --out-dir /tmp/figs
#
# R environment: envs/r-figures.yml. Figures were produced with R 4.3.3 and
# the package versions listed there; other versions can shift pixels.
#
# Figure IDs (= output file names), in render order:
#   fig2  fig3  fig4  fig5  fig6  fig7  fig8  fig9
#   fig_s1  fig_s2  fig_s3  fig_s4  fig_s5  fig_s6  fig_s7  fig_s8  fig_s9
#   fig_s10  fig_s12  fig_s13  fig_s17  fig_s18  fig_s19  fig_s20  fig_s21
#   fig_s22
# Not drawn here: fig1 (study design; results/figures/fig1.png is not made by
# any script in this repo, and the earlier Graphviz source is archived in
# results/figures/archive/fig1_graphviz_20260922/), fig_s14 (authored
# schematic), and fig_s11, fig_s15, fig_s16, which are drawn in Python (see
# the box before Fig. S17).
#
# Builder function names (build_fig2a_*, build_fig5a_*, ...) keep their legacy
# numbering and no longer track the display number. The figure ID is
# authoritative because it becomes the output file name.


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Load required packages, paths and custom functions              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(patchwork)
  library(optparse)
  library(ggplotify)
  library(pheatmap)
  library(RColorBrewer)
  library(reshape2)
  library(tibble)
})

# Repository root = parent of the directory that holds this file. Rscript
# passes the script as --file=<path>; `R -f <path>` passes it after -f. When
# neither is present (e.g. source() from an interactive session) fall back to
# $VBENCH_ROOT, then to the working directory if it contains scripts/figures.R.
ROOT <- local({
  args <- commandArgs(trailingOnly = FALSE)
  pre  <- args[seq_len(match("--args", args, nomatch = length(args) + 1L) - 1L)]
  self <- sub("^--file=", "", grep("^--file=", pre, value = TRUE))
  if (length(self) == 0) {
    i <- which(pre == "-f")
    if (length(i) > 0 && i[1] < length(pre)) self <- pre[i[1] + 1L]
  }
  self <- self[basename(self) == "figures.R"]
  if (length(self) > 0) {
    dirname(dirname(normalizePath(self[1], mustWork = TRUE)))
  } else if (nzchar(Sys.getenv("VBENCH_ROOT"))) {
    normalizePath(Sys.getenv("VBENCH_ROOT"), mustWork = TRUE)
  } else if (file.exists(file.path("scripts", "figures.R"))) {
    normalizePath(getwd(), mustWork = TRUE)
  } else {
    stop("Cannot find the repository root: run `Rscript scripts/figures.R`, ",
         "or set VBENCH_ROOT")
  }
})

TABLES_DIR  <- file.path(ROOT, "results", "tables")      # published tables (TSV)
FIGURES_DIR <- file.path(ROOT, "results", "figures")     # default output dir

source(file.path(ROOT, "scripts", "functions.R"))

# Command line: --out-dir DIR (default FIGURES_DIR), --only ID[,ID...]
opt <- parse_args(OptionParser(
  usage = "Rscript scripts/figures.R [--out-dir DIR] [--only ID[,ID...]]",
  option_list = list(
    make_option(c("-o", "--out-dir"), type = "character", default = FIGURES_DIR,
                help = "Output directory for the <id>.png files [default: %default]"),
    make_option(c("-f", "--only"), type = "character", default = NULL,
                help = "Comma-separated figure IDs to render (e.g. fig3,fig_s9)")
  )
))
OUT_DIR <- normalizePath(opt$`out-dir`, mustWork = FALSE)
ONLY    <- if (is.null(opt$only)) NULL else trimws(strsplit(opt$only, ",")[[1]])



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. 2A-B  - Overall classification performance (Track A)       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A: MCC + AUPRC grouped bar chart across all lengths.
# B: available-score precision-recall curves for the top 6 tools by Table S7
#    AUPRC, with implemented decision rules marked separately. Curves end where
#    usable scores end; they are not a graphical integral of full-cohort AUPRC.
# Rendered at print size (7.1 in, scale = 1), like fig3/fig4/fig5/fig9, so the
# PT_MIN = 8 floor holds.

fig2 <- function(out_dir) {

  # The fill keys are the internal length-bin codes (L500 ..); the legend should
  # show the fragment length the reader is being asked to think about. Keyed by
  # the palette names so the mapping cannot drift out of order.
  LENGTH_BIN_LABELS <- c(L500   = "500 bp",   L1000  = "1,000 bp",
                         L1500  = "1,500 bp", L3000  = "3,000 bp",
                         L5000  = "5,000 bp", L10000 = "10,000 bp")

  # -- Panel A --
  build_fig2a_mcc_auprc_ranking <- function(root = ROOT) {
    df <- .load_fig2_trackA_overall(root)
    df_long <- df %>%
      select(tool, length_bin, note, MCC, AUPRC) %>%
      pivot_longer(cols = c(MCC, AUPRC), names_to = "metric", values_to = "value") %>%
      # MCC is the primary ranking metric (Methods), so it occupies the top facet.
      mutate(metric = factor(metric, levels = c("MCC", "AUPRC")),
             # A tool flagged `no_output` parsed to zero predictions at this length
             # (Jaeger, which requires fragments >= 2,048 bp, at L500/L1000/L1500;
             # VIBRANT at L500). For that case benchmark_fragments.py takes its
             # `n_parsed == 0` branch and hard-codes MCC = 0 / AUPRC = "NA"; the
             # MCC 0 is a placeholder, not a measurement, yet it draws as a bar
             # indistinguishable from Sourmash's zeros, which ARE measured
             # (Sourmash parsed output and called nothing viral at its threshold,
             # note = "thr=0.5"). Blank the placeholder cells and label the slot
             # "n/a" instead -- the same convention as the recall heatmaps in
             # fig4 and fig5.
             value  = ifelse(!is.na(note) & note == "no_output", NA_real_, value))

    # Position the "n/a" marks by hand rather than with position_dodge(). One
    # dodged slot here is 0.8/6 category units = ~1.1 mm at the 180 mm print
    # width, while "n/a" at the PT_MIN = 8 pt floor is ~2.9 mm; three per-slot
    # labels for Jaeger (L500/L1000/L1500) therefore overprint each other. One
    # mark per tool per facet, centred on that tool's masked slots, is the most
    # that is legible: for Jaeger it spans exactly the three blanked slots, and
    # every tool with a mark has ALL of its blanked lengths under it.
    dodge_w <- 0.8
    n_bins  <- nlevels(df_long$length_bin)
    tool_lv <- levels(droplevels(df_long$tool))
    na_lab <- df_long %>%
      filter(!is.na(note), note == "no_output") %>%
      mutate(x = match(as.character(tool), tool_lv) - dodge_w / 2 +
                 (as.integer(length_bin) - 0.5) * dodge_w / n_bins) %>%
      group_by(tool, metric) %>%
      summarise(x = mean(x), .groups = "drop")

    ggplot(df_long, aes(x = tool, y = value, fill = length_bin)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7,
               colour = "black", linewidth = 0.2, na.rm = TRUE) +
      geom_text(data = na_lab, aes(x = x, y = -0.04, label = "n/a"),
                inherit.aes = FALSE, vjust = 0.5,
                size = GEOM_TEXT_MIN, colour = "grey35") +
      facet_wrap(~metric, scales = "free_y", ncol = 1) +
      scale_fill_manual(values = LENGTH_COLORS, labels = LENGTH_BIN_LABELS,
                        name = "Fragment\nlength") +
      labs(x = NULL, y = "Score") +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
  }

  # -- Panel B --
  # Curves cover only finite, reported scores; missing predictions never acquire
  # an invented trajectory. Recall always uses ALL ground-truth positives.
  # Dots show the implemented decision rules from Table S7. DeepVirFinder's rule
  # also filters p-values, so its dot need not lie on the score-only curve.
  # These partial curves are not a graphical integral of the full-cohort AUPRC.
  build_fig2b_pr_curves <- function(root = ROOT, len = "L1500",
                                    n_tools = 6) {
    score_file <- file.path(root, "results/spike_in_benchmark",
                            "pr_curve_scores.tsv.gz")
    if (!file.exists(score_file)) {
      stop("Missing ", score_file, "\n  Regenerate with:\n",
           "  python3 scripts/07_evaluation/dump_pr_curve_scores.py \\\n",
           "    --manifest data/spike_in/fragments/fragment_manifest.tsv \\\n",
           "    --results-root results/spike_in_benchmark --lengths ", len, " \\\n",
           "    --output ", score_file)
    }
    d <- read_tsv(score_file, show_col_types = FALSE) %>%
      filter(length_bin == len)
    if (nrow(d) == 0) stop("No dumped scores for ", len)

    ov <- .load_fig2_trackA_overall(root) %>% filter(length_bin == len)

    # Rank tools by the published AUPRC; curves for the top n_tools only. All 14
    # AUPRC values remain visible as bars in panel A, so nothing is lost.
    ranked <- ov %>%
      filter(!is.na(AUPRC)) %>%
      arrange(desc(AUPRC)) %>%
      mutate(tool = as.character(tool))
    keep <- head(ranked$tool, n_tools)

    # Tie-safe curve: evaluate ONLY at distinct score thresholds (one point per
    # tied block), as sklearn's precision_recall_curve does. Evaluating per row
    # instead makes the curve depend on how ties happen to be ordered, which for
    # the pre-filtering tools (geNomad, VirSorter2, VIBRANT) invents a flat
    # high-precision run out to recall 1.0 that is an artefact of fragment-ID sort
    # order, not tool behaviour.
    curves <- d %>%
      filter(tool %in% keep) %>%
      group_by(tool) %>%
      mutate(n_pos = sum(label)) %>%
      # A reported zero is a valid score. Conversely, VirSorter2 has reported
      # rows with NaN scores; those cannot support a score-threshold curve.
      filter(reported == 1, is.finite(score)) %>%
      arrange(desc(score), .by_group = TRUE) %>%
      mutate(tp = cumsum(label), k = row_number()) %>%
      group_by(tool, score) %>%
      slice_max(k, n = 1, with_ties = FALSE) %>%     # last row of each tie block
      ungroup() %>%
      mutate(precision = tp / k, recall = tp / n_pos) %>%
      arrange(tool, k) %>%
      select(tool, recall, precision)

    # Conventional empty-prediction starting point.
    curves <- bind_rows(
      tibble(tool = keep, recall = 0, precision = 1),
      curves
    ) %>%
      mutate(tool = factor(tool, levels = keep))

    # Precision/recall from the implemented decision rules, including any
    # additional filters, rather than a shared numerical score threshold.
    ops <- ov %>%
      filter(as.character(tool) %in% keep) %>%
      transmute(tool = factor(as.character(tool), levels = keep),
                recall, precision)

    prevalence <- ov %>%
      slice(1) %>%
      transmute(p = (TP + FN) / (TP + FN + TN + FP)) %>%
      pull(p)

    tool_cols <- TOOL_COLORS[keep]
    ggplot(curves, aes(x = recall, y = precision, colour = tool)) +
      geom_hline(yintercept = prevalence, linetype = "dotted", colour = "grey55") +
      geom_path(linewidth = 0.6) +
      geom_point(data = ops, size = 2.2, stroke = 0) +
      scale_colour_manual(values = tool_cols, name = "Tool") +
      guides(colour = guide_legend(override.aes = list(shape = NA))) +
      scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.01))) +
      scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
      labs(x = "Recall", y = "Precision") +
      theme_pub() +
      theme(legend.position = "right",
            legend.key.width = unit(0.6, "cm"))
  }

  # -- Assemble and save --
  p <- build_fig2a_mcc_auprc_ranking(ROOT) /
       build_fig2b_pr_curves(ROOT, "L1500")
  p <- p + plot_layout(heights = c(1.7, 1)) +
    plot_annotation(tag_levels = "A") &
    theme_print()
  save_fig(p, file.path(out_dir, "fig2.png"), w = 7.1, h = 7.6, scale = 1)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. 3  - Metrics vs fragment length (Track A, Q1)              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MCC / AUPRC / recall across the six fragment lengths, attention-based vs
# feature-based tools.
#
# fig3, fig4 and fig5 are rendered AT print size (7.1 in = 180 mm double-column)
# with scale = 1, so every declared point size is the printed point size and
# the PT_* floor in functions.R actually holds. The previous combined figure
# used an 18 x 12 @ 0.8 canvas = 14.4 in wide; shrunk into a 180 mm column
# that is x0.49, which printed the 10 pt axis text at ~4.9 pt.
#
# These were one four-panel figure carrying three separate claims, cited by the
# Results out of order. They were split on the question boundary and then
# ordered by DATASET: fig2 and fig3 both use the Track A reference panel (fig3
# extends fig2 from one fragment length to six), and fig5 introduces the
# separate ANI-novel panel.
#
# fig3 is a SINGLE panel (Q1 only). It carried the category heatmap as panel B
# until 2026-09-15, which made one figure straddle two questions and two
# different tool groupings -- colour-by-approach in A, blocks in B. The heatmap
# is now fig4 in its own right. No plot_annotation() here: a one-panel figure
# takes no "A" tag.

fig3 <- function(out_dir) {

  build_fig3a_metrics_vs_length <- function(root = ROOT) {
    # Fig. 3A visualises the H3 test from the manuscript subsection
    # "Attention-based tools outperform feature-based tools on short contigs".
    # We therefore restrict the panel to the two families named in that
    # subsection and colour curves by family (per-tool linetype keeps
    # individual tools resolvable — e.g. VirSorter2's L1500->L3000 jump).
    fig4a_families <- c("Attention/transformer", "Feature-based ML")
    keep_tools <- names(TOOL_FAMILY)[TOOL_FAMILY %in% fig4a_families]
    # Preserve TOOL_ORDER for linetype ordering within the family levels.
    keep_tools <- intersect(TOOL_ORDER, keep_tools)

    df_overall <- .load_trackA_overall(root) %>%
      filter(as.character(tool) %in% keep_tools) %>%
      mutate(
        tool   = factor(as.character(tool), levels = keep_tools),
        family = factor(TOOL_FAMILY[as.character(tool)], levels = fig4a_families)
      )

    # Sanity/verification hook: print mean AUPRC by family at L500 & L1000
    # so every render cross-checks the exact numbers cited in the paragraph.
    auprc_check <- df_overall %>%
      filter(length_bin %in% c("L500", "L1000")) %>%
      group_by(length_bin, family) %>%
      summarise(mean_auprc = mean(AUPRC, na.rm = TRUE), .groups = "drop")
    for (len in c("L500", "L1000")) {
      row <- auprc_check %>% filter(length_bin == len)
      if (nrow(row) > 0) {
        message(sprintf(
          "[fig4a] %s mean AUPRC: %s",
          len,
          paste(sprintf("%s = %.3f", row$family, row$mean_auprc),
                collapse = ", ")
        ))
      }
    }

    metrics_long <- df_overall %>%
      select(tool, family, length_bp, MCC, AUPRC, recall) %>%
      pivot_longer(cols = c(MCC, AUPRC, recall), names_to = "metric",
                   values_to = "value") %>%
      # The strip label is the raw column name, so `recall` printed lowercase
      # beside the two uppercase acronyms. Fixed here rather than with a labeller
      # so the panel order is set explicitly too, instead of falling out of
      # alphabetical sorting.
      mutate(metric = factor(metric,
                             levels = c("AUPRC", "MCC", "recall"),
                             labels = c("AUPRC", "MCC", "Recall")))

    family_colours <- FAMILY_PALETTE[fig4a_families]

    # Shape is a REDUNDANT third encoding of tool, alongside linetype. It is
    # kept because in the MCC and Recall facets the four green feature-based
    # curves cross repeatedly, and longdash/dotted/dotdash alone are hard to
    # follow through a crossing or in greyscale.
    #
    # Because it is redundant, it must appear in the SAME key as linetype, not
    # in a second "Tool" legend. ggplot merges two guides only when their title
    # and their guide spec match, so `tool_guide` below is handed to the
    # linetype AND the shape scale, identically. Do not set the shape scale to
    # guide = "none": that suppresses the shape guide, blocks the merge, and
    # leaves six marker styles on the curves with nothing defining them.
    #
    # The previous attempt to fix this -- override.aes = list(shape = ...) on
    # the linetype guide alone -- is a silent no-op and was removed. ggplot
    # draws a layer into a key only if that layer maps the guide's aesthetic;
    # geom_point() does not map linetype, so it never enters the Tool key and
    # there is no point glyph for the override to act on. Verified against the
    # built gtable: the Tool keys carried one glyph (the line), while the
    # Approach keys carried two. The colour override on the same call DID work,
    # which is why the breakage looked partial.
    tool_shapes <- setNames(c(16, 17, 15, 18, 8, 4), keep_tools)
    tool_guide <- guide_legend(
      # nrow = 2 keeps the six tool keys inside the print width when the legend
      # sits below the panel.
      order = 2, nrow = 2, byrow = TRUE,
      override.aes = list(colour = unname(family_colours[TOOL_FAMILY[keep_tools]]))
    )

    ggplot(metrics_long, aes(x = length_bp, y = value,
                             colour = family, linetype = tool,
                             group = tool)) +
      geom_line(linewidth = 0.6) +
      geom_point(aes(shape = tool), size = 1.4) +
      facet_wrap(~metric, scales = "free_y", ncol = 3) +
      scale_x_log10(labels = scales::comma_format()) +
      # "Approach" (not "Family"/"Method") is the shared legend title for
      # FAMILY_PALETTE across fig3B, fig4A, fig6C/D and fig_s9. Renamed from
      # "Family" on 2026-09-02 to stop it reading as taxonomic family.
      scale_colour_manual(
        values = family_colours, name = "Approach",
        # shape = NA strips the marker from the Approach keys. Both layers map
        # colour, so without it these keys draw a point too -- and the filled
        # circle it defaults to is ViraLM's marker in the panel, i.e. the only
        # shape in any legend would be attached to the wrong variable. Same idiom
        # as build_fig2b_pr_curves() in fig2.
        guide = guide_legend(order = 1, nrow = 2, byrow = TRUE,
                             override.aes = list(shape = NA))
      ) +
      scale_linetype_manual(
        values = setNames(c("solid", "longdash",
                            "solid", "longdash", "dotted", "dotdash"),
                          keep_tools),
        name = "Tool", guide = tool_guide
      ) +
      scale_shape_manual(
        values = tool_shapes,
        name = "Tool", guide = tool_guide
      ) +
      labs(x = "Fragment length (bp)", y = "Score") +
      theme_pub() +
      # Legend below the panel, not beside it: at print width the right-hand
      # legend cost ~17 % of the canvas and squeezed the three facets.
      theme(legend.position    = "bottom",
            legend.box         = "horizontal",
            legend.key.width   = unit(0.85, "cm"),
            legend.key.height  = unit(0.40, "cm"),
            legend.box.spacing = unit(0.15, "cm"),
            legend.margin      = margin(t = 0, r = 6, b = 0, l = 0),
            # Right margin widened from ggplot's 5.5 pt default. The last x break
            # is "10,000", centred on a tick close to the panel edge, so roughly
            # half the label sits outside the panel; at 5.5 pt the device clipped
            # its final "0". 16 pt clears the widest half-label at PT_MIN. The
            # other three sides keep the default.
            plot.margin        = margin(t = 5.5, r = 16, b = 5.5, l = 5.5))
  }

  p <- build_fig3a_metrics_vs_length(ROOT) + theme_print()
  save_fig(p, file.path(out_dir, "fig3.png"), w = 7.1, h = 3.2, scale = 1)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. 4  - Recall by virus category and declared scope (Q2)      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fig4 answers Q2 ("Published scope alone does not describe eukaryote-infecting
# virus recall"), so its tools are blocked by DECLARED SCOPE, not by the
# marker/reference split that fig5 uses for Q3.
#
# ONE plot at two lengths rather than stacked A/B panels: the tool labels, the
# colour bar and the scope headers are each drawn once, and a reader compares a
# tool across lengths by reading down a column. It also puts Jaeger's n/a
# column at 1,500 bp directly above its real values at 3,000 bp, which shows
# why the second length is there. 3,000 bp is the shortest Track A length at
# which all 14 tools are defined.

fig4 <- function(out_dir) {

  build_fig3b_category_recall <- function(root = ROOT,
                                          len = "L1500") {
    # `len` may be ONE length or several. With several, length becomes a facet ROW
    # and the panel shows the same tool at each length in one column, which is the
    # comparison the Q2 claim rests on. One plot rather than stacked A/B panels:
    # the tool labels, the colour bar and the scope headers are each drawn once,
    # and a reader compares a tool across lengths by reading down a column instead
    # of hopping between panels.
    df_cat <- .load_trackA_bycategory(root)
    df_overall <- .load_trackA_overall(root)

    # Runs that returned no classification are unavailable, NOT measured zero, and
    # the key is (tool, length): Jaeger is n/a at 1,500 bp but evaluable at 3,000,
    # so a tool-only mask would blank a row that has real data.
    no_output_keys <- df_overall %>%
      filter(!is.na(note), note == "no_output") %>%
      transmute(key = paste(tool, length_bin)) %>%
      pull(key)

    # Render-time verification hook (same idiom as the AUPRC check in
    # build_fig3a_metrics_vs_length): re-derives, per length, the tools whose
    # specificity falls below 0.50. The caption does not quote these values -- they
    # live in the Results sentence that cites this figure, and in Table S34 -- so
    # this hook is what keeps the prose and the data in step. A tool that misclassifies more than half the bacterial
    # negatives is close to calling everything viral, so a high recall row for it
    # tells you little -- that caveat lives in the caption, not on the panel.
    for (L in len) {
      chk <- df_overall %>%
        filter(length_bin == L, !is.na(specificity), specificity < 0.50) %>%
        arrange(specificity)
      if (nrow(chk) > 0) {
        message(sprintf(
          "[fig3b] %s specificity < 0.50 (Results text + Table S34 must match): %s", L,
          paste(sprintf("%s = %.3f", chk$tool, chk$specificity), collapse = ", ")))
      }
    }

    all_cats <- c(PHAGE_CATS, EUKARYOTIC_CATS)
    tool_order <- intersect(TOOL_ORDER, unique(as.character(df_cat$tool)))
    len_labels <- setNames(
      paste0(scales::comma(as.integer(sub("L", "", len))), " bp"), len)

    sub <- df_cat %>%
      filter(length_bin %in% len, category %in% all_cats) %>%
      mutate(
        cat_label = CAT_DISPLAY[as.character(category)],
        # Transposed: tools on x, virus categories on y. Keeps tiles near-square
        # at full page width and lets the genus/virus names read horizontally
        # instead of rotated 45 deg. rev() so the phage categories stay at the
        # top in PHAGE_CATS order, eukaryotic below.
        cat_label  = factor(cat_label, levels = rev(CAT_DISPLAY[all_cats])),
        tool       = factor(as.character(tool), levels = TOOL_ORDER),
        # DECLARED TARGET SCOPE drives the column blocks, because this figure
        # answers Q2 ("Published scope alone does not describe eukaryote-infecting
        # virus recall"). Grouping by TOOL_APPROACH instead would show it against
        # Q3's axis, which belongs to the homology figure.
        #
        # Scope is encoded by POSITION, not by axis-label colour as it was until
        # 2026-09-15. The colour version clashed with Fig. S17, which colours the
        # same 14 tools by tool identity: all 14 changed colour between the two
        # figures and 5 landed on a hue that denoted a different approach in the
        # other.
        target     = factor(TOOL_TARGET[as.character(tool)], levels = TARGET_LEVELS),
        length_bin = factor(as.character(length_bin), levels = len),
        recall     = ifelse(paste(tool, length_bin) %in% no_output_keys,
                            NA_real_, recall)
      )
    if (nrow(sub) == 0) stop("No category data for length(s) ", paste(len, collapse = ", "))
    stopifnot(!any(is.na(sub$target)))       # every plotted tool must be classified
    stopifnot(!any(is.na(sub$length_bin)))   # every requested length must have data

    p <- ggplot(sub, aes(x = tool, y = cat_label, fill = recall)) +
      geom_tile(colour = "white", linewidth = 0.4) +
      # 2 dp in BOTH recall heatmaps (this one and
      # build_fig4b_recall_vs_similarity). 1 dp collapsed distinctions the
      # Results quote directly -- e.g. geNomad 0.06 and VIBRANT 0.08 both
      # rendering "0.1". Per-category n is stated in the caption.
      geom_text(aes(label = ifelse(is.na(recall), "n/a",
                                   sprintf("%.2f", recall))),
                size = GEOM_TEXT_MIN) +
      # Identical to the fill scale in build_fig4b_recall_vs_similarity() so both
      # recall heatmaps are read on exactly the same colour ramp.
      scale_fill_gradient2(low = "#d73027", mid = "white", high = "#1a9850",
                           midpoint = 0.5, limits = c(0, 1), name = "Recall",
                           na.value = "grey85") +
      labs(x = NULL, y = NULL) +
      theme_pub() +
      theme(axis.text.x     = element_text(size = PT_MIN, angle = 45, hjust = 1),
            axis.text.y     = element_text(size = PT_MIN),
            # Strip fill, border and text come from theme_pub(); clip = "off" only
            # stops a descender being shaved and lets a label overhang a narrow
            # block. Copied from the homology heatmap.
            strip.clip      = "off",
            # 0.9 lines: wide enough that the block boundary reads as a break in
            # the matrix rather than as another white tile border.
            panel.spacing.x = unit(0.9, "lines"),
            panel.spacing.y = unit(0.5, "lines"))

    # space = "free_x" keeps every tile the same width across the 4/9/1 split;
    # without it the one-tile DB-dep. block would stretch to a third of the width
    # and read as wider-is-more.
    if (length(len) > 1) {
      p + facet_grid(
        length_bin ~ target, scales = "free_x", space = "free_x",
        labeller = labeller(target     = target_strip_labels(tool_order),
                            length_bin = len_labels))
    } else {
      p + facet_grid(~ target, scales = "free_x", space = "free_x",
                     labeller = as_labeller(target_strip_labels(tool_order)))
    }
  }

  p <- build_fig3b_category_recall(ROOT, c("L1500", "L3000")) +
    theme_print()
  save_fig(p, file.path(out_dir, "fig4.png"), w = 7.1, h = 6.8, scale = 1)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. 5  - Recall vs database similarity (ANI-novel panel, Q3)   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Heatmap: tool x {homology-detectable / homology-free} at L1500 and L3000.
#
# Terminology (aligned with manuscript Methods):
#   ANI-novel panel     = 15 vaginal phages, <95 % ANI to MetaVR v5
#   homology-detectable = 13 genomes retaining 84-95 % BLASTn relatives in MetaVR
#   homology-free       = 2 genomes with no MetaVR hit at any threshold
#   (distinct from Track C "viral dark matter" / "dark-matter contigs")
#
# Input: data/novel_spike_discovery/full_scale/dark_vs_relatives_stratification.tsv
#   (TSV column names retained: stratum in {has_relatives, dark} -- machine-readable
#    schema kept stable; only the user-facing labels are updated.)
#
# fig5 is a SINGLE panel. It held two before: a bar chart whose 28 values all
# reappeared in the heatmap, then briefly a delta panel. The delta was dropped
# because the two substrata are different genomes (13 and 2), so a diverging
# bar chart gave a between-genome difference the visual grammar of a
# within-subject effect, with n = 2 driving every sign. The claim is stronger
# stated in absolute terms, which the heatmap shows directly: the marker- and
# reference-dependent tools and the composition-based tools do not overlap on
# the homology-free row.

fig5 <- function(out_dir) {

  .load_fig5 <- function(root) {
    read_tsv(file.path(root,
             "data/novel_spike_discovery/full_scale/dark_vs_relatives_stratification.tsv"),
             show_col_types = FALSE)
  }

  # Tools that returned no classification at a given length, from the Track A
  # `note` column. Jaeger is flagged at L1500 because it requires fragments
  # >= 2,048 bp by construction, so the stratification TSV scores it TP = 0 /
  # FN = all. Plotting that as a measured recall of 0.00 would claim Jaeger ran
  # and detected nothing; it did not run. The category heatmap already masks this
  # case as "n/a" -- this keeps the two heatmaps telling the same story.
  .no_output_tools <- function(root, len) {
    f <- file.path(root, "results/spike_in_benchmark",
                   paste0(len, "_fragments"), "benchmark_results.tsv")
    if (!file.exists(f)) return(character(0))
    df <- read_tsv(f, show_col_types = FALSE)
    if (!"note" %in% names(df)) return(character(0))
    df %>% filter(!is.na(note), note == "no_output") %>%
      pull(tool) %>% as.character()
  }

  # Native ggplot heatmap (replaces the pheatmap version) so the legends can sit
  # horizontally on top like Panel A, the panel fills its full width to match A,
  # and the two tool groups the figure is about are separated on the x axis.
  #
  # Tools are split into two labelled blocks by TOOL_APPROACH (marker/reference-
  # dependent vs composition-based) rather than being coloured by the five-way
  # TOOL_FAMILY taxonomy used elsewhere. The result this figure carries is the
  # two-group one -- on the homology-free 1,500 bp row the groups do not overlap
  # -- and that is only legible if the groups are adjacent. Interleaved in the
  # cross-figure display order, the same 56 cells hide it. The five-way family
  # taxonomy is still carried by fig1, fig6 and Table S21; it was also drawn here
  # as label colour with no legend anywhere in the figure, so it could not be
  # decoded from fig4 alone.
  #
  # Within each block tools keep their global TOOL_ORDER, so the only thing that
  # moved relative to fig2/fig3 is the block boundary. TOOL_ORDER is a fixed
  # display order, NOT a metric ranking (see the note on it in functions.R), so
  # neither this figure nor its legend may describe it as an MCC ordering.
  build_fig4b_recall_vs_similarity <- function(root = ROOT) {
    df <- .load_fig5(root) %>%
      mutate(col_id = paste(stratum, length, sep = "_"))
    tool_order <- intersect(TOOL_ORDER, unique(df$tool))

    col_levels <- c("has_relatives_L1500", "dark_L1500",
                    "has_relatives_L3000", "dark_L3000")
    # Now ROW labels (the heatmap is transposed, see below), so they render
    # horizontally and can stay two-line without eating column width.
    col_labels <- c("Hom-detect.\nL1500", "Hom-free\nL1500",
                    "Hom-detect.\nL3000", "Hom-free\nL3000")
    keep <- col_levels %in% unique(df$col_id)
    col_levels <- col_levels[keep]; col_labels <- col_labels[keep]

    # Mask tool x length cells where the tool produced no output at that length
    # (Jaeger @ L1500), so they read "n/a" rather than a measured 0.00.
    masked <- lapply(unique(df$length),
                     function(l) tibble(length = l,
                                        tool   = .no_output_tools(root, l))) %>%
      bind_rows()

    # Dagger on geNomad, so the block header is not read as a claim that it scores
    # homology alone. See the TOOL_DUAL_BRANCH note in functions.R for why the flag
    # is separate-branch architecture and not "also scores non-homology features".
    lab_levels <- tool_approach_label(tool_order)

    d <- df %>%
      filter(tool %in% tool_order, col_id %in% col_levels) %>%
      mutate(recall = ifelse(paste(tool, length) %in%
                               paste(masked$tool, masked$length),
                             NA_real_, recall),
             # Transposed: tools on x (matching the tool order used everywhere
             # else), strata x length on y. Four columns stretched across a 7.1 in
             # page would have given 1.5 in tiles; this way the tiles are
             # near-square and the stratum labels read horizontally instead of
             # rotated 45 deg.
             approach = factor(TOOL_APPROACH[as.character(tool)],
                               levels = APPROACH_LEVELS),
             tool_lab = factor(tool_approach_label(tool), levels = lab_levels),
             col_id   = factor(col_id, levels = rev(col_levels)))

    stopifnot(!any(is.na(d$approach)))   # every plotted tool must be classified

    ggplot(d, aes(x = tool_lab, y = col_id)) +
      geom_tile(aes(fill = recall), colour = "white", linewidth = 0.4) +
      # 2 dp, matching build_fig3b_category_recall(); see the note there.
      geom_text(aes(label = ifelse(is.na(recall), "n/a",
                                   sprintf("%.2f", recall))),
                size = GEOM_TEXT_MIN) +
      # Identical to the fill scale in build_fig3b_category_recall() so the two
      # recall heatmaps are read on exactly the same colour ramp across figures.
      scale_fill_gradient2(low = "#d73027", mid = "white", high = "#1a9850",
                           midpoint = 0.5, limits = c(0, 1), name = "Recall",
                           na.value = "grey85") +
      scale_y_discrete(labels = rev(col_labels)) +
      # space = "free_x" keeps every tile the same width across the two blocks
      # despite the 6 vs 8 split; without it the six marker/reference tiles would
      # stretch to half the page and read as wider-is-more.
      facet_grid(~ approach, scales = "free_x", space = "free_x",
                 labeller = as_labeller(approach_strip_labels(tool_order))) +
      labs(x = NULL, y = NULL,
           caption = paste0("\u2020 scores sequence with a separate alignment-free ",
                            "branch as well as a marker database")) +
      theme_pub() +
      theme(
        legend.position  = "right",
        axis.text.x      = element_text(size = PT_MIN, angle = 45, hjust = 1),
        axis.text.y      = element_text(size = PT_MIN),
        panel.grid       = element_blank(),
        # Strip fill, border and text now come from theme_pub(); only the clip
        # behaviour is figure-specific. The labels are sized to fit their own
        # panel, so clip = "off" just stops a descender being shaved.
        strip.clip       = "off",
        # Wide enough that the block boundary reads as a break in the matrix
        # rather than as another white tile border.
        panel.spacing.x  = unit(0.9, "lines"),
        plot.caption     = element_text(size = PT_MIN, hjust = 0,
                                        colour = "grey30")
      )
  }

  p <- build_fig4b_recall_vs_similarity(ROOT) + theme_print()
  save_fig(p, file.path(out_dir, "fig5.png"), w = 7.1, h = 3.4, scale = 1)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. 6A-B  - Assembly-based benchmark (Track B)                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A: viral source genomes with >= 1 contig >= 1,500 bp by coverage.
# B: CST-IV-B minus CST-I delta MCC at 10x, nine selected tools
#    (build_fig5d_delta_ci() in functions.R, shared with Fig. S12).

fig6 <- function(out_dir) {

  # -- Panel A --
  build_fig5a_assembly_dynamics <- function(root = ROOT) {
    df_gt <- .load_trackB_assemblies(root)
    if (nrow(df_gt) == 0) stop("No Track B assembly ground truth files found")
    viral_counts <- df_gt %>%
      filter(label == "viral", length >= 1500) %>%
      group_by(background, CST, coverage) %>%
      summarise(n_viral = n_distinct(source_genome), .groups = "drop")

    ggplot(viral_counts, aes(x = factor(coverage), y = n_viral, fill = CST)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6,
               colour = "black", linewidth = 0.3) +
      scale_y_continuous(breaks = seq(0, 14, 2),
                         expand = expansion(mult = c(0, 0.05))) +
      scale_fill_manual(values = CST_COLORS) +
      labs(x = "Spike-in coverage (x)",
           y = "Viral source genomes assembled\n(>= 1 contig >= 1,500 bp)",
           fill = "Background") +
      theme_pub()
  }

  # -- Assemble and save --
  p <- build_fig5a_assembly_dynamics(ROOT) |
       build_fig5d_delta_ci(ROOT)
  p <- p + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 14))
  save_fig(p, file.path(out_dir, "fig6.png"), w = 16, h = 7)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. 7A-D  - Track C ground truth and ranking                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A: evidence-line co-occurrence heatmap.
# B: Phi29 enrichment ratio distribution (Track C) with GMM components.
# C: Track C MCC ranking with AUPRC.
# D: dark-matter detection rate per tool.

fig7 <- function(out_dir) {

  # -- Panel A --
  EVIDENCE_LABELS <- c(
    E1a = "Protein\nhomology", E1b = "Nucleotide\nhomology",
    E2  = "Structural\nsignatures", E3 = "Prophage\nintegration",
    E4  = "CRISPR\nspacers", E5 = "RCA cross-\nvalidation"
  )

  build_fig6a_evidence_cooccurrence <- function(root = ROOT) {
    df_gt <- .load_gt(root)
    ev_codes <- names(EVIDENCE_LABELS)
    for (ev in ev_codes) {
      df_gt[[ev]] <- as.integer(grepl(ev, df_gt$evidence_lines, fixed = TRUE))
    }
    cooccur <- matrix(0, nrow = length(ev_codes), ncol = length(ev_codes),
                      dimnames = list(ev_codes, ev_codes))
    for (i in seq_along(ev_codes)) {
      for (j in seq_along(ev_codes)) {
        cooccur[i, j] <- sum(df_gt[[ev_codes[i]]] & df_gt[[ev_codes[j]]])
      }
    }
    cooccur_pct <- cooccur / nrow(df_gt) * 100
    cooccur_pct[upper.tri(cooccur_pct)] <- NA
    cooccur_long <- reshape2::melt(cooccur_pct, na.rm = TRUE) %>%
      rename(row = Var1, col = Var2, pct = value)

    ggplot(cooccur_long, aes(x = col, y = row, fill = pct)) +
      geom_tile(colour = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.1f", pct)), size = 3.4) +
      scale_fill_gradient(low = "white", high = "#d62728", name = "% contigs") +
      scale_x_discrete(labels = EVIDENCE_LABELS) +
      scale_y_discrete(labels = EVIDENCE_LABELS) +
      # Both axes carry the same six evidence lines; the quantity is the fill /
      # cell value, not either axis. Naming the axes stops readers hunting for a
      # magnitude on them (raised on the assembled fig6, 2026-09-02).
      labs(x = "Evidence line", y = "Evidence line") +
      theme_pub() +
      theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 10),
            legend.title = element_text(size = 11),
            legend.text = element_text(size = 10),
            legend.key.size = unit(0.55, "cm"))
  }

  # -- Panel B --
  build_fig6c_enrichment_ratio <- function(root = ROOT) {
    rca_df <- .load_rca(root)
    rca_finite <- rca_df %>%
      mutate(enrichment_R = suppressWarnings(as.numeric(enrichment_R))) %>%
      filter(!is.na(enrichment_R), enrichment_R > 0, is.finite(enrichment_R)) %>%
      mutate(log_R = log10(enrichment_R))

    # 3-component Gaussian mixture on log10(R), fit in
    # scripts/07_evaluation/plot_fig_s13_track_c_gmm.py (sklearn, random_state = 0,
    # deterministic); parameters and the per-tool robustness check are in Table S23.
    .gmm <- data.frame(comp   = c("depleted", "baseline", "enriched"),
                       mean   = c(-1.909, -0.229, 2.111),
                       sd     = c(0.458, 0.486, 1.166),
                       weight = c(0.643, 0.340, 0.017))
    .gmm_cut <- log10(37.945)  # implied TP cutoff at P(enriched) > 0.95
    .nbins <- 80
    .binw  <- (max(rca_finite$log_R) - min(rca_finite$log_R)) / .nbins
    .nobs  <- nrow(rca_finite)
    .xg    <- seq(min(rca_finite$log_R), max(rca_finite$log_R), length.out = 400)
    comp_curves <- do.call(rbind, lapply(seq_len(nrow(.gmm)), function(i)
      data.frame(comp  = .gmm$comp[i],
                 log_R = .xg,
                 count = .gmm$weight[i] * dnorm(.xg, .gmm$mean[i], .gmm$sd[i]) * .nobs * .binw)))
    comp_curves$comp <- factor(comp_curves$comp,
                               levels = c("depleted", "baseline", "enriched"))

    ggplot(rca_finite, aes(x = log_R)) +
      # The true-negative band (0.5 < R < 2) is shaded rather than drawn as two
      # dotted lines, which frees the two labels those lines needed.
      annotate("rect", xmin = log10(0.5), xmax = log10(2),
               ymin = -Inf, ymax = Inf, fill = "#bdbdbd", alpha = 0.35) +
      # The TP side is shaded to the panel edge for the same reason: R > 10 is an
      # open-ended region, not a band with a right edge, and leaving it as a bare
      # dashed line made the TN band look like the only labelled region.
      annotate("rect", xmin = log10(10), xmax = Inf,
               ymin = -Inf, ymax = Inf, fill = "#d62728", alpha = 0.10) +
      # Observed distribution in neutral grey; colour is reserved for the mixture
      # components, which are the structure the panel is making a claim about.
      geom_histogram(bins = .nbins, fill = "grey88", colour = NA) +
      geom_area(data = comp_curves, aes(x = log_R, y = count, fill = comp),
                position = "identity", alpha = 0.60) +
      geom_line(data = comp_curves, aes(x = log_R, y = count, colour = comp),
                linewidth = 0.7, show.legend = FALSE) +
      scale_fill_manual(values = c(depleted = "#737373", baseline = "#4292c6",
                                   enriched = "#41ae76"),
                        name = "GMM component") +
      scale_colour_manual(values = c(depleted = "#404040", baseline = "#08519c",
                                     enriched = "#006d2c"),
                          guide = "none") +
      geom_vline(xintercept = log10(10), linetype = "dashed",
                 colour = "#d62728", linewidth = 0.8) +
      geom_vline(xintercept = .gmm_cut, linetype = "solid",
                 colour = "#54278f", linewidth = 0.8) +
      # Only the two cutoffs that carry the panel's argument are labelled: the
      # applied TP rule and the data-driven GMM cutoff it is checked against.
      annotate("text", x = -0.35, y = Inf, label = "TN band",
               vjust = 1.6, hjust = 0.5, size = 3.6, colour = "grey25") +
      annotate("text", x = 1.95, y = Inf, label = "TP band\nR > 10",
               vjust = 1.5, hjust = 0, size = 4, colour = "#d62728") +
      annotate("text", x = 1.95, y = Inf, label = "GMM cutoff\nR = 37.9",
               vjust = 4.6, hjust = 0, size = 4, colour = "#54278f") +
      labs(x = expression(log[10](R)), y = "Number of contigs") +
      theme_pub() +
      theme(legend.position = "right",
            legend.direction = "vertical",
            legend.title = element_text(size = 11),
            legend.text = element_text(size = 10),
            legend.key.width = unit(0.65, "cm"),
            legend.key.height = unit(0.45, "cm")) +
      guides(fill = guide_legend(ncol = 1, title.position = "top"))
  }

  # -- Panel C --
  build_fig7a_track_c_mcc_ranking <- function(root = ROOT) {
    metrics_file <- file.path(root,
      "results/test_real/coassembly/benchmark/overall_metrics.tsv")
    if (!file.exists(metrics_file)) {
      stop("Missing Track C overall_metrics.tsv at: ", metrics_file)
    }
    df <- read_tsv(metrics_file, show_col_types = FALSE)
    stopifnot(all(c("tool", "MCC", "AUPRC", "TP", "FP", "TN", "FN") %in% colnames(df)))
    # Derive displayed MCC from counts: rounding an already rounded four-decimal
    # score can change its three-decimal label (e.g. 0.599465 -> 0.5995 -> 0.600).
    df <- df %>% mutate(
      mcc_denom = sqrt((as.double(TP) + FP) * (TP + FN) * (TN + FP) * (TN + FN)),
      MCC = if_else(mcc_denom > 0, (TP * TN - FP * FN) / mcc_denom, 0)
    )
    df <- df %>%
      filter(tool %in% TOOL_ORDER) %>%
      mutate(family = TOOL_FAMILY[tool])
    df <- df %>% arrange(desc(MCC))
    tool_ord <- df$tool
    df <- df %>%
      mutate(tool = factor(tool, levels = rev(tool_ord)))

    df_long <- df %>%
      select(tool, family, MCC, AUPRC) %>%
      pivot_longer(cols = c(MCC, AUPRC), names_to = "metric", values_to = "value") %>%
      mutate(metric = factor(metric, levels = c("MCC", "AUPRC")))

    p <- ggplot(df_long, aes(x = value, y = tool, fill = family)) +
      geom_col(colour = "black", linewidth = 0.25, width = 0.7) +
      # Negative bars are tiny (min -0.049); drawing their labels leftward would
      # need ~0.3 axis units of otherwise-empty margin, so anchor them at zero and
      # let the bar direction plus the minus sign carry the sign.
      #
      # Sourmash's AUPRC is NA in overall_metrics.tsv: evaluate_rca_benchmark.py
      # writes "NA" whenever the tool's score vector has fewer than two distinct
      # values, and Sourmash made no calls at all on Track C (TP = FP = 0). Its
      # AUPRC is therefore undefined, not zero. geom_col silently drops that row,
      # leaving the AUPRC facet showing 13 of the 14 tools with no sign that one
      # is missing. Label the slot "n/a" -- same text style as the numeric labels
      # -- so the absence is visible and cannot be read as a measured 0.
      geom_text(aes(label = ifelse(is.na(value), "n/a", sprintf("%.3f", value)),
                    x = ifelse(is.na(value), 0, pmax(value, 0))),
                hjust = -0.1, size = 2.8) +
      facet_wrap(~metric, nrow = 1) +
      scale_fill_manual(values = FAMILY_PALETTE, name = "Approach") +
      scale_x_continuous(limits = c(-0.07, 0.86),
                         breaks = c(0, 0.2, 0.4, 0.6),
                         expand = expansion(mult = c(0, 0.03))) +
      labs(x = "Score (Track C: 29 TP, 931 TN)", y = NULL) +
      theme_pub() +
      theme(axis.text.y = element_text(
              colour = FAMILY_PALETTE[TOOL_FAMILY[rev(tool_ord)]]),
            legend.position = "right",
            legend.direction = "vertical",
            legend.title = element_text(size = 11),
            legend.text = element_text(size = 10),
            legend.key.size = unit(0.55, "cm")) +
      guides(fill = guide_legend(ncol = 1, title.position = "top"))
    p
  }

  # -- Panel D --
  # Headline dark-matter detection-rate bars only (single panel, no inner tag).
  # Used as a panel of the main multi-evidence figure (Fig. 6D); the per-contig
  # high-consensus view lives in Table S14.
  build_dark_matter_rates <- function(root = ROOT) {
    rates_path <- file.path(TABLES_DIR, "table_s15_dark_matter_detection_rates.tsv")
    if (!file.exists(rates_path)) stop("Missing table_s15_dark_matter_detection_rates.tsv")
    df_rates <- read_tsv(rates_path, show_col_types = FALSE) %>%
      arrange(detection_rate) %>%
      mutate(
        tool = factor(tool, levels = tool),
        family = TOOL_FAMILY[as.character(tool)]
      )
    p <- ggplot(df_rates, aes(x = detection_rate, y = tool)) +
      geom_col(aes(fill = family),
               colour = "black", linewidth = 0.25, width = 0.7) +
      geom_text(aes(label = paste0(n_detected, "/", n_dark_matter)),
                hjust = -0.1, size = 3) +
      geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey50") +
      scale_fill_manual(values = FAMILY_PALETTE, name = "Approach") +
      scale_x_continuous(limits = c(0, 1.1), labels = percent_format()) +
      labs(x = "Dark-matter detection rate", y = NULL) +
      theme_pub() +
      theme(legend.position = "right",
            legend.direction = "vertical",
            legend.title = element_text(size = 11),
            legend.text = element_text(size = 10),
            legend.key.size = unit(0.55, "cm")) +
      guides(fill = guide_legend(ncol = 1, title.position = "top"))
    colour_axis_by_family(p, axis = "y")
  }

  # -- Assemble and save --
  # 6B is wrapped so patchwork does not row-align it to 6A: 6A's rotated
  # two-line tick labels make a tall axis row, and aligning to it stranded
  # 6B's log10(R) title far below its own axis.
  top <- build_fig6a_evidence_cooccurrence(ROOT) |
         wrap_elements(full = build_fig6c_enrichment_ratio(ROOT))
  # 6C and 6D are the same family-coloured tool bars; collect their
  # identical Family legends into one shared vertical legend.
  bottom <- (build_fig7a_track_c_mcc_ranking(ROOT) |
             build_dark_matter_rates(ROOT)) +
            plot_layout(guides = "collect")
  p <- top / bottom
  p <- p + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 14))
  save_fig(p, file.path(out_dir, "fig7.png"), w = 15, h = 11)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. 8A-B  - External validation, 30-sample MiTCH cohort        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A: pooled per-tool MCC with sample-block bootstrap 95% CI, alongside the
#    matched primary-cohort (no-E5) MCC; see Table S22 for rank correlations.
# B: observational diversity contrast (high- vs low-diversity), delta MCC + CI.
# Data: results/test_real/viralm_investigation/mitch_scaleup30/scaleup30_*.tsv

fig8 <- function(out_dir) {

  # Tools shown in the main-figure ranking/diversity panels (top tier + mid tier);
  # the full 14-tool matrix is Fig. S13.
  MITCH_MAIN_TOOLS <- c("geNomad", "Jaeger", "VirSorter2", "VIBRANT", "ViraLM",
                        "PPR-Meta", "MetaPhinder", "DeepVirFinder")

  # -- Panel A --
  # Validation MCC (with CI) next to discovery MCC, ordered by validation MCC.
  build_fig8a_ranking <- function(root = ROOT) {
    rk <- .load_scaleup(root, "scaleup30_ranking_ci.tsv")
    ord <- rk %>% arrange(MCC) %>% pull(tool)   # ascending so highest plots at top
    long <- bind_rows(
      rk %>% transmute(tool, MCC, CI_low, CI_high, cohort = "Validation (MiTCH)"),
      rk %>% transmute(tool, MCC = disc_no_e5_MCC,
                       CI_low = NA_real_, CI_high = NA_real_, cohort = "Primary (UChoose)")
    ) %>% mutate(tool = factor(tool, levels = ord),
                 cohort = factor(cohort, levels = c("Validation (MiTCH)", "Primary (UChoose)")))

    ggplot(long, aes(x = MCC, y = tool, colour = cohort)) +
      geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.25,
                     linewidth = 0.5, na.rm = TRUE) +
      geom_point(size = 2.4) +
      scale_colour_manual(values = c("Validation (MiTCH)" = "#08519c",
                                     "Primary (UChoose)" = "#bdbdbd"), name = NULL) +
      labs(x = "Pooled MCC", y = NULL) +
      theme_pub()
  }

  # -- Panel B --
  # Observational diversity contrast (high minus low), delta MCC + 95% CI.
  build_fig8c_diversity <- function(root = ROOT) {
    div <- .load_scaleup(root, "scaleup30_diversity_contrast.tsv") %>%
      mutate(resolved = ifelse(as.character(delta_excludes_0) %in% c("TRUE", "True"),
                               "resolved", "not resolved")) %>%
      filter(tool %in% MITCH_MAIN_TOOLS) %>%
      mutate(tool = factor(tool, levels = rev(MITCH_MAIN_TOOLS)))

    ggplot(div, aes(x = delta_high_minus_low, y = tool)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_errorbarh(aes(xmin = delta_CIl, xmax = delta_CIh), height = 0.25,
                     linewidth = 0.5) +
      geom_point(aes(shape = resolved), size = 2.6) +
      # Display only: internal levels stay "resolved"/"not resolved" (set from
      # delta_excludes_0 above); the key states the criterion instead, matching
      # build_fig5d_delta_ci() in functions.R.
      scale_shape_manual(values = c("resolved" = 19, "not resolved" = 1),
                         labels = c("resolved"     = "CI excludes 0",
                                    "not resolved" = "CI includes 0"),
                         name = NULL) +
      labs(x = "ΔMCC (high − low diversity)", y = NULL) +
      theme_pub()
  }

  # -- Assemble and save --
  p <- build_fig8a_ranking(ROOT) |
       build_fig8c_diversity(ROOT)
  p <- p + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 14))
  save_fig(p, file.path(out_dir, "fig8.png"), w = 15, h = 8)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. 9A-B  - Compute cost does not predict accuracy             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A: Track A 1,500 bp accuracy-vs-cost space. x = wall-clock runtime (log),
#    y = MCC, point area = peak resident memory (log), colour = methodological
#    approach. The stepped line is the Pareto frontier: tools that no other
#    tool beats on BOTH axes. Everything below/right of it is dominated.
# B: per-contig runtime on the two Track B 10x backgrounds. The CST-IV-B
#    co-assembly holds 6.15x more contigs than the CST-I one (5,149 vs 837),
#    so raw wall-time ratios confound community state with input size.
#    Normalising per contig removes the CPU effect entirely (median ratio
#    0.98) and reveals the real signal: tools with a large fixed startup
#    cost (GPU model load, geNomad's marker database) get CHEAPER per contig
#    as the input grows.
#
# Scope of the dominance claim: one hardware configuration (AMD Zen4 CPU /
# NVIDIA L40s GPU), one allocation, one track, one fragment length, published
# default settings. It is ordinal on this configuration, not universal.
#
# Data (in TABLES_DIR): table_s16_resource_usage_L1500.tsv (runtime, RSS, VRAM)
#                       table_s34_track_a_overall_metrics_full.tsv (MCC, AUPRC)
#                       table_s8_track_b_metrics.tsv (contig counts)
#
# Rendered AT print size (7.1 in = 180 mm, scale = 1) like fig2-fig5, so
# declared point sizes are printed point sizes and the PT_MIN floor holds.
# Panel B exists to defuse a confound: the raw CST-IV-B/CST-I wall-time ratio
# (median 6.0x for CPU tools) is almost exactly the 6.15x difference in contig
# count between the two co-assemblies, so it is input size, not community
# complexity. Plotting per-contig cost is what makes that visible.

fig9 <- function(out_dir) {

  # Track B 10x background identities (coverage sweep backgrounds; Methods).
  BG_CST_I    <- "UC115_V2"
  BG_CST_IV_B <- "UC093_V3"

  .load_resource <- function(root) {
    f <- file.path(TABLES_DIR, "table_s16_resource_usage_L1500.tsv")
    if (!file.exists(f)) stop("Missing table_s16_resource_usage_L1500.tsv")
    read_tsv(f, comment = "#", show_col_types = FALSE)
  }

  .load_track_a <- function(root) {
    f <- file.path(TABLES_DIR, "table_s34_track_a_overall_metrics_full.tsv")
    if (!file.exists(f)) stop("Missing table_s34_track_a_overall_metrics_full.tsv")
    read_tsv(f, show_col_types = FALSE)
  }

  # Contig counts actually handed to the tools for each 10x co-assembly.
  .trackb_contig_counts <- function(root) {
    f <- file.path(TABLES_DIR, "table_s8_track_b_metrics.tsv")
    if (!file.exists(f)) stop("Missing table_s8_track_b_metrics.tsv")
    read_tsv(f, show_col_types = FALSE) %>%
      filter(coverage == 10, background %in% c(BG_CST_I, BG_CST_IV_B)) %>%
      distinct(background, n_total)
  }

  # Non-dominated set on (low runtime, high MCC).
  .pareto <- function(df) {
    keep <- vapply(seq_len(nrow(df)), function(i) {
      !any(df$wall_time_sec <= df$wall_time_sec[i] & df$MCC >= df$MCC[i] &
             (df$wall_time_sec < df$wall_time_sec[i] | df$MCC > df$MCC[i]))
    }, logical(1))
    df[keep, ] %>% arrange(wall_time_sec)
  }

  # -- Panel A --
  # Panel A: accuracy-vs-cost space with the Pareto frontier.
  build_fig9a_cost_accuracy <- function(root = ROOT) {
    res <- .load_resource(root) %>%
      filter(condition == "TrackA_L1500", exit_code == 0)   # drops Jaeger no-op
    dat <- res %>%
      inner_join(.load_track_a(root) %>% select(tool, MCC), by = "tool") %>%
      mutate(
        family = TOOL_FAMILY[tool],
        hw     = ifelse(gpu_required == "yes", "GPU", "CPU")
      )

    front <- .pareto(dat)
    n_dom <- nrow(dat) - nrow(front)

    ggplot(dat, aes(x = wall_time_sec, y = MCC)) +
      geom_step(data = front, direction = "vh",
                colour = "grey45", linetype = "22", linewidth = 0.45) +
      geom_point(aes(size = peak_rss_mb, fill = family, shape = hw),
                 colour = "grey20", stroke = 0.35, alpha = 0.9) +
      ggrepel::geom_text_repel(aes(label = tool), size = GEOM_TEXT_MIN,
                               min.segment.length = 0.15, seed = 42,
                               box.padding = 0.38, max.overlaps = 20,
                               segment.colour = "grey55", segment.size = 0.3) +
      annotate("text", x = 28, y = -0.09, hjust = 0, size = GEOM_TEXT_MIN,
               colour = "grey30", fontface = "italic",
               label = sprintf("%d of %d tools dominated", n_dom, nrow(dat))) +
      scale_x_log10(breaks = c(30, 60, 300, 600, 3600),
                    labels = c("30 s", "1 min", "5 min", "10 min", "1 h")) +
      # Labelled in MB, the unit Table S16 records, so the figure does not commit
      # to a decimal-vs-binary GB conversion (the full recommendation table
      # reports geNomad as ~18.5 GB, i.e. peak_rss_mb / 1024).
      scale_size_continuous(trans = "log10", range = c(1.8, 8),
                            breaks = c(100, 1000, 10000),
                            labels = c("100", "1,000", "10,000"),
                            name = "Peak RAM (MB)") +
      scale_fill_manual(values = FAMILY_PALETTE, name = "Approach") +
      scale_shape_manual(values = c(CPU = 21, GPU = 24), name = NULL) +
      guides(fill  = guide_legend(order = 1, override.aes = list(shape = 21, size = 3.2)),
             shape = guide_legend(order = 2, override.aes = list(size = 3.2)),
             size  = guide_legend(order = 3)) +
      labs(x = "Wall-clock runtime, 1,500 bp panel (log scale)",
           y = "MCC (Track A, 1,500 bp)") +
      theme_pub() +
      theme(legend.key.height = unit(0.8, "lines"))
  }

  # -- Panel B --
  # Panel B: does the CST-IV-B slowdown survive normalisation by contig count?
  # Each tool gets its raw wall-time ratio and its per-contig ratio on one axis.
  # The raw ratios sit on the contig-count line (6.15x), i.e. they are fully
  # explained by input size; the per-contig ratios collapse onto 1.0 for CPU
  # tools (no residual complexity effect) and fall below it for tools with a
  # large fixed startup cost.
  build_fig9b_per_contig_scaling <- function(root = ROOT) {
    nb <- .trackb_contig_counts(root)
    n_i   <- nb$n_total[nb$background == BG_CST_I]
    n_ivb <- nb$n_total[nb$background == BG_CST_IV_B]
    size_ratio <- n_ivb / n_i

    w <- .load_resource(root) %>%
      filter(condition %in% c("TrackB_CST-I_10x", "TrackB_CST-IV-B_10x")) %>%
      select(tool, gpu_required, condition, wall_time_sec) %>%
      pivot_wider(names_from = condition, values_from = wall_time_sec) %>%
      mutate(
        raw  = .data[["TrackB_CST-IV-B_10x"]] / .data[["TrackB_CST-I_10x"]],
        norm = raw / size_ratio,
        hw   = ifelse(gpu_required == "yes", "GPU tools", "CPU tools")
      ) %>%
      arrange(norm) %>%
      mutate(tool = factor(tool, levels = tool))

    long <- w %>%
      select(tool, hw, raw, norm) %>%
      pivot_longer(c(raw, norm), names_to = "measure", values_to = "ratio") %>%
      mutate(measure = factor(measure, levels = c("raw", "norm"),
                              labels = c("Raw wall-time ratio",
                                         "Normalised per contig")))

    # Median normalised ratio, pinned to the left margin of each facet's top row
    # (the region left of 0.3 is empty except for the bottom-row tools).
    med <- w %>% group_by(hw) %>%
      summarise(m = median(norm), top = as.character(tool[which.max(norm)]),
                .groups = "drop") %>%
      mutate(lab = sprintf("median %.2f", m),
             tool = factor(top, levels = levels(w$tool)), x = 0.105)

    # Reference-line captions: drawn once, in the CPU facet, above the top row.
    # No `tool` column -- carrying one would inject a phantom row into the facet.
    refs <- data.frame(
      hw  = "CPU tools",
      x   = c(1, size_ratio),
      h   = c(1.05, -0.05),
      lab = c("1.0\nper-contig cost\nunchanged",
              sprintf("%.2f\ncontig-count\nratio", size_ratio))
    )

    ggplot(long, aes(x = ratio, y = tool)) +
      geom_vline(xintercept = 1, colour = "grey60", linewidth = 0.45) +
      geom_vline(xintercept = size_ratio, colour = "grey60",
                 linetype = "22", linewidth = 0.45) +
      geom_line(aes(group = tool), colour = "grey78", linewidth = 0.45) +
      geom_point(aes(colour = measure), size = 2.1) +
      geom_text(data = refs, aes(x = x, y = Inf, label = lab, hjust = h),
                inherit.aes = FALSE, vjust = 1.05, size = GEOM_TEXT_MIN,
                colour = "grey35", fontface = "italic", lineheight = 0.95) +
      geom_text(data = med, aes(x = x, y = tool, label = lab), inherit.aes = FALSE,
                hjust = 0, size = GEOM_TEXT_MIN, colour = "#08519c") +
      facet_grid(hw ~ ., scales = "free_y", space = "free_y") +
      scale_x_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
                    labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10"),
                    expand = expansion(mult = c(0.04, 0.26))) +
      scale_y_discrete(expand = expansion(add = c(0.6, 1.9))) +
      scale_colour_manual(values = c("Raw wall-time ratio" = "grey55",
                                     "Normalised per contig" = "#08519c"),
                          name = NULL) +
      labs(x = "CST-IV-B / CST-I runtime ratio at 10x (log scale)", y = NULL) +
      theme_pub() +
      theme(panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3),
            strip.background = element_blank(),
            strip.text.y = element_text(angle = -90, face = "bold"),
            legend.position = "top",
            legend.margin = margin(b = -4))
  }

  # -- Assemble and save --
  p <- build_fig9a_cost_accuracy(ROOT) /
       build_fig9b_per_contig_scaling(ROOT)
  p <- p + plot_layout(heights = c(1.25, 1)) +
    plot_annotation(tag_levels = "A") &
    theme_print()
  save_fig(p, file.path(out_dir, "fig9.png"), w = 7.1, h = 7.4, scale = 1)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. S1  - Contig lengths after the 1,500 bp filter             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Contig length distribution across the 13 shotgun metagenome assemblies.

fig_s1 <- function(out_dir) {

  build_fig_s1_contig_length_distribution <- function(root = ROOT,
                                                      min_length = 1500) {
    gt_dir <- file.path(root, "results/test_real/ground_truth")
    # Restrict to the 13 primary UChoose samples (SAMPLE_CST in functions.R).
    # The ground_truth directory also holds MiTCH and RCA-pool folders, which
    # are not part of this figure's cohort.
    samples <- file.path(gt_dir, names(SAMPLE_CST))
    data_list <- list()
    for (sdir in samples) {
      sample <- basename(sdir)
      gt_file <- file.path(sdir, "ground_truth_with_kraken2.tsv")
      if (!file.exists(gt_file)) stop("Fig. S1: missing ground truth for ", sample)
      df <- read_tsv(gt_file, show_col_types = FALSE)
      if (!"length" %in% colnames(df)) stop("Fig. S1: no length column for ", sample)
      data_list[[sample]] <- df %>%
        filter(length >= min_length) %>%
        mutate(sample = sample, CST = SAMPLE_CST[sample])
    }
    data <- bind_rows(data_list)

    cst_rank <- c("CST-I" = 1, "CST-III" = 2, "CST-IV" = 3)
    sample_order <- data %>%
      distinct(sample, CST) %>%
      arrange(cst_rank[CST], sample) %>%
      pull(sample)
    data <- data %>% mutate(sample = factor(sample, levels = sample_order))

    ggplot(data, aes(x = sample, y = length, fill = CST)) +
      geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2, linewidth = 0.4) +
      geom_hline(yintercept = min_length, linetype = "dashed",
                 colour = "grey50", linewidth = 0.6) +
      scale_y_log10(labels = scales::comma_format()) +
      scale_fill_manual(values = CST_COLORS) +
      labs(
        x = NULL,
        y = "Contig length (bp)",
        fill = "CST",
        subtitle = sprintf("%s contigs >= %s bp across %d samples",
                           format(nrow(data), big.mark = ","),
                           format(min_length, big.mark = ","),
                           length(sample_order))
      ) +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
  }

  p <- build_fig_s1_contig_length_distribution(ROOT)
  save_fig(p, file.path(out_dir, "fig_s1.png"), w = 12, h = 5)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. S2  - Precision vs recall at 3,000 bp (Track A)            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fig_s2 <- function(out_dir) {

  build_fig2b_precision_recall <- function(root = ROOT, len = "L1500") {
    df <- .load_fig2_trackA_overall(root)
    sub <- df %>% filter(length_bin == len)
    if (nrow(sub) == 0) stop("No data for length ", len)
    sub$family <- TOOL_FAMILY[as.character(sub$tool)]

    ggplot(sub, aes(x = recall, y = precision, colour = family)) +
      geom_point(size = 3, stroke = 0.5) +
      ggrepel::geom_text_repel(aes(label = tool), size = 3, max.overlaps = 20) +
      # Legend title is "Approach" everywhere the FAMILY_PALETTE is used (fig4A,
      # fig6C/D, fig_s9) -- keep it uniform here too.
      scale_colour_manual(values = FAMILY_PALETTE, name = "Approach") +
      scale_x_continuous(limits = c(0, 1)) +
      scale_y_continuous(limits = c(0, 1)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey60") +
      labs(x = "Recall", y = "Precision") +
      theme_pub()
  }

  p <- build_fig2b_precision_recall(ROOT, "L3000")
  save_fig(p, file.path(out_dir, "fig_s2.png"), w = 8, h = 6)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. S3  - Pairwise McNemar significance heatmap at 1,500 bp    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Demoted from former main Fig. 3C.

fig_s3 <- function(out_dir) {

  build_fig_s3_mcnemar_heatmap <- function(root = ROOT,
                                           len = "L1500") {
    input_path <- file.path(TABLES_DIR,
      "table_s6_mcnemar_significant_pairs.tsv")
    if (!file.exists(input_path)) stop("Missing table_s6_mcnemar_significant_pairs.tsv")
    df <- read_tsv(input_path, show_col_types = FALSE)
    df_l <- df %>% filter(grepl(paste0("^", len), length_bin))

    tools <- intersect(TOOL_ORDER, unique(c(df_l$tool_a, df_l$tool_b)))
    # Cell value is -log10(q) clipped at CAP so the visible gradient is not
    # crushed by floating-point underflow (q == 0 -> -log10(q) = Inf, otherwise
    # would saturate at ~300 and make the heatmap appear binary).
    CAP <- 20
    mat <- matrix(0, nrow = length(tools), ncol = length(tools),
                  dimnames = list(tools, tools))
    qmat <- matrix(NA_real_, nrow = length(tools), ncol = length(tools),
                   dimnames = list(tools, tools))
    for (i in seq_len(nrow(df_l))) {
      a <- df_l$tool_a[i]; b <- df_l$tool_b[i]
      if (a %in% tools && b %in% tools) {
        q  <- df_l$q_value[i]
        sig <- isTRUE(df_l$significant_bh[i])
        val <- if (sig) min(-log10(max(q, 1e-300)), CAP) else 0
        mat[a, b] <- val
        mat[b, a] <- val
        qmat[a, b] <- q
        qmat[b, a] <- q
      }
    }
    diag(mat) <- NA

    # Selective annotation: print q-value text only for cells where the gradient
    # actually carries information (significant but q > 1e-6, i.e., not at the
    # underflow ceiling). All other cells show no number.
    label_mat <- matrix("", nrow = length(tools), ncol = length(tools),
                        dimnames = list(tools, tools))
    for (i in seq_len(nrow(df_l))) {
      a <- df_l$tool_a[i]; b <- df_l$tool_b[i]
      if (a %in% tools && b %in% tools) {
        q   <- df_l$q_value[i]
        sig <- isTRUE(df_l$significant_bh[i])
        if (sig && is.finite(q) && q > 1e-6) {
          lbl <- if (q >= 0.001) sprintf("%.3f", q) else sprintf("%.0e", q)
          label_mat[a, b] <- lbl
          label_mat[b, a] <- lbl
        }
      }
    }

    row_annot <- data.frame(Family = TOOL_FAMILY[tools], row.names = tools)
    ann_colors <- list(Family = FAMILY_PALETTE)

    # Force the colour breaks to span 0..CAP exactly so the legend and gradient
    # line up; without explicit breaks pheatmap rescales to data range, which
    # collapses the gradient when most cells are at the cap.
    breaks <- seq(0, CAP, length.out = 101)

    ph <- pheatmap(
      mat,
      color = colorRampPalette(c("white", "#9ecae1", "#08519c"))(100),
      breaks = breaks,
      na_col = "grey90",
      cluster_rows = FALSE, cluster_cols = FALSE,
      display_numbers = label_mat, number_color = "black",
      fontsize_number = 7,
      border_color = "white",
      annotation_row = row_annot,
      annotation_colors = ann_colors,
      fontsize = 9, silent = TRUE,
      legend_breaks = c(0, 1.3, 5, 10, 15, 20),
      legend_labels = c("NS", "0.05", "1e-5", "1e-10", "1e-15", "<=1e-20")
    )
    ggplotify::as.ggplot(ph)
  }

  p <- build_fig_s3_mcnemar_heatmap(ROOT, "L1500")
  save_fig(p, file.path(out_dir, "fig_s3.png"), w = 9, h = 8)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. S4  - Threshold calibration at 1,500 bp                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Default-threshold MCC vs best-threshold MCC per tool, showing how much of
# each tool's apparent gap can be recovered by re-tuning its score threshold.
# Source: table_s41_threshold_calibration.tsv in TABLES_DIR.

fig_s4 <- function(out_dir) {

  build_fig_s4_threshold_calibration <- function(root = ROOT,
                                                  len = "L1500") {
    input_path <- file.path(TABLES_DIR,
      "table_s41_threshold_calibration.tsv")
    if (!file.exists(input_path)) stop("Missing table_s41_threshold_calibration.tsv")
    df <- read_tsv(input_path, show_col_types = FALSE)
    df <- df %>%
      filter(grepl(paste0("^", len), length_bin),
             tool %in% TOOL_ORDER) %>%
      mutate(
        tool     = factor(tool, levels = rev(TOOL_ORDER)),
        family   = TOOL_FAMILY[as.character(tool)]
      )

    if (nrow(df) == 0) stop("No calibration data for ", len)

    df_long <- df %>%
      select(tool, family, default_MCC, best_MCC, best_threshold, default_threshold) %>%
      pivot_longer(c(default_MCC, best_MCC),
                   names_to = "variant", values_to = "MCC") %>%
      mutate(variant = recode(variant,
                              default_MCC = "Benchmark rule",
                              best_MCC    = "In-sample optimum"),
             variant = factor(variant,
                              levels = c("Benchmark rule", "In-sample optimum")))

    ggplot(df_long, aes(x = MCC, y = tool, fill = variant)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6,
               colour = "black", linewidth = 0.25) +
      geom_text(data = df %>% filter(is.na(best_MCC)),
                aes(x = 0.025, y = tool, label = "n/a"),
                inherit.aes = FALSE, hjust = 0, size = 3) +
      scale_fill_manual(values = c("Benchmark rule"  = "#999999",
                                   "In-sample optimum" = "#2980b9"),
                        name = NULL) +
      scale_x_continuous(limits = c(-0.05, 1), expand = expansion(mult = c(0, 0.05))) +
      labs(
        x = sprintf("MCC on Track A at %s bp",
                    format(as.integer(sub("^L", "", len)), big.mark = ",")),
        y = NULL,
        subtitle = sprintf("Benchmark vs in-sample optimized MCC (%d/%d tools threshold-evaluable)",
                           sum(!is.na(df$best_MCC)), nrow(df))
      ) +
      theme_pub() +
      theme(
        axis.text.y      = element_text(
          colour = FAMILY_PALETTE[TOOL_FAMILY[rev(TOOL_ORDER)]],
          size = 8
        ),
        legend.position  = "top",
        plot.subtitle    = element_text(size = 9)
      )
  }

  p <- build_fig_s4_threshold_calibration(ROOT, "L1500")
  save_fig(p, file.path(out_dir, "fig_s4.png"), w = 9, h = 7)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. S5  - Eukaryote-infecting virus recall across lengths      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Manuscript citation: "herpesviruses were the hardest target, with only
# ViraLM (100%) and TransGINmer (81%) exceeding 50% recall at L1500
# (Supplementary Fig. S5)".
# Source: results/spike_in_benchmark/L*_fragments/benchmark_results_by_category.tsv

fig_s5 <- function(out_dir) {

  build_fig_s5_eukaryotic_recall <- function(root = ROOT) {
    frames <- list()
    for (len in ALL_LENGTHS) {
      catf <- file.path(root, "results/spike_in_benchmark",
                        paste0(len, "_fragments"),
                        "benchmark_results_by_category.tsv")
      overf <- file.path(root, "results/spike_in_benchmark",
                         paste0(len, "_fragments"),
                         "benchmark_results.tsv")
      if (!file.exists(catf)) next
      dfc <- read_tsv(catf, show_col_types = FALSE) %>%
        mutate(length_bin = len)
      # Flag tool-length combos with no output and propagate as NA
      no_out <- character(0)
      if (file.exists(overf)) {
        dfo <- read_tsv(overf, show_col_types = FALSE)
        if ("note" %in% colnames(dfo)) {
          no_out <- dfo %>%
            filter(!is.na(note), note == "no_output") %>%
            pull(tool) %>% as.character()
        }
      }
      dfc <- dfc %>%
        mutate(recall = ifelse(as.character(tool) %in% no_out,
                               NA_real_, recall))
      frames[[len]] <- dfc
    }
    df <- bind_rows(frames) %>%
      filter(category %in% EUKARYOTIC_CATS, tool %in% TOOL_ORDER) %>%
      mutate(
        tool       = factor(tool, levels = rev(TOOL_ORDER)),
        length_bin = factor(length_bin, levels = ALL_LENGTHS),
        cat_label  = CAT_DISPLAY[as.character(category)],
        cat_label  = factor(cat_label, levels = CAT_DISPLAY[EUKARYOTIC_CATS])
      )

    panels <- lapply(levels(df$cat_label), function(cl) {
      sub <- df %>% filter(cat_label == cl)
      ggplot(sub, aes(x = length_bin, y = tool, fill = recall)) +
        geom_tile(colour = "white", linewidth = 0.5) +
        geom_text(aes(label = ifelse(is.na(recall), "n/a",
                                     sprintf("%.1f", recall))),
                  size = 2.6) +
        scale_fill_gradient2(low = "#d73027", mid = "white", high = "#1a9850",
                             midpoint = 0.5, limits = c(0, 1), name = "Recall",
                             na.value = "grey85") +
        labs(x = NULL, y = NULL, title = cl) +
        theme_pub() +
        theme(
          axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
          axis.text.y = element_text(
            size = 8,
            colour = FAMILY_PALETTE[TOOL_FAMILY[rev(TOOL_ORDER)]]
          ),
          legend.position = "right"
        )
    })

    wrap_plots(panels, nrow = 1, guides = "collect") &
      theme(legend.position = "right")
  }

  p <- build_fig_s5_eukaryotic_recall(ROOT)
  save_fig(p, file.path(out_dir, "fig_s5.png"), w = 16, h = 7)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. S6  - Track B recall by spike-in coverage                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fig_s6 <- function(out_dir) {

  build_fig_s7_trackb_recall_by_coverage <- function(root = ROOT) {
    cov_file <- file.path(TABLES_DIR, "table_s10_track_b_recall_by_coverage.tsv")
    if (!file.exists(cov_file)) stop("Missing table_s10_track_b_recall_by_coverage.tsv")
    df_cov <- read_tsv(cov_file, show_col_types = FALSE)
    stopifnot(all(c("tool", "coverage", "recall") %in% colnames(df_cov)))
    df_cov <- df_cov %>%
      mutate(tool = factor(tool, levels = TOOL_ORDER),
             CST = ifelse(grepl("115", background), "CST-I", "CST-IV-B"))

    ggplot(df_cov, aes(x = factor(coverage), y = recall,
                       colour = tool, group = tool)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2) +
      facet_wrap(~CST) +
      scale_colour_manual(values = TOOL_COLORS, name = "Tool") +
      labs(x = "Spike-in coverage (x)", y = "Recall") +
      theme_pub() +
      theme(legend.text = element_text(size = 7))
  }

  p <- build_fig_s7_trackb_recall_by_coverage(ROOT)
  save_fig(p, file.path(out_dir, "fig_s6.png"), w = 14, h = 6)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. S7  - RCA-validated recall by virus category               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fig_s7 <- function(out_dir) {

  build_fig_s8_rca_category <- function(root = ROOT) {
    input_path <- file.path(TABLES_DIR,
      "table_s3_track_c_by_virus_category.tsv")
    if (!file.exists(input_path)) stop("Missing table_s3_track_c_by_virus_category.tsv")
    df <- read_tsv(input_path, show_col_types = FALSE)

    cat_col <- intersect(c("virus_category", "category", "viral_category"),
                         colnames(df))
    if (length(cat_col) == 0) stop("No category column in table_s8")
    recall_col <- intersect(c("rca_validated_recall", "recall",
                              "agg_rca_recall", "rca_recall"),
                            colnames(df))
    if (length(recall_col) == 0) stop("No recall column in table_s8")

    # Aggregate across samples: weight by n_rca_positive (true positives seen)
    # so per-tool × category recall is the micro-recall.
    df <- df %>%
      rename(cat = !!cat_col[1], rec = !!recall_col[1]) %>%
      mutate(rec = as.numeric(rec))

    tool_order <- intersect(TOOL_ORDER, unique(df$tool))

    if ("cell_A" %in% colnames(df) && "cell_D" %in% colnames(df)) {
      agg <- df %>%
        filter(tool %in% tool_order) %>%
        group_by(tool, cat) %>%
        summarise(rec = sum(cell_A, na.rm = TRUE) /
                        pmax(sum(cell_A, na.rm = TRUE) +
                             sum(cell_D, na.rm = TRUE), 1),
                  .groups = "drop")
    } else {
      agg <- df %>%
        filter(tool %in% tool_order) %>%
        group_by(tool, cat) %>%
        summarise(rec = mean(rec, na.rm = TRUE), .groups = "drop")
    }

    agg <- agg %>%
      mutate(tool = factor(tool, levels = rev(tool_order)))

    ggplot(agg, aes(x = rec, y = tool, colour = cat)) +
      geom_point(size = 3.5) +
      geom_line(aes(group = tool), colour = "grey70", linewidth = 0.4) +
      # Raw column values ("eukaryotic_virus", "phage") were reaching the legend
      # verbatim, snake_case and all. Relabelled here rather than by mutating the
      # data so the join keys stay untouched.
      scale_colour_brewer(palette = "Set1", name = "Virus category",
                          labels = c(eukaryotic_virus = "Eukaryote-infecting virus",
                                     phage             = "Phage")) +
      scale_x_continuous(limits = c(0, 1), labels = percent_format()) +
      labs(x = "RCA-validated recall (13 samples, 11 participants)", y = NULL) +
      theme_pub() +
      theme(
        axis.text.y = element_text(
          colour = SCOPE_LABEL_COLORS[TOOL_SCOPE[rev(tool_order)]]
        )
      )
  }

  p <- build_fig_s8_rca_category(ROOT)
  save_fig(p, file.path(out_dir, "fig_s7.png"), w = 8, h = 7)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. S8  - Track C master-contig lengths with N50               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Matches manuscript "pooled N50 22,105 bp ... (Fig. S8)". The 16,672 bp
# figure an earlier comment carried until 2026-08-11 was UC115_V2's single-
# participant co-assembly N50, not the pooled value this panel computes.

fig_s8 <- function(out_dir) {

  .compute_N50 <- function(lengths) {
    s <- sort(lengths, decreasing = TRUE)
    half <- sum(s) / 2
    cs <- cumsum(s)
    s[which(cs >= half)[1]]
  }

  build_fig7b_track_c_contig_length <- function(root = ROOT) {
    gt_file <- file.path(root, "results/test_real/coassembly/rca_ground_truth.tsv")
    if (!file.exists(gt_file)) {
      stop("Missing rca_ground_truth.tsv at: ", gt_file)
    }
    df <- read_tsv(gt_file, show_col_types = FALSE)
    stopifnot(all(c("contig_length", "label", "patient") %in% colnames(df)))

    df <- df %>%
      mutate(
        contig_length = as.numeric(contig_length),
        label = factor(label,
                       levels = c("TP", "TN", "dark_matter", "excluded"),
                       labels = c("Viral TP (R>10)",
                                  "Bacterial TN",
                                  "Dark matter",
                                  "Excluded"))
      ) %>%
      filter(!is.na(contig_length), contig_length > 0)

    n50_overall <- .compute_N50(df$contig_length)
    n50_tp      <- .compute_N50(df$contig_length[df$label == "Viral TP (R>10)"])
    median_len  <- median(df$contig_length, na.rm = TRUE)

    label_palette <- c("Viral TP (R>10)" = "#d62728",
                       "Bacterial TN"    = "#2980b9",
                       "Dark matter"     = "#7f3b08",
                       "Excluded"        = "grey70")

    ggplot(df, aes(x = contig_length, fill = label)) +
      geom_density(alpha = 0.45, linewidth = 0.3) +
      geom_vline(xintercept = n50_overall, linetype = "dashed",
                 colour = "black", linewidth = 0.5) +
      geom_vline(xintercept = n50_tp, linetype = "dashed",
                 colour = "#d62728", linewidth = 0.5) +
      annotate("text",
               x = n50_overall, y = Inf,
               label = sprintf("N50 (all) = %s bp",
                               format(round(n50_overall), big.mark = ",")),
               vjust = 1.5, hjust = -0.05, size = 2.8, colour = "black") +
      annotate("text",
               x = n50_tp, y = Inf,
               label = sprintf("N50 (TP) = %s bp",
                               format(round(n50_tp), big.mark = ",")),
               vjust = 3.0, hjust = -0.05, size = 2.8, colour = "#d62728") +
      scale_x_log10(labels = scales::comma_format()) +
      scale_fill_manual(values = label_palette, name = NULL) +
      labs(
        x = "Master-contig length (bp, log scale)",
        y = "Density",
        subtitle = sprintf("Median length = %s bp;  %s contigs across 13 samples (11 participants)",
                           format(round(median_len), big.mark = ","),
                           format(nrow(df), big.mark = ","))
      ) +
      theme_pub() +
      theme(legend.position = "top",
            plot.subtitle = element_text(size = 9))
  }

  p <- build_fig7b_track_c_contig_length(ROOT)
  save_fig(p, file.path(out_dir, "fig_s8.png"), w = 10, h = 5)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. S9  - RCA threshold-sensitivity recall heatmap             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Demoted from the former main Fig. 7A.

fig_s9 <- function(out_dir) {

  build_fig_s10_threshold_sensitivity <- function(root = ROOT) {
    grid_file <- file.path(root,
      "results/test_real/rca_sensitivity_grid/grid_results.tsv")
    if (!file.exists(grid_file)) stop("Missing grid_results.tsv")
    df_grid <- read_tsv(grid_file, show_col_types = FALSE)

    top_tools <- c("geNomad", "ViraLM", "VirSorter2")
    recall_cols <- paste0(top_tools, "_recall")
    recall_cols <- recall_cols[recall_cols %in% colnames(df_grid)]
    if (length(recall_cols) == 0) stop("No recall columns found")

    df_plot <- df_grid %>%
      select(min_breadth_pct, min_read_pairs, all_of(recall_cols)) %>%
      pivot_longer(cols = all_of(recall_cols), names_to = "tool",
                   values_to = "recall") %>%
      mutate(tool = sub("_recall$", "", tool))

    ggplot(df_plot, aes(x = factor(min_breadth_pct),
                        y = factor(min_read_pairs),
                        fill = recall)) +
      geom_tile(colour = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.0f%%", recall * 100)), size = 3) +
      facet_wrap(~tool, ncol = 3) +
      scale_fill_gradient2(low = "#d73027", mid = "#ffffbf", high = "#1a9850",
                           midpoint = 0.5, limits = c(0, 1), name = "Recall") +
      labs(x = "Min. breadth (%)", y = "Min. read pairs") +
      theme_pub()
  }

  p <- build_fig_s10_threshold_sensitivity(ROOT)
  save_fig(p, file.path(out_dir, "fig_s9.png"), w = 14, h = 5)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. S10  - Track C enrichment-ratio (R) threshold sensitivity  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# How tool MCC and MCC-based ranking shift across a 3 x 3 grid of R-cutoff
# combinations used to define TP and TN contigs:
#   TP thresholds:  R > 5, R > 10 (default), R > 20
#   TN bands:       0.33-3.0, 0.5-2.0 (default), 0.8-1.25
# A: heatmap -- rows = tools (ordered by default MCC), columns = threshold
#    combinations, cells = MCC. Star marks the default combination.
# B: ranking stability -- tool rank by MCC across combinations (1 = best),
#    one line per tool that appears in the top 5 at any combination.
# Input: table_s18_track_c_R_sensitivity_mcc.tsv in TABLES_DIR.

fig_s10 <- function(out_dir) {

  build_fig_s11_r_threshold_sensitivity <- function(root = ROOT) {
    mcc_path <- file.path(TABLES_DIR,
      "table_s18_track_c_R_sensitivity_mcc.tsv")
    if (!file.exists(mcc_path)) {
      stop("Missing table_s18_track_c_R_sensitivity_mcc.tsv")
    }
    df <- read_tsv(mcc_path, show_col_types = FALSE)

    # Pretty combination label:  "R>10 / TN 0.5-2.0"
    df <- df %>%
      mutate(
        combo_label = sprintf("R>%g\n%g-%g", r_tp_min, r_bg_low, r_bg_high),
        is_default  = as.logical(is_default)
      )

    # Annotate default label with a star
    default_label <- df$combo_label[df$is_default]
    df$combo_label[df$is_default] <- paste0(default_label, " *")

    # Preserve grid order (sort by TP then TN low)
    combo_levels <- df %>%
      arrange(r_tp_min, r_bg_low) %>%
      pull(combo_label)
    df$combo_label <- factor(df$combo_label, levels = combo_levels)

    # Long form: tool × combo → MCC
    tool_cols_in_file <- intersect(TOOL_ORDER, names(df))
    long <- df %>%
      pivot_longer(all_of(tool_cols_in_file),
                   names_to = "tool", values_to = "MCC") %>%
      mutate(MCC = suppressWarnings(as.numeric(MCC)))

    # Order tools by default-combo MCC (descending)
    default_row <- df %>% filter(is_default)
    default_mcc <- vapply(tool_cols_in_file,
                          function(t) as.numeric(default_row[[t]]),
                          numeric(1))
    tool_order_by_default <- names(sort(default_mcc, decreasing = TRUE))
    long$tool <- factor(long$tool, levels = rev(tool_order_by_default))

    # ---- Panel A: Heatmap --------------------------------------------------
    p_a <- ggplot(long, aes(x = combo_label, y = tool, fill = MCC)) +
      geom_tile(colour = "white", linewidth = 0.6) +
      geom_text(aes(label = sprintf("%.2f", MCC)), size = 2.8, colour = "black") +
      scale_fill_gradient2(low = "#f7f7f7", mid = "#fdae61", high = "#b2182b",
                           midpoint = 0.35,
                           limits = c(-0.1, max(long$MCC, na.rm = TRUE) + 0.05),
                           name = "MCC") +
      labs(x = "TP threshold (top line) / TN band (bottom line)  [* = default]",
           y = NULL, tag = "A") +
      theme_pub() +
      theme(axis.text.x = element_text(size = 8),
            legend.position = "right")

    p_a <- colour_axis_by_family(p_a, axis = "y")

    # ---- Panel B: Ranking stability ---------------------------------------
    ranks <- long %>%
      group_by(combo_label) %>%
      mutate(rank = rank(-MCC, ties.method = "min")) %>%
      ungroup()

    top_tools <- ranks %>%
      group_by(tool) %>%
      summarise(best_rank = min(rank, na.rm = TRUE), .groups = "drop") %>%
      filter(best_rank <= 5) %>%
      pull(tool) %>%
      as.character()

    ranks_top <- ranks %>% filter(tool %in% top_tools)
    ranks_top$tool <- factor(ranks_top$tool, levels = top_tools)

    p_b <- ggplot(ranks_top, aes(x = combo_label, y = rank,
                                 group = tool, colour = tool)) +
      geom_line(linewidth = 1.0, alpha = 0.85) +
      geom_point(size = 2.6) +
      scale_y_reverse(breaks = 1:max(ranks_top$rank, na.rm = TRUE)) +
      scale_colour_manual(values = TOOL_COLORS[top_tools], name = "Tool") +
      labs(x = "Threshold combination",
           y = "Rank (1 = best MCC)",
           tag = "B") +
      theme_pub() +
      theme(axis.text.x = element_text(size = 8))

    default_idx <- which(levels(ranks_top$combo_label) ==
                         as.character(df$combo_label[df$is_default]))
    if (length(default_idx) == 1) {
      p_b <- p_b + geom_vline(xintercept = default_idx,
                              linetype = "dashed", colour = "grey50",
                              linewidth = 0.5)
    }

    p_a / p_b + plot_layout(heights = c(1.2, 1))
  }

  p <- build_fig_s11_r_threshold_sensitivity(ROOT)
  save_fig(p, file.path(out_dir, "fig_s10.png"), w = 11, h = 9)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. S12  - Track B CST contrast, all 14 tools                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta MCC with 95% CI for every tool (the Fig. 6B builder, all tools).

fig_s12 <- function(out_dir) {

  build_fig_s13_trackb_all_tools <- function(root = ROOT) {
    all_tools <- intersect(TOOL_ORDER, unique(.load_trackB_deltas(root)$tool))
    build_fig5d_delta_ci(root, tools = all_tools) +
      labs(title = "Track B CST contrast, all 14 tools")
  }

  p <- build_fig_s13_trackb_all_tools(ROOT)
  save_fig(p, file.path(out_dir, "fig_s12.png"), w = 9, h = 8)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. S13  - MiTCH per-sample MCC matrix                         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# All 14 tools x 30 samples.

fig_s13 <- function(out_dir) {

  build_fig_s14_mitch_per_sample <- function(root = ROOT) {
    ps <- .load_scaleup(root, "scaleup30_per_sample_MCC.tsv")
    long <- ps %>%
      pivot_longer(-sample, names_to = "tool", values_to = "MCC") %>%
      mutate(tool = factor(tool, levels = TOOL_ORDER))
    ggplot(long, aes(x = tool, y = sample, fill = MCC)) +
      geom_tile(colour = "white", linewidth = 0.2) +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                           midpoint = 0, limits = c(-1, 1)) +
      labs(x = NULL, y = NULL, fill = "MCC") +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = element_text(size = 6))
  }

  p <- build_fig_s14_mitch_per_sample(ROOT)
  save_fig(p, file.path(out_dir, "fig_s13.png"), w = 12, h = 10)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Figs. S11, S15, S16  - drawn in Python, not in this script      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run from the repository root. <TABLES_DIR> and <FIGURES_DIR> stand for the
# two path constants defined at the top of this file.
#
#   Fig. S11:
#     python scripts/10_vmgc_crosscheck/05_plot_vmgc_overlap.py \
#         --overlap-table <TABLES_DIR>/table_s24_vmgc_benchmark_overlap.tsv \
#         --out <FIGURES_DIR>/fig_s11.png
#   Fig. S15:
#     python scripts/07_evaluation/viralm_investigation/06_plot_gt_construction.py
#   Fig. S16:
#     python scripts/07_evaluation/viralm_investigation/07_tier_evidence_breakdown.py
#
# Fig. S14 is the preserved authored ground-truth construction schematic;
# 06_plot_gt_construction.py --with-flow reproduces its numeric content.



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. S17  - Phage vs eukaryote-infecting virus recall           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mean recall per tool, faceted by fragment length. Rendered at print size.

fig_s17 <- function(out_dir) {

  build_fig3c_phage_vs_eukaryotic <- function(root = ROOT) {
    df_cat <- .load_trackA_bycategory(root)
    no_output <- .load_trackA_overall(root) %>%
      filter(!is.na(note), note == "no_output") %>%
      select(tool, length_bin)

    # A run returning no classification is unavailable, rather than measured
    # zero recall. Keep the same per-category means for every evaluable run.
    bias <- df_cat %>%
      anti_join(no_output, by = c("tool", "length_bin")) %>%
      mutate(virus_group = case_when(
        category %in% PHAGE_CATS      ~ "Phage",
        category %in% EUKARYOTIC_CATS ~ "Eukaryotic",
        TRUE                          ~ NA_character_
      )) %>%
      filter(!is.na(virus_group)) %>%
      group_by(tool, length_bin, virus_group) %>%
      summarise(mean_recall = mean(recall, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = virus_group, values_from = mean_recall)

    length_labels <- setNames(
      paste0(scales::comma(as.integer(sub("L", "", ALL_LENGTHS))), " bp"),
      ALL_LENGTHS
    )

    ggplot(bias, aes(x = Phage, y = Eukaryotic, colour = tool)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  colour = "grey65", linewidth = 0.4) +
      geom_point(size = 2.2) +
      facet_wrap(~length_bin, ncol = 3,
                 labeller = as_labeller(length_labels)) +
      scale_colour_manual(values = TOOL_COLORS, name = NULL, drop = FALSE) +
      scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
      scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
      coord_fixed() +
      guides(colour = guide_legend(ncol = 4, byrow = TRUE,
                                   override.aes = list(size = 2.5))) +
      labs(x = "Mean phage recall", y = "Mean eukaryote-infecting virus recall") +
      theme_pub() + theme_print() +
      theme(legend.position = "bottom",
            legend.key.width = unit(0.35, "cm"),
            legend.key.height = unit(0.35, "cm"),
            panel.spacing = unit(0.30, "cm"),
            plot.margin = margin(6, 8, 6, 6))
  }

  p <- build_fig3c_phage_vs_eukaryotic(ROOT)
  save_fig(p, file.path(out_dir, "fig_s17.png"), w = 6.3, h = 5.6, scale = 1)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. S18  - Track B viral contig length distribution            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fig_s18 <- function(out_dir) {

  build_fig5b_length_distribution <- function(root = ROOT) {
    df_gt <- .load_trackB_assemblies(root)
    if (nrow(df_gt) == 0) stop("No Track B assembly ground truth files found")
    viral_len <- df_gt %>%
      filter(label == "viral", length >= 500)

    ggplot(viral_len, aes(x = factor(coverage), y = length, fill = CST)) +
      geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3, linewidth = 0.4) +
      geom_hline(yintercept = 1500, linetype = "dashed",
                 colour = "grey50", linewidth = 0.5) +
      scale_y_log10(labels = scales::comma_format()) +
      scale_fill_manual(values = CST_COLORS) +
      labs(x = "Spike-in coverage (x)",
           y = "Viral contig length (bp)",
           fill = "Background") +
      theme_pub()
  }

  p <- build_fig5b_length_distribution(ROOT)
  save_fig(p, file.path(out_dir, "fig_s18.png"), w = 10, h = 5)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. S19  - Track B per-background MCC by CST arm               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Per-background MCC dots, the two arms dodged into separate sub-rows per
# tool, with a bold arm-mean bar.

fig_s19 <- function(out_dir) {

  .load_trackB_per_background <- function(root) {
    f <- file.path(root, "results/expansion/track_b/trackB_expanded_per_background_excl.tsv")
    if (!file.exists(f)) stop("Missing trackB_expanded_per_background_excl.tsv")
    manifest <- .load_trackB_manifest(root)
    pb <- read_tsv(f, show_col_types = FALSE) %>%
      inner_join(manifest, by = c("background", "stratum"))
    expected <- tidyr::crossing(background = manifest$background, tool = TOOL_ORDER)
    if (nrow(anti_join(expected, pb, by = c("background", "tool"))) > 0 ||
        anyDuplicated(pb[c("background", "tool")])) {
      stop("Track B figure requires exactly one result per retained background and tool")
    }
    pb %>% mutate(arm = stratum)
  }

  build_fig5c_cst_comparison <- function(root = ROOT,
                                         tools = .trackB_main_tools(root)) {
    pb <- .load_trackB_per_background(root) %>%
      filter(tool %in% tools) %>%
      mutate(tool = factor(tool, levels = rev(tools)),
             arm  = factor(arm, levels = c("CST-I", "CST-IV-B")))
    arm_mean <- pb %>% group_by(tool, arm) %>%
      summarise(mean_MCC = mean(MCC), .groups = "drop")
    dw <- 0.75

    p <- ggplot(pb, aes(x = MCC, y = tool, colour = arm, group = arm)) +
      geom_hline(yintercept = seq(0.5, length(tools) + 0.5, by = 1),
                 colour = "grey90", linewidth = 0.3) +
      # Jitterdodge flips its axes for this horizontal layout: width separates
      # points along the tool rows, while height must be zero to preserve MCC.
      geom_point(position = position_jitterdodge(jitter.height = 0, jitter.width = 0.09,
                                                 dodge.width = dw, seed = 1),
                 size = 2, alpha = 0.7) +
      geom_point(data = arm_mean,
                 aes(x = mean_MCC, y = tool, colour = arm, group = arm),
                 position = position_dodge(width = dw),
                 shape = 124, size = 11, show.legend = FALSE) +
      scale_colour_manual(values = .ARM_COLORS, name = "Background") +
      scale_x_continuous(limits = c(NA, 1)) +
      labs(x = "MCC at 10x (per background)", y = NULL,
           subtitle = .trackB_sample_sizes(root)) +
      theme_pub()

    p
  }

  p <- build_fig5c_cst_comparison(ROOT)
  save_fig(p, file.path(out_dir, "fig_s19.png"), w = 12, h = 7)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. S20  - Ground-truth tier distribution per sample           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fig_s20 <- function(out_dir) {

  build_fig6b_tier_by_sample <- function(root = ROOT) {
    df_gt <- .load_gt(root)
    # Benchmark-ready set contains only tier 0/1/2 (tier 3 excluded by construction);
    # panel B shows the high-confidence viral positives (Tier 1 + Tier 2) per sample.
    tier_counts <- df_gt %>%
      filter(tier %in% c(1, 2)) %>%
      count(sample, CST, tier) %>%
      mutate(tier = factor(tier))
    cst_rank <- c("CST-I" = 1, "CST-III" = 2, "CST-IV" = 3)
    sample_order <- tier_counts %>%
      distinct(sample, CST) %>%
      arrange(cst_rank[CST], sample) %>%
      pull(sample)
    tier_counts <- tier_counts %>%
      mutate(sample = factor(sample, levels = sample_order))

    ggplot(tier_counts, aes(x = sample, y = n, fill = tier)) +
      geom_col(colour = "white", linewidth = 0.3) +
      scale_fill_manual(values = TIER_COLORS, name = "Tier") +
      labs(x = NULL, y = "Viral contigs") +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  }

  p <- build_fig6b_tier_by_sample(ROOT)
  save_fig(p, file.path(out_dir, "fig_s20.png"), w = 10, h = 6)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. S21  - Viral category composition of TP contigs            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fig_s21 <- function(out_dir) {

  build_fig6d_category_composition <- function(root = ROOT) {
    df_gt <- .load_gt(root)
    rca_df <- .load_rca(root)
    gt_viral <- df_gt %>% filter(tier %in% c(1, 2))
    gt_cats  <- gt_viral %>% count(category, name = "count_gt")
    rca_tp <- rca_df %>% filter(label == "TP")
    rca_cats <- rca_tp %>%
      mutate(category = case_when(
        provirus == "Yes" ~ "prophage",
        viral_genes > 0   ~ "free_phage",
        TRUE              ~ "other"
      )) %>%
      count(category, name = "count_rca")

    cats_plot <- c("free_phage", "prophage", "eukaryotic_virus")
    cat_labels <- c("Free\nphage", "Prophage", "Eukaryote-infecting\nvirus")
    comp <- tibble(category = cats_plot) %>%
      left_join(gt_cats, by = "category") %>%
      left_join(rca_cats, by = "category") %>%
      replace_na(list(count_gt = 0, count_rca = 0)) %>%
      pivot_longer(cols = starts_with("count_"), names_to = "source",
                   values_to = "count") %>%
      mutate(
        source = ifelse(source == "count_gt", "Multi-evidence", "Track C (RCA)"),
        category = factor(category, levels = cats_plot)
      )

    ggplot(comp, aes(x = category, y = log10(count + 1), fill = source)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6,
               colour = "black", linewidth = 0.3) +
      scale_x_discrete(labels = cat_labels) +
      scale_fill_manual(values = c("Multi-evidence" = "#3182bd",
                                   "Track C (RCA)" = "#e6550d"),
                        name = NULL) +
      labs(x = NULL, y = expression(log[10](count + 1))) +
      theme_pub()
  }

  p <- build_fig6d_category_composition(ROOT)
  save_fig(p, file.path(out_dir, "fig_s21.png"), w = 8, h = 6)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig. S22  - MiTCH best-performing tools by CST                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Top-3 tools per well-populated CST, parsed from "tool:value" cells.

fig_s22 <- function(out_dir) {

  build_fig8b_per_cst <- function(root = ROOT) {
    pc <- .load_scaleup(root, "scaleup30_per_CST.tsv") %>%
      filter(CST %in% c("IV-B", "I", "III"))
    long <- pc %>%
      pivot_longer(c(top1, top2, top3), names_to = "rank", values_to = "entry") %>%
      separate(entry, into = c("tool", "mcc"), sep = ":", convert = TRUE) %>%
      mutate(
        cst_lab = sprintf("CST-%s (n=%d)", CST, n_samples),
        cst_lab = factor(cst_lab,
                         levels = c(sprintf("CST-IV-B (n=%d)", pc$n_samples[pc$CST == "IV-B"]),
                                    sprintf("CST-I (n=%d)",   pc$n_samples[pc$CST == "I"]),
                                    sprintf("CST-III (n=%d)", pc$n_samples[pc$CST == "III"]))),
        rank = factor(rank, levels = c("top1", "top2", "top3"),
                      labels = c("1st", "2nd", "3rd"))
      )

    ggplot(long, aes(x = rank, y = mcc, fill = tool)) +
      geom_col(width = 0.7, colour = "black", linewidth = 0.3) +
      geom_text(aes(label = tool), vjust = -0.4, size = 2.7) +
      facet_wrap(~cst_lab) +
      scale_fill_manual(values = TOOL_COLORS, guide = "none") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
      labs(x = "Within-CST rank", y = "MCC") +
      theme_pub()
  }

  p <- build_fig8b_per_cst(ROOT)
  save_fig(p, file.path(out_dir, "fig_s22.png"), w = 10, h = 6)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Run  - render the selected figures                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Each figure runs in its own function scope inside tryCatch: a failure is
# reported and the remaining figures still render. Exit status is 1 if any
# figure failed.

FIGURES <- list(
  fig2    = fig2,    fig3    = fig3,    fig4    = fig4,    fig5    = fig5,
  fig6    = fig6,    fig7    = fig7,    fig8    = fig8,    fig9    = fig9,
  fig_s1  = fig_s1,  fig_s2  = fig_s2,  fig_s3  = fig_s3,  fig_s4  = fig_s4,
  fig_s5  = fig_s5,  fig_s6  = fig_s6,  fig_s7  = fig_s7,  fig_s8  = fig_s8,
  fig_s9  = fig_s9,  fig_s10 = fig_s10, fig_s12 = fig_s12, fig_s13 = fig_s13,
  fig_s17 = fig_s17, fig_s18 = fig_s18, fig_s19 = fig_s19, fig_s20 = fig_s20,
  fig_s21 = fig_s21, fig_s22 = fig_s22
)

targets <- names(FIGURES)
if (!is.null(ONLY)) {
  unknown <- setdiff(ONLY, targets)
  if (length(unknown) > 0) {
    stop("Unknown --only target(s): ", paste(unknown, collapse = ", "),
         "\nAvailable: ", paste(targets, collapse = ", "))
  }
  targets <- intersect(targets, ONLY)   # keep render order
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
message(sprintf("[figures] root=%s", ROOT))
message(sprintf("[figures] out =%s", OUT_DIR))

n_fail <- 0L
for (id in targets) {
  t0 <- Sys.time()
  res <- tryCatch({
    out_file <- FIGURES[[id]](OUT_DIR)
    list(status = "ok",
         note = sprintf("%5.1fs  %s",
                        as.numeric(difftime(Sys.time(), t0, units = "secs")),
                        basename(out_file)))
  }, error = function(e) {
    list(status = "fail", note = conditionMessage(e))
  })
  if (res$status != "ok") n_fail <- n_fail + 1L
  message(sprintf("[figures] %-8s %-5s  %s", id, res$status, res$note))
}
message(sprintf("[figures] %d ok, %d failed", length(targets) - n_fail, n_fail))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Session info                                                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print(sessionInfo())
if (n_fail > 0 && !interactive()) quit(save = "no", status = 1)
