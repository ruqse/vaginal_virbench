# 15_plot_dark_matter_gene_maps.R — gene maps + evidence grid for the 48
# Track C dark-matter contigs.
#
# Inputs (written by 13_* and 14_* in this directory):
#   results/vmgc_crosscheck/dark_matter_metavr_refseq/gene_map_genes.tsv
#   results/vmgc_crosscheck/dark_matter_metavr_refseq/gene_map_tool_calls.tsv
#   results/vmgc_crosscheck/dark_matter_metavr_refseq/gene_map_topology.tsv
#   results/vmgc_crosscheck/dark_matter_metavr_refseq/dark_matter_integrated_evidence.tsv
# Output:
#   docs/figures/fig_dark_matter_gene_maps.{png,pdf}
#
# Run (Rscript is broken in farukRenv; use R):
#   module load Miniconda3/25.7.0-2 && conda activate farukRenv
#   export FONTCONFIG_PATH=$CONDA_PREFIX/etc/fonts
#   R --no-echo --no-save --file=scripts/10_vmgc_crosscheck/15_plot_dark_matter_gene_maps.R

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(patchwork)
})

PROJ <- Sys.getenv("VBENCH_ROOT")
if (!nzchar(PROJ)) stop("set VBENCH_ROOT (see config/paths.example.sh)")
D    <- file.path(PROJ, "results/vmgc_crosscheck/dark_matter_metavr_refseq")
OUTD <- file.path(PROJ, "docs/figures")

# Fixed tool display order, shared with all benchmark figures (pub_style.R).
TOOL_ORDER <- c("geNomad", "ViraLM", "TransGINmer", "VirSorter2", "Jaeger",
                "HVSeeker", "PPR-Meta", "VIBRANT", "DeepVirFinder", "Seeker",
                "VirSorter", "MetaPhinder", "VirFinder", "Sourmash")

# Categorical slots 1-3 of the reference palette (validated all-pairs, light);
# "no domain" is a neutral, not a categorical hue.
GENE_COL <- c("Plasmid replication / mobilisation" = "#2a78d6",
              "Phage-type domain"                  = "#eb6834",
              "Other Pfam domain"                  = "#1baf7a",
              "No Pfam domain"                     = "#d9d8d2")
EV_COL <- c("Strong" = "#52514e", "Partial or weak" = "#b9b8b2", "None" = "#fcfcfb")
INK  <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e4e3dd"

genes <- read_tsv(file.path(D, "gene_map_genes.tsv"), show_col_types = FALSE)
tools <- read_tsv(file.path(D, "gene_map_tool_calls.tsv"), show_col_types = FALSE)
topo  <- read_tsv(file.path(D, "gene_map_topology.tsv"), show_col_types = FALSE)
ev    <- read_tsv(file.path(D, "dark_matter_integrated_evidence.tsv"), show_col_types = FALSE)

# --- Row order: evidence group, then length (longest first) -----------------
ev <- ev |>
  mutate(
    group = case_when(
      pfam_plasmid_rep_mob == 1 ~ "Plasmid Rep/Mob domain",
      genomad_plasmid_call == 1 | genomad_plasmid_score >= 0.7 ~ "geNomad plasmid evidence only",
      TRUE ~ "No plasmid evidence"),
    group = factor(group, levels = c("Plasmid Rep/Mob domain",
                                     "geNomad plasmid evidence only",
                                     "No plasmid evidence")),
    label = sub("_length_.*$", "", contig_id) |> sub(pattern = "__", replacement = " ")
  ) |>
  arrange(group, desc(length_bp))

# Numeric y with a header row above each group (top of plot = first row).
rows <- list(); y <- 0; headers <- list()
for (g in levels(ev$group)) {
  sub_ev <- filter(ev, group == g)
  y <- y + 1.6
  headers[[g]] <- tibble(group = g, y = y, n = nrow(sub_ev))
  for (i in seq_len(nrow(sub_ev))) { y <- y + 1; rows[[length(rows) + 1]] <- tibble(contig_id = sub_ev$contig_id[i], y = y) }
}
rows <- bind_rows(rows); headers <- bind_rows(headers)
ev <- left_join(ev, rows, by = "contig_id")
YMAX <- max(rows$y) + 0.7

# --- Panel A: gene arrows ------------------------------------------------------
H <- 0.36  # half-height of an arrow head
arrow_poly <- function(id, s, e, strand, yc) {
  len  <- e - s + 1
  head <- min(160, 0.45 * len)
  if (strand >= 0) {
    x <- c(s, e - head, e - head, e, e - head, e - head, s)
  } else {
    x <- c(e, s + head, s + head, s, s + head, s + head, e)
  }
  yy <- yc + c(-H * 0.72, -H * 0.72, -H, 0, H, H * 0.72, H * 0.72)
  tibble(gene = id, x = x, y = yy)
}
genes <- genes |> left_join(rows, by = "contig_id") |>
  mutate(category = factor(category, levels = names(GENE_COL)),
         top_family = sub(";.*$", "", pfam_families))
polys <- genes |> rowwise() |>
  do(arrow_poly(.$gene, .$start, .$end, .$strand, .$y)) |> ungroup() |>
  left_join(select(genes, gene, category), by = "gene")
glabs <- genes |> filter(!is.na(pfam_families), pfam_families != "") |>
  mutate(xm = (start + end) / 2, inside = (end - start + 1) >= 110 * nchar(top_family))

pA <- ggplot() +
  geom_segment(data = ev, aes(x = 1, xend = length_bp, y = y, yend = y),
               colour = "#8d8c86", linewidth = 0.35) +
  geom_polygon(data = polys, aes(x, y, group = gene, fill = category),
               colour = "#fcfcfb", linewidth = 0.25) +
  geom_text(data = filter(glabs, inside), aes(xm, y - 0.02, label = top_family),
            size = 1.6, colour = "white", fontface = "bold") +
  geom_text(data = filter(glabs, !inside), aes(xm, y - H - 0.04, label = top_family),
            vjust = 0, size = 1.7, colour = INK2) +
  geom_text(data = headers, aes(x = -40, y = y, label = sprintf("%s (n = %d)", group, n)),
            hjust = 0, fontface = "bold", size = 2.6, colour = INK) +
  scale_fill_manual(values = GENE_COL, name = "Gene (Pfam 35.0)", drop = FALSE) +
  scale_x_continuous(name = "Position (bp)", breaks = seq(0, 7000, 1000),
                     labels = scales::comma, expand = expansion(mult = c(0, 0.01)),
                     limits = c(-60, 7000)) +
  scale_y_reverse(breaks = ev$y, labels = ev$label, expand = c(0, 0),
                  limits = c(YMAX, 0.6)) +
  theme_minimal(base_size = 8) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = GRID, linewidth = 0.25),
        axis.text.y = element_text(size = 6, colour = INK, family = "mono"),
        axis.title.y = element_blank(), axis.text.x = element_text(colour = INK2),
        legend.position = "bottom", legend.title = element_text(face = "bold"),
        plot.background = element_rect(fill = "white", colour = NA)) +
  guides(fill = guide_legend(nrow = 2))

# --- Panel B: evidence grid ----------------------------------------------------
lev <- function(strong, partial) factor(ifelse(strong, "Strong", ifelse(partial, "Partial or weak", "None")),
                                         levels = names(EV_COL))
dbev <- ev |>
  left_join(topo, by = "contig_id") |>
  transmute(contig_id, y,
            PLSDB    = lev(plsdb_status == "plsdb_near_identical", plsdb_status == "plsdb_partial"),
            MetaVR   = lev(metavr_status == "metavr_species", metavr_status == "metavr_relative"),
            RefSeq   = lev(refseq_status == "refseq_homolog", refseq_status == "refseq_weak"),
            `geNomad plasmid` = lev(genomad_plasmid_call == 1, genomad_plasmid_score >= 0.7),
            DTR      = lev(!is.na(topology) & topology == "DTR", FALSE)) |>
  pivot_longer(-c(contig_id, y), names_to = "column", values_to = "level") |>
  mutate(block = "Evidence")
tev <- tools |> inner_join(select(ev, contig_id, y), by = "contig_id") |>
  transmute(contig_id, y, column = tool, level = lev(call == 1, FALSE), block = "Viral call (default threshold)")
grid <- bind_rows(dbev, tev) |>
  mutate(column = factor(column, levels = c("PLSDB", "MetaVR", "RefSeq", "geNomad plasmid", "DTR", TOOL_ORDER)),
         block = factor(block, levels = c("Evidence", "Viral call (default threshold)")))
counts <- tev |> group_by(column) |> summarise(n = sum(level == "Strong"), .groups = "drop") |>
  mutate(column = factor(column, levels = levels(grid$column)),
         block = factor("Viral call (default threshold)", levels = levels(grid$block)))

pB <- ggplot(grid, aes(column, y)) +
  geom_tile(aes(fill = level), width = 0.86, height = 0.80, colour = "#c9c8c2", linewidth = 0.2) +
  geom_text(data = counts, aes(column, YMAX + 0.05, label = n), vjust = 0, size = 1.8, colour = INK2) +
  facet_grid(. ~ block, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = EV_COL, name = "Evidence", drop = FALSE) +
  scale_x_discrete(position = "top") +
  scale_y_reverse(expand = c(0, 0), limits = c(YMAX + 0.9, 0.6)) +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, colour = INK, size = 6.5),
        strip.text = element_text(face = "bold", size = 7, colour = INK),
        strip.placement = "outside", legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        plot.background = element_rect(fill = "white", colour = NA))

# Keep both panels on the same y range (panel B adds a counts row at the bottom).
pA <- pA + scale_y_reverse(breaks = ev$y, labels = ev$label, expand = c(0, 0), limits = c(YMAX + 0.9, 0.6))

fig <- pA + pB + plot_layout(widths = c(2.3, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11))

dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUTD, "fig_dark_matter_gene_maps.png"), fig, width = 11, height = 11.5, dpi = 300, bg = "white")
ggsave(file.path(OUTD, "fig_dark_matter_gene_maps.pdf"), fig, width = 11, height = 11.5, bg = "white")
cat("wrote", file.path(OUTD, "fig_dark_matter_gene_maps.{png,pdf}"), "\n")
