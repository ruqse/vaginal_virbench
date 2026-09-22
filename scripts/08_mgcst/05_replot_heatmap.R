#!/usr/bin/env Rscript
# Regenerate mgCST heatmap with readable sample labels.
#
# The upstream run_VISTA.R hardcodes fontsize_col = 1 (designed for 53+ samples).
# With 13 samples we have space for proper labels. This script reads the existing
# VISTA relabund output and reproduces the heatmap with improved sizing.
#
# Usage:
#   module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
#   Rscript scripts/08_mgcst/05_replot_heatmap.R
#
# Input:  results/mgcst/vista/relabund_w_mgCSTs_*.csv (from VISTA pipeline)
# Output: results/mgcst/mgCST_heatmap.pdf

library(pheatmap)

proj <- Sys.getenv("VBENCH_ROOT")

# Find the date-stamped relabund file (robust to any run date)
candidates <- Sys.glob(file.path(proj, "results/mgcst/vista/relabund_w_mgCSTs_*.csv"))
if (length(candidates) == 0) stop("No relabund_w_mgCSTs_*.csv found in results/mgcst/vista/")
infile <- candidates[length(candidates)]  # latest if multiple
outfile <- file.path(proj, "results/mgcst/mgCST_heatmap.pdf")

cat("Input:", infile, "\n")

# --- Read VISTA relabund output ---
relabund <- read.csv(infile, row.names = 1, check.names = FALSE)
cat("Loaded:", nrow(relabund), "samples x", ncol(relabund), "columns\n")

# Find where taxa end and mgCST Yue-Clayton theta columns begin
mgcst1_idx <- which(colnames(relabund) == "mgCST 1")
if (length(mgcst1_idx) == 0) stop("Cannot find 'mgCST 1' column")
n <- mgcst1_idx - 1  # number of taxa columns
cat("Taxa columns: 1 to", n, "\n")

# Get assigned mgCST from last column
mgcst_col <- relabund[["mgCST"]]
cat("mgCST assignments:", paste(unique(mgcst_col), collapse = ", "), "\n")

# --- Reconstruct heatmap data (same logic as run_VISTA.R lines 404-415) ---

# Taxa relative abundances, sorted by total abundance
relabund_taxa <- relabund[, 1:n]
relabund_taxa <- relabund_taxa[, order(colSums(relabund_taxa), decreasing = TRUE)]

# Add mgCST label (numeric only)
relabund_taxa$mgCST <- gsub("mgCST ", "", mgcst_col)

# Sort rows by mgCST number
relabund_taxa <- relabund_taxa[order(as.numeric(relabund_taxa[["mgCST"]])), ]

# Replace underscores with spaces in taxa names (cosmetic)
colnames(relabund_taxa) <- gsub("_", " ", colnames(relabund_taxa))

# Top 50 most abundant taxa as heatmap matrix (taxa = rows, samples = columns)
ntaxa <- min(50, n)
mat <- t(as.matrix(relabund_taxa[, 1:ntaxa]))

# --- mgCST color palette (identical to run_VISTA.R lines 393-400) ---
mgCST_palette <- data.frame(
  mgCST = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15",
            "16","17","18","19","20","21","22","23","24","25","","NA"),
  color = c("#FE0308","#F54C5E","#F07084","#EC94A5","#F0BCCC","#F6D3DA",
            "#86C61A","#B4DB29",
            "#F68A11","#FF981C","#FFA435",
            "#FAE727","#FBEA3F","#FBED58",
            "#E1C775",
            "#589682","#6BA290",
            "#2C31A0","#3C44A8","#444DAC","#676EBC","#6B7EC0","#829CCD",
            "#C7FFC7","#8c8c8c","white","white"),
  stringsAsFactors = FALSE
)

# Annotation: mgCST per sample
annotation_col <- data.frame(mgCST = relabund_taxa[["mgCST"]])
rownames(annotation_col) <- rownames(relabund_taxa)

# Annotation colors: only mgCSTs present in our data
present <- mgCST_palette$mgCST %in% annotation_col$mgCST
annotation_colors <- list(
  mgCST = setNames(mgCST_palette$color[present], mgCST_palette$mgCST[present])
)

# Heatmap color ramp (same as run_VISTA.R line 403)
colfunc <- colorRampPalette(c("khaki", "limegreen", "darkslategray1",
                               "mediumblue", "magenta", "red"))

# --- Plot with readable labels ---
cat("Saving heatmap to:", outfile, "\n")

pdf(outfile, width = 10, height = 12)
pheatmap(
  mat,
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  color         = colfunc(100),
  main          = paste("mgCST Heatmap\nnSamples =", nrow(relabund_taxa)),
  fontsize_row  = 7,
  fontsize_col  = 8,       # was 1 in run_VISTA.R -- now readable
  angle_col     = 45,      # angled for better legibility
  cellwidth     = 30,      # fixed width per sample
  annotation_col    = annotation_col,
  annotation_colors = annotation_colors,
  legend        = TRUE,
  border_color  = "grey90"  # light grid lines
)
dev.off()

cat("Done.\n")
