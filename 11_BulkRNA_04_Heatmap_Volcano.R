# =============================================================================
# 11_BulkRNA_04_Heatmap_Volcano.R
# Marker heatmaps and DEG volcano plots
# =============================================================================

library(pheatmap)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)

setwd("D:\\Eve_Bulk\\R\\Pairwise")


# =============================================================================
# 1. Marker heatmaps
# =============================================================================

marker_genes <- c(
  "COL17A1", "IGFBP5", "TP73", "PDPN", "NGFR",
  "KRT15", "KRT14", "MKI67", "CDK1", "KRT6B",
  "KRT6C", "KRT4", "ITGA6", "TP63", "CDH1"
)

plot_marker_heatmap <- function(
    input_file,
    output_file,
    group_labels,
    group_colors,
    width = 6,
    height = 5
) {

  counts <- read.table(
    input_file,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  genes_present <- intersect(marker_genes, rownames(counts))
  mat <- as.matrix(counts[genes_present, , drop = FALSE])

  mat_z <- t(scale(t(mat)))
  mat_z[mat_z > 2] <- 2
  mat_z[mat_z < -2] <- -2

  annotation_col <- data.frame(
    condition = factor(group_labels, levels = unique(group_labels))
  )

  rownames(annotation_col) <- colnames(mat_z)

  annotation_colors <- list(
    condition = group_colors
  )

  my_palette <- colorRampPalette(
    c("#3D6B98", "#F0F0F0", "#C0203C")
  )(100)

  brks <- seq(-2, 2, length.out = length(my_palette) + 1)

  pdf(output_file, width = width, height = height)

  pheatmap(
    mat_z,
    color = my_palette,
    breaks = brks,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    border_color = NA
  )

  dev.off()
}


# P0 NGFR- vs NGFR+
p0 <- read.table(
  "P0 N- vs N+.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

plot_marker_heatmap(
  "P0 N- vs N+.txt",
  "[Heatmap]P0 N- vs N+.pdf",
  c(
    rep("NGFR_neg", sum(grepl("NGFR_neg", colnames(p0)))),
    rep("NGFR_pos", sum(grepl("NGFR_pos", colnames(p0))))
  ),
  c(
    "NGFR_neg" = "#e5ce81",
    "NGFR_pos" = "#459943"
  )
)


# P1 control NGFR- vs NGFR+
p1 <- read.table(
  "P1_Ctrl N- vs N+.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

plot_marker_heatmap(
  "P1_Ctrl N- vs N+.txt",
  "[Heatmap]P1_Ctrl N- vs N+.pdf",
  ifelse(grepl("NGFR_neg", colnames(p1)), "NGFR_neg", "NGFR_pos"),
  c(
    "NGFR_neg" = "#e5ce81",
    "NGFR_pos" = "#459943"
  )
)


# NGFR+ P0 vs P1
nplus <- read.table(
  "[Ctrl] N+ P0 vs P1.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

plot_marker_heatmap(
  "[Ctrl] N+ P0 vs P1.txt",
  "[Heatmap][Ctrl] N+ P0 vs P1.pdf",
  ifelse(grepl("_P0_", colnames(nplus)), "P0", "P1"),
  c(
    "P0" = "#e5ce81",
    "P1" = "#459943"
  )
)


# NGFR- P0 vs P1
nminus <- read.table(
  "[Ctrl] N- P0 vs P1.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

plot_marker_heatmap(
  "[Ctrl] N- P0 vs P1.txt",
  "[Heatmap][Ctrl] N- P0 vs P1.pdf",
  ifelse(grepl("_P0_", colnames(nminus)), "P0", "P1"),
  c(
    "P0" = "#e5ce81",
    "P1" = "#459943"
  )
)


# P1 all cells: control vs IL1A
all_il1a <- read.table(
  "P1 Ctr vs IL1A.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

plot_marker_heatmap(
  "P1 Ctr vs IL1A.txt",
  "[Heatmap]P1 Ctr vs IL1A.pdf",
  ifelse(grepl("Ctrl|Ctr", colnames(all_il1a)), "Ctrl", "IL1A"),
  c(
    "Ctrl" = "#e5ce81",
    "IL1A" = "#459943"
  )
)


# P1 NGFR+ control vs IL1A
npos_il1a <- read.table(
  "P1 N+ Ctr vs IL1A.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

plot_marker_heatmap(
  "P1 N+ Ctr vs IL1A.txt",
  "[Heatmap]P1 N+ Ctr vs IL1A.pdf",
  ifelse(grepl("Ctrl|Ctr", colnames(npos_il1a)), "Ctrl", "IL1A"),
  c(
    "Ctrl" = "#e5ce81",
    "IL1A" = "#459943"
  )
)


# P1 NGFR- control vs IL1A
nneg_il1a <- read.table(
  "P1 N- Ctr vs IL1A.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

plot_marker_heatmap(
  "P1 N- Ctr vs IL1A.txt",
  "[Heatmap]P1 N- Ctr vs IL1A.pdf",
  ifelse(grepl("Ctrl|Ctr", colnames(nneg_il1a)), "Ctrl", "IL1A"),
  c(
    "Ctrl" = "#e5ce81",
    "IL1A" = "#459943"
  )
)


# =============================================================================
# 2. Volcano plots
# =============================================================================
# The primary significance definition follows the original plotting code:
# |log2FoldChange| > 0.585 and adjusted P value < 0.05.

plot_volcano <- function(
    result_file,
    output_file,
    title,
    label_up = character(),
    label_down = character(),
    lfc_cutoff = 0.585,
    padj_cutoff = 0.05
) {

  res <- read.table(
    result_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  res <- res %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      direction = case_when(
        log2FoldChange > lfc_cutoff & padj < padj_cutoff ~ "Up",
        log2FoldChange < -lfc_cutoff & padj < padj_cutoff ~ "Down",
        TRUE ~ "NS"
      ),
      neglog10_padj = -log10(padj),
      label = ifelse(
        geneNames %in% c(label_up, label_down),
        geneNames,
        NA_character_
      )
    )

  cat(
    title,
    "\nUpregulated:", sum(res$direction == "Up"),
    "\nDownregulated:", sum(res$direction == "Down"),
    "\n"
  )

  p <- ggplot(
    res,
    aes(x = log2FoldChange, y = neglog10_padj, color = direction)
  ) +
    geom_point(size = 1, alpha = 0.75) +
    geom_vline(
      xintercept = c(-lfc_cutoff, lfc_cutoff),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_text_repel(
      aes(label = label),
      na.rm = TRUE,
      size = 3,
      max.overlaps = Inf
    ) +
    scale_color_manual(
      values = c(
        "Down" = "#3D6B98",
        "NS" = "gray75",
        "Up" = "#C0203C"
      )
    ) +
    labs(
      title = title,
      x = expression(log[2] * " Fold Change"),
      y = expression(-log[10] * " adjusted P value"),
      color = NULL
    ) +
    theme_classic(base_size = 13)

  ggsave(
    output_file,
    p,
    width = 7,
    height = 6
  )

  invisible(p)
}


# Selected labels are retained from the original volcano-plot code.

plot_volcano(
  "{Pairwise}[DESeq2]P0_N+_vs_N-.txt",
  "[Volcano]P0 N+ vs N-.pdf",
  "P0 NGFR+ vs NGFR-",
  label_up = c("NGFR", "PDPN", "IL1R2", "COL17A1", "P73", "KRT14", "KRT15"),
  label_down = c("KRT13", "KRT6B", "KRT6C", "KRT4")
)

plot_volcano(
  "{Pairwise}[DESeq2]P1_Ctrl N- vs N+.txt",
  "[Volcano]P1 N+ vs N-.pdf",
  "P1 Control NGFR+ vs NGFR-",
  label_up = c("NGFR", "PDPN", "IL1R2", "COL17A1", "P73", "KRT14", "KRT15"),
  label_down = c("KRT13", "KRT6B", "KRT6C", "KRT4")
)

plot_volcano(
  "{Pairwise}[DESeq2]N+ P1 vs P0.txt",
  "[Volcano]N+ P1 vs P0.pdf",
  "NGFR+ P1 vs P0",
  label_up = c("WNT7A", "DKK1", "CXCL5", "IL1B", "PTCHD4", "DUSP4", "NFKB2", "KRT14"),
  label_down = c("WIF1", "PAPPA", "IL1R2", "IGFBP5", "WNT5B", "Notch4", "TP73", "COL17A1", "KRT15")
)

plot_volcano(
  "{Pairwise}[DESeq2]N- P1 vs P0.txt",
  "[Volcano]N- P1 vs P0.pdf",
  "NGFR- P1 vs P0",
  label_up = c("IL1B", "HES5", "KRT17", "IL1A", "DLL3", "FGF1", "FGFR1", "EREG", "AREG", "KRT6B", "KRT6C"),
  label_down = c("KRT4", "KRT78", "KRT13", "AQP5")
)

plot_volcano(
  "{Pairwise}[DESeq2]P1 IL1A vs Ctr.txt",
  "[Volcano]P1 IL1A vs Ctr.pdf",
  "P1 IL1A vs Control"
)

plot_volcano(
  "{Pairwise}[DESeq2]P1 N+ IL1A vs Ctr.txt",
  "[Volcano]P1 N+ IL1A vs Ctr.pdf",
  "P1 NGFR+ IL1A vs Control",
  label_up = c("EMILIN2", "CXCL1", "IL13RA2", "CXCL3", "CXCL6", "CCL20", "CXCL8", "IL1A"),
  label_down = c("LAMB4", "LSP1", "USP35", "HSPB8")
)

plot_volcano(
  "{Pairwise}[DESeq2]P1 N- IL1A vs Ctr.txt",
  "[Volcano]P1 N- IL1A vs Ctr.pdf",
  "P1 NGFR- IL1A vs Control",
  label_up = c("SPRR2A", "SPRR2D", "DEFB4A", "S100A7", "Il1RL1", "CCL20", "CXCL1", "CXCL8", "IL1A", "TNF", "NFKBIA"),
  label_down = c("IGFL1", "LYRM9")
)
