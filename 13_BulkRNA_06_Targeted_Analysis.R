# =============================================================================
# 13_BulkRNA_06_Targeted_Analysis.R
# Targeted IL1-family visualization and bulk-derived signature export
# =============================================================================

library(openxlsx)
library(pheatmap)
library(gridExtra)
library(dplyr)

setwd("D:\\Eve_Bulk\\R\\Pairwise")


# =============================================================================
# 1. IL1-family heatmap across P0/P1 and NGFR states
# =============================================================================
#
# This section reproduces the targeted bulk-RNA IL1-family comparison from the
# original workflow. It requires "P0+P1 N- vs N+.txt" if that matrix has already
# been generated in the project.

if (file.exists("P0+P1 N- vs N+.txt")) {

  cts <- read.table(
    "P0+P1 N- vs N+.txt",
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  counts <- as.matrix(cts)

  genes <- c(
    "IL1A", "IL1B", "IL1RN",
    "IL1R1", "IL1R2", "IL1RAP"
  )

  genes <- intersect(genes, rownames(counts))
  sub_counts <- counts[genes, , drop = FALSE]

  zmat <- t(scale(t(sub_counts)))
  zmat[zmat > 2] <- 2
  zmat[zmat < -2] <- -2

  grp_P0_neg <- grep("_P0_NGFR_neg$", colnames(zmat))
  grp_P0_pos <- grep("_P0_NGFR_pos$", colnames(zmat))
  grp_P1_neg <- grep("_P1_Ctrl_NGFR_neg$", colnames(zmat))
  grp_P1_pos <- grep("_P1_Ctrl_NGFR_pos$", colnames(zmat))

  avg_z <- cbind(
    P0_neg = rowMeans(zmat[, grp_P0_neg, drop = FALSE], na.rm = TRUE),
    P0_pos = rowMeans(zmat[, grp_P0_pos, drop = FALSE], na.rm = TRUE),
    P1_neg = rowMeans(zmat[, grp_P1_neg, drop = FALSE], na.rm = TRUE),
    P1_pos = rowMeans(zmat[, grp_P1_pos, drop = FALSE], na.rm = TRUE)
  )

  my_cols <- colorRampPalette(
    c("#3D6B98", "#F0F0F0", "#C0203C")
  )(100)

  my_breaks <- seq(
    -1.5,
    1.5,
    length.out = length(my_cols) + 1
  )

  pdf("IL1_Heatmap.pdf", width = 3.5, height = 4)

  pheatmap(
    avg_z,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = my_cols,
    breaks = my_breaks,
    border_color = NA,
    show_rownames = TRUE,
    show_colnames = TRUE,
    main = "Average z-score (IL1 family)"
  )

  dev.off()


  # Separate pairwise panels used in the original workflow.
  hm1 <- pheatmap(
    avg_z[, c("P0_neg", "P0_pos"), drop = FALSE],
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = my_cols,
    breaks = my_breaks,
    legend = FALSE,
    silent = TRUE
  )

  hm2 <- pheatmap(
    avg_z[, c("P1_neg", "P1_pos"), drop = FALSE],
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = my_cols,
    breaks = my_breaks,
    legend = FALSE,
    silent = TRUE
  )

  hm3 <- pheatmap(
    avg_z[, c("P0_neg", "P1_neg"), drop = FALSE],
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = my_cols,
    breaks = my_breaks,
    legend = FALSE,
    silent = TRUE
  )

  hm4 <- pheatmap(
    avg_z[, c("P0_pos", "P1_pos"), drop = FALSE],
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = my_cols,
    breaks = my_breaks,
    legend = TRUE,
    silent = TRUE
  )

  pdf("IL1_family_pairwise_heatmaps.pdf", width = 8, height = 4)

  grid.arrange(
    hm1$gtable,
    hm2$gtable,
    hm3$gtable,
    hm4$gtable,
    nrow = 1
  )

  dev.off()
}


# =============================================================================
# 2. Export the top 50 bulk-derived signatures used in Figure 2I and Figure 3I
# =============================================================================
#
# Figure 2I:
#   P0 tissue-derived NGFR-positive and NGFR-negative signatures
#
# Figure 3I:
#   P1 control organoid-derived NGFR-positive and NGFR-negative signatures
#
# The actual projection of these signatures onto the scRNA-seq UMAP is kept in
# Fig2.R and 05_Figure3.R. This file only performs the bulk-RNA-side extraction.


extract_directional_top50 <- function(
    result_file,
    positive_output,
    negative_output
) {

  res <- read.table(
    result_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  res <- res %>%
    filter(
      !is.na(log2FoldChange),
      !is.na(pvalue)
    )

  positive <- res %>%
    filter(log2FoldChange > 0) %>%
    arrange(desc(log2FoldChange), pvalue) %>%
    slice_head(n = 50)

  negative <- res %>%
    filter(log2FoldChange < 0) %>%
    arrange(log2FoldChange, pvalue) %>%
    slice_head(n = 50)

  write.xlsx(
    positive,
    positive_output,
    rowNames = FALSE
  )

  write.xlsx(
    negative,
    negative_output,
    rowNames = FALSE
  )

  invisible(list(
    positive = positive,
    negative = negative
  ))
}


# Figure 2I source comparison: P0 NGFR+ vs NGFR-
Fig2I_signatures <- extract_directional_top50(
  "{Pairwise}[DESeq2]P0_N+_vs_N-.txt",
  "[Signature Top50]P0_NGFR_pos.xlsx",
  "[Signature Top50]P0_NGFR_neg.xlsx"
)


# Figure 3I source comparison: P1 control NGFR+ vs NGFR-
Fig3I_signatures <- extract_directional_top50(
  "{Pairwise}[DESeq2]P1_Ctrl N- vs N+.txt",
  "[Signature Top50]P1_NGFR_pos.xlsx",
  "[Signature Top50]P1_NGFR_neg.xlsx"
)
