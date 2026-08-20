# =============================================================================
# 09_BulkRNA_02_PCA_QC.R
# Pairwise PCA plots and DESeq2 library-size QC
# =============================================================================

library(ggplot2)
library(DESeq2)

setwd("D:\\Eve_Bulk\\R\\Pairwise")


# =============================================================================
# Helper: PCA from a pairwise raw-count matrix
# =============================================================================

plot_pairwise_pca <- function(
    input_file,
    output_file,
    group_function,
    group_colors,
    width = 6,
    height = 6
) {

  counts <- read.table(
    input_file,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE
  )

  # Follow the original PCA workflow: log2(count + 1), remove zero-variance
  # genes, then run PCA with gene scaling.
  log_counts <- log2(counts + 1)

  gene_vars <- apply(log_counts, 1, var)
  log_counts <- log_counts[gene_vars > 0, , drop = FALSE]

  pca_res <- prcomp(t(log_counts), scale. = TRUE)

  pca_df <- data.frame(
    Sample = rownames(pca_res$x),
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    Group = group_function(rownames(pca_res$x))
  )

  p <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
    geom_point(size = 3) +
    geom_text(aes(label = Sample), vjust = -0.5, size = 3) +
    scale_color_manual(values = group_colors) +
    labs(
      x = paste0(
        "PC1 (",
        round(summary(pca_res)$importance[2, 1] * 100, 1),
        "%)"
      ),
      y = paste0(
        "PC2 (",
        round(summary(pca_res)$importance[2, 2] * 100, 1),
        "%)"
      )
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black")
    )

  ggsave(output_file, p, width = width, height = height)

  invisible(list(pca = pca_res, plot = p))
}


# =============================================================================
# 1. P0: NGFR-negative vs NGFR-positive
# =============================================================================

plot_pairwise_pca(
  "P0 N- vs N+.txt",
  "[PCA] P0 N- vs N+.pdf",
  group_function = function(x) {
    ifelse(grepl("NGFR_neg", x), "NGFR_neg", "NGFR_pos")
  },
  group_colors = c(
    "NGFR_neg" = "navy",
    "NGFR_pos" = "firebrick3"
  )
)


# =============================================================================
# 2. P1 control: NGFR-negative vs NGFR-positive
# =============================================================================

plot_pairwise_pca(
  "P1_Ctrl N- vs N+.txt",
  "[PCA] P1_Ctrl N- vs N+.pdf",
  group_function = function(x) {
    ifelse(grepl("NGFR_neg", x), "NGFR_neg", "NGFR_pos")
  },
  group_colors = c(
    "NGFR_neg" = "navy",
    "NGFR_pos" = "firebrick3"
  )
)


# =============================================================================
# 3. NGFR-positive control cells: P0 vs P1
# =============================================================================

plot_pairwise_pca(
  "[Ctrl] N+ P0 vs P1.txt",
  "[PCA] [Ctrl] N+ P0 vs P1.pdf",
  group_function = function(x) {
    ifelse(grepl("_P0_", x), "P0", "P1")
  },
  group_colors = c(
    "P0" = "navy",
    "P1" = "firebrick3"
  )
)


# =============================================================================
# 4. NGFR-negative control cells: P0 vs P1
# =============================================================================

plot_pairwise_pca(
  "[Ctrl] N- P0 vs P1.txt",
  "[PCA] [Ctrl] N- P0 vs P1.pdf",
  group_function = function(x) {
    ifelse(grepl("_P0_", x), "P0", "P1")
  },
  group_colors = c(
    "P0" = "navy",
    "P1" = "firebrick3"
  )
)


# =============================================================================
# 5. P1: control vs IL1A, NGFR-positive and NGFR-negative combined
# =============================================================================

plot_pairwise_pca(
  "P1 Ctr vs IL1A.txt",
  "[PCA]P1 Ctr vs IL1A.pdf",
  group_function = function(x) {
    ifelse(
      grepl("Ctrl_NGFR_neg", x), "Ctrl_NGFR_neg",
      ifelse(
        grepl("Ctrl_NGFR_pos", x), "Ctrl_NGFR_pos",
        ifelse(
          grepl("IL1A_NGFR_neg", x), "IL1A_NGFR_neg",
          ifelse(grepl("IL1A_NGFR_pos", x), "IL1A_NGFR_pos", NA)
        )
      )
    )
  },
  group_colors = c(
    "Ctrl_NGFR_neg" = "#6666B3",
    "Ctrl_NGFR_pos" = "navy",
    "IL1A_NGFR_neg" = "#E27676",
    "IL1A_NGFR_pos" = "firebrick3"
  ),
  width = 8
)


# =============================================================================
# 6. P1 NGFR-positive: control vs IL1A
# =============================================================================

plot_pairwise_pca(
  "P1 N+ Ctr vs IL1A.txt",
  "[PCA]P1 N+ Ctr vs IL1A.pdf",
  group_function = function(x) {
    ifelse(grepl("Ctr|Ctrl", x, ignore.case = TRUE), "Ctrl", "IL1A")
  },
  group_colors = c(
    "Ctrl" = "navy",
    "IL1A" = "firebrick3"
  )
)


# =============================================================================
# 7. P1 NGFR-negative: control vs IL1A
# =============================================================================

plot_pairwise_pca(
  "P1 N- Ctr vs IL1A.txt",
  "[PCA]P1 N- Ctr vs IL1A.pdf",
  group_function = function(x) {
    ifelse(grepl("Ctr|Ctrl", x, ignore.case = TRUE), "Ctrl", "IL1A")
  },
  group_colors = c(
    "Ctrl" = "navy",
    "IL1A" = "firebrick3"
  )
)


# =============================================================================
# Helper: raw and normalized DESeq2 library sizes
# =============================================================================

plot_library_sizes <- function(dds, output_file) {

  raw_lib_sizes <- colSums(counts(dds))
  norm_lib_sizes <- colSums(counts(dds, normalized = TRUE))

  pdf(output_file, width = 20, height = 5)

  par(
    mfrow = c(1, 2),
    mar = c(7, 4, 4, 2)
  )

  bp1 <- barplot(
    raw_lib_sizes,
    col = "navy",
    names.arg = rep("", length(raw_lib_sizes)),
    main = "Raw Library Sizes",
    ylab = "Total Reads",
    xlab = ""
  )

  text(
    x = bp1,
    y = par("usr")[3] - 0.05 * max(raw_lib_sizes),
    labels = names(raw_lib_sizes),
    srt = 45,
    adj = 1,
    xpd = TRUE,
    cex = 0.8
  )

  bp2 <- barplot(
    norm_lib_sizes,
    col = "firebrick3",
    names.arg = rep("", length(norm_lib_sizes)),
    main = "Normalized Library Sizes",
    ylab = "Total Reads",
    xlab = ""
  )

  text(
    x = bp2,
    y = par("usr")[3] - 0.05 * max(norm_lib_sizes),
    labels = names(norm_lib_sizes),
    srt = 45,
    adj = 1,
    xpd = TRUE,
    cex = 0.8
  )

  dev.off()
}
