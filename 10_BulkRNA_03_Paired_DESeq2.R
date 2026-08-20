# =============================================================================
# 10_BulkRNA_03_Paired_DESeq2.R
# Paired differential-expression analysis
# =============================================================================

library(DESeq2)
library(openxlsx)

setwd("D:\\Eve_Bulk\\R\\Pairwise")


# =============================================================================
# Helper functions
# =============================================================================

read_count_matrix <- function(input_file) {
  mat <- read.table(
    input_file,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  mat <- as.matrix(mat)

  na_genes <- rowSums(is.na(mat)) > 0

  if (any(na_genes)) {
    message(sum(na_genes), " genes contain NA values and will be removed.")
    mat <- mat[!na_genes, , drop = FALSE]
  }

  round(mat)
}


run_paired_deseq2 <- function(
    counts,
    patient,
    condition,
    condition_levels,
    contrast,
    result_xlsx,
    result_txt,
    normalized_xlsx,
    rlog_xlsx
) {

  coldata <- data.frame(
    row.names = colnames(counts),
    patient = factor(patient),
    condition = factor(condition, levels = condition_levels)
  )

  stopifnot(identical(rownames(coldata), colnames(counts)))

  # Paired design used in the final analysis.
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = ~ patient + condition
  )

  dds <- DESeq(dds)

  res <- results(
    dds,
    contrast = contrast
  )

  lfc <- res$log2FoldChange
  padj <- res$padj
  pvalue <- res$pvalue

  cat(
    "Upregulated genes (|log2FC| > 0.585, padj < 0.05): ",
    sum(lfc > 0.585 & padj < 0.05, na.rm = TRUE),
    "\n"
  )

  cat(
    "Downregulated genes (|log2FC| > 0.585, padj < 0.05): ",
    sum(lfc < -0.585 & padj < 0.05, na.rm = TRUE),
    "\n"
  )

  cat(
    "Upregulated genes (log2FC > 0.585, pvalue < 0.05): ",
    sum(lfc > 0.585 & pvalue < 0.05, na.rm = TRUE),
    "\n"
  )

  cat(
    "Downregulated genes (log2FC < -0.585, pvalue < 0.05): ",
    sum(lfc < -0.585 & pvalue < 0.05, na.rm = TRUE),
    "\n"
  )

  res_df <- as.data.frame(res)
  res_df$geneNames <- rownames(res_df)
  res_df <- res_df[, c("geneNames", setdiff(names(res_df), "geneNames"))]
  res_df <- res_df[order(res_df$padj), ]

  write.xlsx(
    res_df,
    file = result_xlsx,
    rowNames = FALSE
  )

  write.table(
    res_df,
    file = result_txt,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  norm_counts <- counts(dds, normalized = TRUE)

  norm_df <- as.data.frame(norm_counts)
  norm_df$gene <- rownames(norm_df)
  norm_df <- norm_df[, c("gene", setdiff(colnames(norm_df), "gene"))]

  write.xlsx(
    norm_df,
    file = normalized_xlsx,
    rowNames = FALSE
  )

  rld <- rlog(dds, blind = FALSE)
  rlog_mat <- assay(rld)

  rlog_df <- as.data.frame(rlog_mat)
  rlog_df$gene <- rownames(rlog_df)
  rlog_df <- rlog_df[, c("gene", setdiff(colnames(rlog_df), "gene"))]

  write.xlsx(
    rlog_df,
    file = rlog_xlsx,
    rowNames = FALSE
  )

  invisible(list(
    dds = dds,
    result = res,
    result_df = res_df,
    normalized_counts = norm_counts,
    rlog = rld
  ))
}


# =============================================================================
# 1. P0 NGFR-positive vs NGFR-negative
# =============================================================================

counts <- read_count_matrix("P0 N- vs N+.txt")
samples <- colnames(counts)
parts <- strsplit(samples, "_")

patient <- sapply(parts, function(x) paste(x[1], x[2], sep = "_"))
condition <- sapply(parts, function(x) paste(x[3], x[4], sep = "_"))

P0_NGFR <- run_paired_deseq2(
  counts = counts,
  patient = patient,
  condition = condition,
  condition_levels = c("NGFR_neg", "NGFR_pos"),
  contrast = c("condition", "NGFR_pos", "NGFR_neg"),
  result_xlsx = "{Pairwise}[DESeq2]P0_N+_vs_N-.xlsx",
  result_txt = "{Pairwise}[DESeq2]P0_N+_vs_N-.txt",
  normalized_xlsx = "[normalized counts]P0_N+_vs_N-.xlsx",
  rlog_xlsx = "[log counts]P0_N+_vs_N-.xlsx"
)


# =============================================================================
# 2. P1 control NGFR-positive vs NGFR-negative
# =============================================================================

counts <- read_count_matrix("P1_Ctrl N- vs N+.txt")
samples <- colnames(counts)
parts <- strsplit(samples, "_")

patient <- sapply(parts, function(x) paste(x[1], x[2], sep = "_"))
condition <- sapply(parts, function(x) paste(x[4], x[5], sep = "_"))

P1_NGFR <- run_paired_deseq2(
  counts = counts,
  patient = patient,
  condition = condition,
  condition_levels = c("NGFR_neg", "NGFR_pos"),
  contrast = c("condition", "NGFR_pos", "NGFR_neg"),
  result_xlsx = "{Pairwise}[DESeq2]P1_Ctrl N- vs N+.xlsx",
  result_txt = "{Pairwise}[DESeq2]P1_Ctrl N- vs N+.txt",
  normalized_xlsx = "[normalized counts]P1_Ctrl N- vs N+.xlsx",
  rlog_xlsx = "[log counts]P1_Ctrl N- vs N+.xlsx"
)


# =============================================================================
# 3. NGFR-positive control cells: P1 vs P0
# =============================================================================

counts <- read_count_matrix("[Ctrl] N+ P0 vs P1.txt")
samples <- colnames(counts)
parts <- strsplit(samples, "_")

patient <- sapply(parts, `[`, 1)
condition <- sapply(parts, `[`, 2)

Npos_P1_vs_P0 <- run_paired_deseq2(
  counts = counts,
  patient = patient,
  condition = condition,
  condition_levels = c("P0", "P1"),
  contrast = c("condition", "P1", "P0"),
  result_xlsx = "{Pairwise}[DESeq2]N+ P1 vs P0.xlsx",
  result_txt = "{Pairwise}[DESeq2]N+ P1 vs P0.txt",
  normalized_xlsx = "[normalized counts]N+ P1 vs P0.xlsx",
  rlog_xlsx = "[log counts]N+ P1 vs P0.xlsx"
)


# =============================================================================
# 4. NGFR-negative control cells: P1 vs P0
# =============================================================================

counts <- read_count_matrix("[Ctrl] N- P0 vs P1.txt")
samples <- colnames(counts)
parts <- strsplit(samples, "_")

patient <- sapply(parts, `[`, 1)
condition <- sapply(parts, `[`, 2)

Nneg_P1_vs_P0 <- run_paired_deseq2(
  counts = counts,
  patient = patient,
  condition = condition,
  condition_levels = c("P0", "P1"),
  contrast = c("condition", "P1", "P0"),
  result_xlsx = "{Pairwise}[DESeq2]N- P1 vs P0.xlsx",
  result_txt = "{Pairwise}[DESeq2]N- P1 vs P0.txt",
  normalized_xlsx = "[normalized counts]N- P1 vs P0.xlsx",
  rlog_xlsx = "[log counts]N- P1 vs P0.xlsx"
)


# =============================================================================
# 5. P1 control vs IL1A, NGFR-positive and NGFR-negative combined
# =============================================================================

counts <- read_count_matrix("P1 Ctr vs IL1A.txt")
samples <- colnames(counts)

patient <- sub("_.*", "", samples)
condition <- ifelse(grepl("Ctrl|Ctr", samples), "Ctrl", "IL1A")

P1_IL1A_all <- run_paired_deseq2(
  counts = counts,
  patient = patient,
  condition = condition,
  condition_levels = c("Ctrl", "IL1A"),
  contrast = c("condition", "IL1A", "Ctrl"),
  result_xlsx = "{Pairwise}[DESeq2]P1 IL1A vs Ctr.xlsx",
  result_txt = "{Pairwise}[DESeq2]P1 IL1A vs Ctr.txt",
  normalized_xlsx = "[DESeq2]P1 IL1A vs Ctr normalized counts.xlsx",
  rlog_xlsx = "[DESeq2]P1 IL1A vs Ctr rlog counts.xlsx"
)


# =============================================================================
# 6. P1 NGFR-positive: control vs IL1A
# =============================================================================

counts <- read_count_matrix("P1 N+ Ctr vs IL1A.txt")
samples <- colnames(counts)

patient <- sub("_.*", "", samples)
condition <- ifelse(grepl("Ctrl|Ctr", samples), "Ctrl", "IL1A")

P1_Npos_IL1A <- run_paired_deseq2(
  counts = counts,
  patient = patient,
  condition = condition,
  condition_levels = c("Ctrl", "IL1A"),
  contrast = c("condition", "IL1A", "Ctrl"),
  result_xlsx = "{Pairwise}[DESeq2]P1 N+ IL1A vs Ctr.xlsx",
  result_txt = "{Pairwise}[DESeq2]P1 N+ IL1A vs Ctr.txt",
  normalized_xlsx = "[DESeq2]P1 N+ IL1A vs Ctr normalized counts.xlsx",
  rlog_xlsx = "[DESeq2]P1 N+ IL1A vs Ctr rlog counts.xlsx"
)


# =============================================================================
# 7. P1 NGFR-negative: control vs IL1A
# =============================================================================

counts <- read_count_matrix("P1 N- Ctr vs IL1A.txt")
samples <- colnames(counts)

patient <- sub("_.*", "", samples)
condition <- ifelse(grepl("Ctrl|Ctr", samples), "Ctrl", "IL1A")

P1_Nneg_IL1A <- run_paired_deseq2(
  counts = counts,
  patient = patient,
  condition = condition,
  condition_levels = c("Ctrl", "IL1A"),
  contrast = c("condition", "IL1A", "Ctrl"),
  result_xlsx = "{Pairwise}[DESeq2]P1 N- IL1A vs Ctr.xlsx",
  result_txt = "{Pairwise}[DESeq2]P1 N- IL1A vs Ctr.txt",
  normalized_xlsx = "[DESeq2]P1 N- I1A vs Ctr normalized counts.xlsx",
  rlog_xlsx = "[DESeq2]P1 N- IL1A vs Ctr rlog counts.xlsx"
)
