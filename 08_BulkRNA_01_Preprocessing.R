# =============================================================================
# 08_BulkRNA_01_Preprocessing.R
# Bulk RNA-seq count-matrix preparation and pairwise matrix export
# =============================================================================

library(readxl)

setwd("D:\\Eve_Bulk\\R")

# =============================================================================
# 1. Load the raw count matrix
# =============================================================================

raw <- read_excel("Eve_BulkRNA_counts_by_sample_sorted.xlsx", sheet = 1)

counts <- as.data.frame(raw)
rownames(counts) <- counts[[1]]
counts[[1]] <- NULL


# =============================================================================
# 2. Gene-level filtering
# =============================================================================
# Retain genes with:
#   1) total counts > 100 across all samples, and
#   2) non-zero counts in at least 50% of samples.
#
# These criteria follow the final filtering section of the original workflow.

n_samples <- ncol(counts)
gene_sums <- rowSums(counts)
expr_count <- rowSums(counts > 0)

keep <- (gene_sums > 100) & (expr_count >= 0.5 * n_samples)

counts_filtered <- counts[keep, ]

cat("Genes before filtering:", nrow(counts), "\n")
cat("Genes after filtering:", nrow(counts_filtered), "\n")

counts <- counts_filtered
rm(counts_filtered)


# =============================================================================
# 3. Reassign donor IDs
# =============================================================================
# Original numeric donor prefixes are converted to sequential IDs.
# The sample-name suffixes are otherwise retained unchanged.

old_names <- colnames(counts)

orig_ids <- as.numeric(sub("^([0-9]+).*", "\\1", old_names))
uniq_ids <- sort(unique(orig_ids))
new_ids <- match(orig_ids, uniq_ids)

rest <- sub("^[0-9]+", "", old_names)
new_names <- paste0(new_ids, rest)

colnames(counts) <- new_names

print(colnames(counts))


# =============================================================================
# 4. Remove samples with reassigned donor ID 3
# =============================================================================
# This exclusion is retained exactly from the final original workflow.

samples_to_remove <- grep("^3", colnames(counts), value = TRUE)

cat("Samples removed from the analysis:\n")
print(samples_to_remove)

counts <- counts[, !colnames(counts) %in% samples_to_remove, drop = FALSE]

cat("Remaining genes:", nrow(counts), "\n")
cat("Remaining samples:", ncol(counts), "\n")


# =============================================================================
# Helper function
# =============================================================================

write_pairwise_matrix <- function(mat, columns, output_file) {
  out <- mat[, columns, drop = FALSE]

  write.table(
    out,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA
  )

  invisible(out)
}


# =============================================================================
# 5. P0: NGFR-negative vs NGFR-positive
# =============================================================================

p0_cols <- grep("P0", colnames(counts), value = TRUE)
neg_cols <- grep("NGFR_neg", p0_cols, value = TRUE)
pos_cols <- grep("NGFR_pos", p0_cols, value = TRUE)

write_pairwise_matrix(
  counts,
  c(neg_cols, pos_cols),
  "P0 N- vs N+.txt"
)


# =============================================================================
# 6. P1 control organoids: NGFR-negative vs NGFR-positive
# =============================================================================

p1_ctrl_cols <- grep("P1_Ctrl", colnames(counts), value = TRUE)
neg_cols <- grep("NGFR_neg", p1_ctrl_cols, value = TRUE)
pos_cols <- grep("NGFR_pos", p1_ctrl_cols, value = TRUE)

write_pairwise_matrix(
  counts,
  c(neg_cols, pos_cols),
  "P1_Ctrl N- vs N+.txt"
)


# =============================================================================
# 7. NGFR-positive control cells: P0 vs P1
# =============================================================================

pos_cols <- grep("NGFR_pos", colnames(counts), value = TRUE)
pos_cols <- pos_cols[!grepl("IL1A", pos_cols)]

p0_cols <- grep("P0", pos_cols, value = TRUE)
p1_cols <- grep("P1", pos_cols, value = TRUE)

write_pairwise_matrix(
  counts,
  c(p0_cols, p1_cols),
  "[Ctrl] N+ P0 vs P1.txt"
)


# =============================================================================
# 8. NGFR-negative control cells: P0 vs P1
# =============================================================================

neg_cols <- grep("NGFR_neg", colnames(counts), value = TRUE)
neg_cols <- neg_cols[!grepl("IL1A", neg_cols)]

p0_cols <- grep("P0", neg_cols, value = TRUE)
p1_cols <- grep("P1", neg_cols, value = TRUE)

write_pairwise_matrix(
  counts,
  c(p0_cols, p1_cols),
  "[Ctrl] N- P0 vs P1.txt"
)


# =============================================================================
# 9. P1 organoids: control vs IL1A, NGFR-positive and NGFR-negative combined
# =============================================================================

p1_cols <- grep("P1", colnames(counts), value = TRUE)
ctrl_cols <- grep("Ctrl", p1_cols, value = TRUE)
il1a_cols <- grep("IL1A", p1_cols, value = TRUE)

write_pairwise_matrix(
  counts,
  c(ctrl_cols, il1a_cols),
  "P1 Ctr vs IL1A.txt"
)


# =============================================================================
# 10. P1 NGFR-positive cells: control vs IL1A
# =============================================================================

p1_pos_cols <- grep("P1", colnames(counts), value = TRUE)
p1_pos_cols <- p1_pos_cols[!grepl("NGFR_neg", p1_pos_cols)]

ctrl_cols <- grep("Ctr|Ctrl", p1_pos_cols, value = TRUE)
il1a_cols <- grep("IL1A", p1_pos_cols, value = TRUE)

write_pairwise_matrix(
  counts,
  c(ctrl_cols, il1a_cols),
  "P1 N+ Ctr vs IL1A.txt"
)


# =============================================================================
# 11. P1 NGFR-negative cells: control vs IL1A
# =============================================================================

p1_neg_cols <- grep("P1", colnames(counts), value = TRUE)
p1_neg_cols <- p1_neg_cols[!grepl("NGFR_pos", p1_neg_cols)]

ctrl_cols <- grep("Ctr|Ctrl", p1_neg_cols, value = TRUE)
il1a_cols <- grep("IL1A", p1_neg_cols, value = TRUE)

write_pairwise_matrix(
  counts,
  c(ctrl_cols, il1a_cols),
  "P1 N- Ctr vs IL1A.txt"
)
