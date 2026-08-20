# =============================================================================
# Figure 2I
# Projection of NGFR-positive and NGFR-negative bulk RNA-seq signatures
# onto the epithelial single-cell RNA-seq UMAP
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)

# -----------------------------------------------------------------------------
# Load the epithelial-only single-cell object
# -----------------------------------------------------------------------------

test.seu <- readRDS("D:/Eve_Bulk/R/[paper]Health___EpiONLY.RDS")


# =============================================================================
# Figure 2I - NGFR-positive signature
# =============================================================================

# Top 50 genes enriched in NGFR-positive cells from the P0 tissue
# NGFR-positive vs NGFR-negative bulk RNA-seq comparison.
NGFR_pos_top50 <- c(
  "ACKR3",     "NTF3",    "NPPC",     "IQCA1",    "CH25H",
  "DLL1",      "DLK2",    "C1QTNF12", "CYBRD1",   "GEM",
  "CCL2",      "TSLP",    "NGFR",     "MXRA5",    "CYYR1",
  "TGFBI",     "LFNG",    "MOXD1",    "NRG2",     "SYNPO2",
  "C1R",       "RAMP1",   "JAG2",     "CXCL14",   "SERPINF1",
  "RNF165",    "WNT10A",  "PDPN",     "CTNNAL1",  "NNAT",
  "KHDC1-AS1", "IL1R2",   "VIT",      "TNS1",     "SNCA",
  "RNU2-1",    "CLIP2",   "DOCK10",   "KIRREL1",  "LAMB4",
  "CHL1",      "COL14A1", "COL17A1",  "EFEMP1",   "FOXP2",
  "DCN",       "ABCA1",   "BCAM",     "TNS3",     "PHYHIP"
)

# Keep only genes present in the Seurat object.
NGFR_pos_top50_f <- intersect(NGFR_pos_top50, rownames(test.seu))

missing_pos <- setdiff(NGFR_pos_top50, NGFR_pos_top50_f)
if (length(missing_pos) > 0) {
  message(
    "NGFR-positive signature genes not found in the Seurat object and ignored: ",
    paste(missing_pos, collapse = ", ")
  )
}

# Calculate the module score.
test.seu <- AddModuleScore(
  test.seu,
  features = list(NGFR_pos_top50_f),
  name = "NGFR_pos_signature"
)


# =============================================================================
# Figure 2I - NGFR-negative signature
# =============================================================================

# Top 50 genes enriched in NGFR-negative cells from the P0 tissue
# NGFR-positive vs NGFR-negative bulk RNA-seq comparison.
NGFR_neg_top50 <- c(
  "IGHG4",    "TCN1",    "SERPINA1", "PLA2G2A", "ALDH1A1",
  "IGLC3",    "PPP1R1B", "GSTA1",    "KRT1",    "RNASE1",
  "BSPRY",    "NOSTRIN", "CLDN18",   "KALRN",   "CALML5",
  "C10orf99", "SPAG17",  "KRT13",    "IGKC",    "GABRP",
  "DEGS2",    "FAM3B",   "CXCL17",   "CYP2C18", "VSIG1",
  "BBOX1",    "CES1P1",  "MTARC1",   "SCIN",    "FCGBP",
  "MGLL",     "SERPINB3","SLC39A2",  "SPARCL1", "LINC01610",
  "SERPINB4", "PTGS1",   "ALDH3B2",  "IVL",     "MT4",
  "AKR7A3",   "RHCG",    "KRT4",     "ANKRD22", "TSPAN15",
  "IGFL2-AS1","CLDN4",   "GBP6",     "CRABP2",  "DSC2"
)

# Keep only genes present in the Seurat object.
NGFR_neg_top50_f <- intersect(NGFR_neg_top50, rownames(test.seu))

missing_neg <- setdiff(NGFR_neg_top50, NGFR_neg_top50_f)
if (length(missing_neg) > 0) {
  message(
    "NGFR-negative signature genes not found in the Seurat object and ignored: ",
    paste(missing_neg, collapse = ", ")
  )
}

# Calculate the module score.
test.seu <- AddModuleScore(
  test.seu,
  features = list(NGFR_neg_top50_f),
  name = "NGFR_neg_signature"
)


# =============================================================================
# Normalize the two module scores to 0-1
# =============================================================================

norm01 <- function(x) {
  rng <- range(x, na.rm = TRUE)

  if (diff(rng) == 0) {
    return(rep(0, length(x)))
  }

  (x - rng[1]) / diff(rng)
}

test.seu$NGFR_pos_signature_norm <- norm01(
  test.seu$NGFR_pos_signature1
)

test.seu$NGFR_neg_signature_norm <- norm01(
  test.seu$NGFR_neg_signature1
)


# =============================================================================
# Figure 2I - Plot signatures on the epithelial UMAP
# =============================================================================

signature_colors <- c("#F0F0F0", "#C0203C")

p_NGFR_neg <- FeaturePlot(
  test.seu,
  features = "NGFR_neg_signature_norm",
  reduction = "umap",
  cols = signature_colors,
  order = TRUE
) +
  ggtitle("NGFR- Signature") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

p_NGFR_pos <- FeaturePlot(
  test.seu,
  features = "NGFR_pos_signature_norm",
  reduction = "umap",
  cols = signature_colors,
  order = TRUE
) +
  ggtitle("NGFR+ Signature") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# Combine the two signature panels.
Fig2I <- p_NGFR_neg + p_NGFR_pos

print(Fig2I)

ggsave(
  filename = "Figure2I_NGFR_signature_scores.pdf",
  plot = Fig2I,
  device = "pdf",
  width = 6,
  height = 3,
  units = "in"
)
