# =============================================================================
# Figure 3I
# Projection of organoid-derived NGFR-positive and NGFR-negative signatures
# onto the human esophageal epithelial single-cell RNA-seq UMAP
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
# Figure 3I - NGFR-positive organoid signature
# =============================================================================

# Top 50 genes enriched in NGFR-positive cells from the P1 control organoid
# NGFR-positive vs NGFR-negative bulk RNA-seq comparison.
NGFR_pos_top50 <- c(
  "NGFR",     "CCK",     "LSP1",    "FST",     "CMTM7",
  "DLK2",     "RAMP1",   "SPARC",   "RGS3",    "MNS1",
  "LIG1",     "KRT15",   "CTNNAL1", "NRG1",    "SERPINH1",
  "HJURP",    "RAD18",   "DKK1",    "ACSL4",   "ASCL2",
  "TOP2A",    "PSMB9",   "G0S2",    "UBE2C",   "EMP3",
  "ASS1",     "SLC7A8",  "TSPAN4",  "NUF2",    "CXCL2",
  "TNFRSF21", "HMGB2",   "HIRIP3",  "FIGNL1",  "LAMA3",
  "MT2A",     "ALDH3A1", "E2F7",    "DEPDC1",  "CDK1",
  "SESN3",    "KIF23",   "TCEAL3",  "CDCA3",   "CLSPN",
  "AURKB",    "ASPM",    "INHBA",   "ANLN",    "KIF11"
)

# Keep only genes present in the single-cell object.
NGFR_pos_top50_f <- intersect(NGFR_pos_top50, rownames(test.seu))

missing_pos <- setdiff(NGFR_pos_top50, NGFR_pos_top50_f)
if (length(missing_pos) > 0) {
  message(
    "NGFR-positive signature genes not found in the Seurat object and ignored: ",
    paste(missing_pos, collapse = ", ")
  )
}

# Calculate the NGFR-positive module score.
test.seu <- AddModuleScore(
  test.seu,
  features = list(NGFR_pos_top50_f),
  name = "NGFR_pos_signature"
)


# =============================================================================
# Figure 3I - NGFR-negative organoid signature
# =============================================================================

# Top 50 genes enriched in NGFR-negative cells from the P1 control organoid
# NGFR-positive vs NGFR-negative bulk RNA-seq comparison.
NGFR_neg_top50 <- c(
  "IGFL2",     "BSPRY",    "IVL",      "DSG1",      "SULT2B1",
  "MYO16",     "KRT6C",    "FAM25C",   "GDPD3",     "RRAD",
  "CALB2",     "SLPI",     "KLK7",     "LINC00520", "PI3",
  "IL36RN",    "RDH12",    "NEBL",     "FABP5P3",   "KRTDAP",
  "CLIC3",     "DEGS2",    "SPRR1B",   "BBOX1",     "KRT16",
  "ASPG",      "SBSN",     "ALDH3B2",  "IGFL1",     "ZNF750",
  "CDA",       "FAM25A",   "FABP5P2",  "KRT6B",     "HES4",
  "ADAP2",     "FABP5P9",  "NECTIN4",  "LY6D",      "CYP2C18",
  "SERPINB13", "CDH16",    "LY6G6C",   "FABP5P7",   "RHOV",
  "FABP5",     "GGT6",     "C10orf99", "SERPINB3",  "SPRR1A"
)

# Keep only genes present in the single-cell object.
NGFR_neg_top50_f <- intersect(NGFR_neg_top50, rownames(test.seu))

missing_neg <- setdiff(NGFR_neg_top50, NGFR_neg_top50_f)
if (length(missing_neg) > 0) {
  message(
    "NGFR-negative signature genes not found in the Seurat object and ignored: ",
    paste(missing_neg, collapse = ", ")
  )
}

# Calculate the NGFR-negative module score.
test.seu <- AddModuleScore(
  test.seu,
  features = list(NGFR_neg_top50_f),
  name = "NGFR_neg_signature"
)


# =============================================================================
# Normalize module scores to 0-1
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
# Figure 3I - Plot
# =============================================================================

# Epithelial population colors used for the reference UMAP.
epithelial_colors <- c(
  "Basal" = "#EFE2AA",
  "Cycling Suprabasal" = "#FDB462",
  "Suprabasal" = "#EF9A9A",
  "Apical" = "#E15759"
)

# Reference epithelial UMAP.
p_reference <- DimPlot(
  test.seu,
  reduction = "umap",
  group.by = "celltype2",
  cols = epithelial_colors,
  pt.size = 0.5,
  label = FALSE
) +
  ggtitle("Esophageal Epithelium") +
  NoLegend() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# Signature colors used in the original analysis.
signature_colors <- c("#F0F0F0", "#C0203C")

# NGFR-negative organoid signature.
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

# NGFR-positive organoid signature.
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

# Assemble Figure 3I.
Fig3I <- p_reference + p_NGFR_neg + p_NGFR_pos +
  plot_layout(ncol = 3)

print(Fig3I)

ggsave(
  filename = "Figure3I_NGFR_signature_scores.pdf",
  plot = Fig3I,
  device = "pdf",
  width = 9,
  height = 3,
  units = "in"
)
