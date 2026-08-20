# =============================================================================
# Figure 1
# Esophageal epithelial single-cell analysis
# =============================================================================

.libPaths(c("~/SeuratV4", .libPaths()))

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)
library(ComplexHeatmap)
library(circlize)

setwd("E:\\EoE")

# Custom plotting functions used in the project.
# Keep EoE_source.R in the working directory if these helper functions are needed.
source("EoE_source.R")


# -----------------------------------------------------------------------------
# Load the epithelial-only object generated in 01_EoE_scRNA_processing.R
# -----------------------------------------------------------------------------

test.seu2 <- test.seu2 <- subset(test.seu, subset = celltype3 == "Epithelial")

DefaultAssay(test.seu2) <- "RNA"

color_scheme <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"
)


# =============================================================================
# Figure 1A
# UMAP of epithelial populations
# =============================================================================

# Keep the epithelial populations in the order shown in the manuscript.
test.seu2$celltype2 <- factor(
  test.seu2$celltype2,
  levels = c("Basal", "Cycling Suprabasal", "Suprabasal", "Apical")
)

fig1A <- DimPlot(
  test.seu2,
  reduction = "umap",
  group.by = "celltype2",
  cols = color_scheme,
  pt.size = 0.5,
  label = FALSE
) +
  NoLegend() +
  NoAxes() +
  ggtitle("Esophageal Epithelium") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(5, 5, 5, 5)
  )

print(fig1A)

ggsave(
  filename = "Figure1A_Epithelial_UMAP.pdf",
  plot = fig1A,
  device = "pdf",
  width = 8,
  height = 7,
  units = "cm"
)


# =============================================================================
# Figure 1B
# Heatmap of epithelial differentiation markers
# =============================================================================

# Marker order follows the basal-to-apical differentiation axis in Figure 1B.
fig1B_genes <- c(
  "COL17A1",
  "IGFBP3",
  "TP73",
  "PDPN",
  "NGFR",
  "KRT15",
  "KRT14",
  "TP63",
  "KLF4",
  "MKI67",
  "CDK1",
  "KRT6B",
  "KRT6C",
  "TGM1",
  "KRT4",
  "CNFN",
  "KRT78"
)

fig1B_genes <- intersect(fig1B_genes, rownames(test.seu2))

# Scale the marker genes if they are not already present in scale.data.
test.seu2 <- ScaleData(
  test.seu2,
  features = fig1B_genes,
  verbose = FALSE
)

# Downsample each epithelial population for visualization only.
# This does not modify the saved Seurat object.
set.seed(123)
cells_fig1B <- unlist(
  lapply(levels(test.seu2$celltype2), function(ct) {
    cells <- colnames(test.seu2)[test.seu2$celltype2 == ct]
    sample(cells, min(length(cells), 500))
  })
)

fig1B_obj <- subset(test.seu2, cells = cells_fig1B)
Idents(fig1B_obj) <- "celltype2"

fig1B <- DoHeatmap(
  fig1B_obj,
  features = fig1B_genes,
  group.by = "celltype2",
  slot = "scale.data",
  raster = TRUE,
  disp.min = -2,
  disp.max = 2
) +
  scale_fill_gradient2(
    low = "#5B8DB8",
    mid = "white",
    high = "#C75B5B",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish,
    name = "Scaled\nExpr."
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank()
  )

print(fig1B)

ggsave(
  filename = "Figure1B_Epithelial_marker_heatmap.pdf",
  plot = fig1B,
  device = "pdf",
  width = 12,
  height = 8,
  units = "cm"
)
