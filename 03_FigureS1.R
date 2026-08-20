# =============================================================================
# Figure S1
# Single-cell overview and epithelial marker characterization
# =============================================================================

.libPaths(c("~/SeuratV4", .libPaths()))

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)
library(Nebulosa)
library(ggrastr)
library(scales)

setwd("E:\\EoE")

# Load project-specific plotting functions, including Nebulosa helpers and
# density_umap_triple() used for Figure S1I.
source("EoE_source.R")


# -----------------------------------------------------------------------------
# Load the complete processed object generated in 01_EoE_scRNA_processing.R
# -----------------------------------------------------------------------------

test.seu <- readRDS(
  "./[2025](06.1)EoE_raw_100000_NoRemission_REcelltype_25-0.7-30-0.3.RDS"
)

DefaultAssay(test.seu) <- "RNA"

color_scheme <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"
)

color_cc <- c("#E15759", "#F0C420", "#1965B0")
color_ccm <- c("gray95", "#B72230")


# =============================================================================
# Figure S1A
# UMAP of all major cell types
# =============================================================================

figS1A <- DimPlot(
  test.seu,
  reduction = "umap",
  group.by = "celltype",
  cols = color_scheme,
  pt.size = 0.4,
  label = TRUE,
  repel = TRUE
) +
  NoLegend() +
  NoAxes()

print(figS1A)

ggsave(
  filename = "FigureS1A_All_celltypes_UMAP.pdf",
  plot = figS1A,
  device = "pdf",
  width = 9,
  height = 8,
  units = "cm"
)


# =============================================================================
# Figure S1B
# Nebulosa density plots for lineage markers
# =============================================================================

figS1B_genes <- c(
  "COL1A1",
  "CPE",
  "HBB",
  "PLVAP",
  "KRT5",
  "PTPRC"
)

figS1B_genes <- intersect(figS1B_genes, rownames(test.seu))

# Use Nebulosa kernel-density estimation, following the custom density plotting
# approach implemented in EoE_source.R.
figS1B_list <- lapply(figS1B_genes, function(gene) {
  p <- plot_umap_density_ccm(
    seurat_obj = test.seu,
    gene = gene,
    reduction = "umap",
    method = "ks",
    adjust = 1,
    size = 0.5,
    color_ccm = color_ccm,
    raster = TRUE,
    legend_title = "Density"
  ) +
    ggtitle(gene) +
    theme(
      plot.title = element_text(face = "italic", hjust = 0),
      legend.position = "right"
    )

  p
})

figS1B <- wrap_plots(figS1B_list, ncol = 3)

print(figS1B)

ggsave(
  filename = "FigureS1B_Lineage_markers_Nebulosa.pdf",
  plot = figS1B,
  device = "pdf",
  width = 16,
  height = 10,
  units = "cm"
)


# =============================================================================
# Figure S1C
# Dot plot of immune-cell markers
# =============================================================================

figS1C_markers <- c(
  # Macrophage / monocyte markers
  "CD14", "C1QA", "VCAN", "FCN1",

  # Mast-cell markers
  "KIT", "TPSAB1",

  # Dendritic-cell markers
  "CD207", "CD1E",

  # Plasma/B-cell markers
  "IGKC", "IGHA1", "BANK1", "CD19",

  # NK-cell markers
  "XCL1", "NCAM1",

  # T-cell markers
  "CD3G", "CD3D"
)

figS1C_markers <- intersect(figS1C_markers, rownames(test.seu))

figS1C_celltypes <- c(
  "Macrophage",
  "Mast cell",
  "Dendritic cell",
  "Monocyte",
  "Plasma cell",
  "B cell",
  "NK cell",
  "T cell"
)

figS1C_obj <- subset(
  test.seu,
  subset = celltype2 %in% figS1C_celltypes
)

figS1C_obj$celltype2 <- factor(
  figS1C_obj$celltype2,
  levels = figS1C_celltypes
)

figS1C <- DotPlot(
  figS1C_obj,
  features = figS1C_markers,
  group.by = "celltype2",
  dot.scale = 6
) +
  RotatedAxis() +
  scale_color_gradient(
    low = "gray90",
    high = "#B72230",
    name = "Ave.\nExpr.\nScaled"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Immune cell markers"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(figS1C)

ggsave(
  filename = "FigureS1C_Immune_marker_DotPlot.pdf",
  plot = figS1C,
  device = "pdf",
  width = 12,
  height = 9,
  units = "cm"
)


# =============================================================================
# Figure S1D
# Epithelial marker expression across the UMAP
# =============================================================================

figS1D_genes <- c(
  "PDPN",
  "COL17A1",
  "TP63",
  "KRT6B",
  "TGM1",
  "FLG",
  "NGFR",
  "TP73",
  "KRT14",
  "KLF4",
  "KRT4",
  "KRT78"
)

figS1D_genes <- intersect(figS1D_genes, rownames(test.seu))

# Use the complete object, as in the supplementary figure.
figS1D <- FeaturePlot(
  test.seu,
  features = figS1D_genes,
  reduction = "umap",
  ncol = 6,
  order = TRUE,
  pt.size = 0.15,
  cols = c("gray95", "#B72230")
) &
  NoAxes() &
  theme(
    plot.title = element_text(face = "italic", size = 10)
  )

print(figS1D)

ggsave(
  filename = "FigureS1D_Epithelial_marker_FeaturePlot.pdf",
  plot = figS1D,
  device = "pdf",
  width = 20,
  height = 8,
  units = "cm"
)


# =============================================================================
# Figure S1E
# Cell-cycle phase on the UMAP
# =============================================================================

# Cell-cycle scores were calculated during processing using CellCycleScoring().
test.seu$Phase <- factor(
  test.seu$Phase,
  levels = c("G2M", "S", "G1")
)

figS1E <- DimPlot(
  test.seu,
  group.by = "Phase",
  reduction = "umap",
  pt.size = 0.2,
  cols = color_cc,
  label = FALSE
) +
  NoAxes()

print(figS1E)

ggsave(
  filename = "FigureS1E_CellCycle_UMAP.pdf",
  plot = figS1E,
  device = "pdf",
  width = 8,
  height = 8,
  units = "cm"
)


# =============================================================================
# Figure S1F
# Cell-cycle phase distribution across epithelial populations
# =============================================================================

figS1F_df <- test.seu@meta.data %>%
  filter(celltype2 %in% c(
    "Basal",
    "Cycling Suprabasal",
    "Suprabasal",
    "Apical"
  )) %>%
  mutate(
    epithelial_group = recode(
      celltype2,
      "Basal" = "Basal",
      "Cycling Suprabasal" = "Prol. Cells",
      "Suprabasal" = "SB",
      "Apical" = "Apical"
    ),
    epithelial_group = factor(
      epithelial_group,
      levels = c("Basal", "Prol. Cells", "SB", "Apical")
    ),
    Phase = factor(Phase, levels = c("G2M", "S", "G1"))
  ) %>%
  count(epithelial_group, Phase, name = "n") %>%
  group_by(epithelial_group) %>%
  mutate(Percentage = n / sum(n) * 100) %>%
  ungroup()

figS1F <- ggplot(
  figS1F_df,
  aes(x = epithelial_group, y = Percentage, fill = Phase)
) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = color_cc) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Percentage [%]",
    fill = NULL,
    title = "Cell Cycle Phase\nDistribution"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "none"
  )

print(figS1F)

ggsave(
  filename = "FigureS1F_CellCycle_distribution.pdf",
  plot = figS1F,
  device = "pdf",
  width = 7,
  height = 8,
  units = "cm"
)




# =============================================================================
# Figure S1I
# Density overlap of epithelial marker pairs
# =============================================================================

# The manuscript panel contains two marker-pair comparisons:
#   top row:    NGFR versus KRT14
#   bottom row: KRT14 versus KRT4
#
# density_umap_triple() is defined in EoE_source.R and uses Nebulosa density
# when available, with KNN smoothing as a fallback.

# Use the epithelial-only object for the overlap plots.
test.seu2 <- readRDS(
  "E:\\EoE\\[2025](06.2-Epi)EoE_raw_100000_NoRemission_REcelltype_25-0.7-30-0.3.RDS"
)

# Figure S1I, top row: NGFR / KRT14 / overlap.
density_umap_triple(
  test.seu2,
  "NGFR",
  "KRT14",
  radialize = TRUE,
  radial_scale = 0.5,
  point_size = 0.8,
  col_g1 = "#155696",
  col_g2 = "#CCA71B",
  col_min = "#22763F",
  pdf_file = "FigureS1I_NGFR_KRT14_overlap.pdf",
  width = 16,
  height = 5,
  knn_k = 100
)

# Figure S1I, bottom row: KRT14 / KRT4 / overlap.
density_umap_triple(
  test.seu2,
  "KRT14",
  "KRT4",
  radialize = TRUE,
  radial_scale = 0.5,
  point_size = 0.8,
  col_g1 = "#155696",
  col_g2 = "#CCA71B",
  col_min = "#22763F",
  pdf_file = "FigureS1I_KRT14_KRT4_overlap.pdf",
  width = 16,
  height = 5,
  knn_k = 100
)
