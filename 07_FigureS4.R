# =============================================================================
# Figure S4
# Intraepithelial IL1 signaling in EoE
# Related to Figure 4
# =============================================================================

.libPaths(c("~/SeuratV4", .libPaths()))

library(Seurat)
library(qs)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(scales)

setwd("E:\\2025Plot")
source("E:/2025Plot/EoE_source.R")

setwd("E:\\2025Plot\\IL1A paper")


# -----------------------------------------------------------------------------
# Single-cell objects used for Figure S4
# -----------------------------------------------------------------------------

test.seu <- readRDS("E:/2025Plot/IL1A paper/[paper]EoE_all.RDS")

# Epithelial-only object used for Figure S4A-E and S4G.
test.seu4 <- subset(test.seu, subset = celltype3 == "Epithelial")
test.seu4 <- subset(
  test.seu4,
  subset = UMAP_1 >= -7 & UMAP_2 <= 13.5 & UMAP_2 >= 0.5
)

test.seu4$celltype2 <- factor(
  test.seu4$celltype2,
  levels = c("Basal", "Cycling Suprabasal", "Suprabasal", "Apical")
)

color_scheme <- c(
  "#EFE2AA", "#EBAEA9", "#83B4EF", "#BFA6C9",
  "#7EC0C3", "#F5E0BA", "#8ECFF8", "#AAD0AC",
  "#DBC9B3", "#89B780", "#EED0E0", "#F5D8D0",
  "#95A6DA", "#CB95BB", "#AED0DF"
)

color_cc <- c(
  "G1"  = "#E1575970",
  "S"   = "#F0C42070",
  "G2M" = "#1965B070"
)

cond_cols <- c(
  "Health" = "gray75",
  "EoE" = "#b40001",
  "Remission" = "#3D6B98"
)

# =============================================================================
# Figure S4A
# Epithelial UMAP shown separately for healthy and EoE samples
# =============================================================================

figS4A <- DimPlot(
  test.seu4,
  reduction = "umap",
  group.by = "celltype2",
  split.by = "condition",
  cols = color_scheme,
  pt.size = 0.4,
  label = FALSE,
  ncol = 2
) +
  NoLegend()

figS4A

ggsave(
  filename = "FigureS4A_Epithelial_UMAP_Health_vs_EoE.pdf",
  plot = figS4A,
  device = "pdf",
  width = 20,
  height = 10,
  units = "cm"
)


# =============================================================================
# Figure S4B
# Relative epithelial subpopulation distribution in healthy and EoE samples
# =============================================================================

df_pct <- test.seu4@meta.data %>%
  filter(condition %in% c("Health", "EoE")) %>%
  group_by(condition, celltype2) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(perc = count / sum(count) * 100) %>%
  ungroup()

df_pct$celltype2 <- factor(
  df_pct$celltype2,
  levels = c("Basal", "Cycling Suprabasal", "Suprabasal", "Apical")
)

figS4B <- ggplot(
  df_pct,
  aes(x = condition, y = perc, fill = celltype2)
) +
  geom_bar(
    stat = "identity",
    position = position_stack(reverse = TRUE),
    width = 0.7
  ) +
  geom_text(
    aes(label = paste0(round(perc, 1), "%")),
    position = position_stack(reverse = TRUE, vjust = 0.5),
    color = "black",
    size = 3
  ) +
  scale_fill_manual(
    values = c("#EFE2AA", "#EBAEA9", "#83B4EF", "#BFA6C9")
  ) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_classic() +
  labs(
    x = "Condition",
    y = "Percentage (%)",
    fill = "Cell Type"
  )

figS4B

ggsave(
  filename = "FigureS4B_Epithelial_celltype_distribution.pdf",
  plot = figS4B,
  device = "pdf",
  width = 14,
  height = 12,
  units = "cm"
)


# =============================================================================
# Figure S4C
# Cell-cycle phase distribution in epithelial clusters
# =============================================================================

df_cycle <- test.seu4@meta.data %>%
  filter(condition %in% c("Health", "EoE")) %>%
  group_by(celltype2, condition, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(celltype2, condition) %>%
  mutate(perc = n / sum(n) * 100) %>%
  ungroup()

offset <- 0.2

df_cycle <- df_cycle %>%
  mutate(
    celltype2 = factor(
      celltype2,
      levels = c("Basal", "Cycling Suprabasal", "Suprabasal", "Apical")
    ),
    Phase = factor(
      Phase,
      levels = c("G1", "S", "G2M")
    ),
    x = as.numeric(celltype2) +
      ifelse(condition == "Health", -offset, offset)
  )

df_cycle_label <- df_cycle %>%
  group_by(celltype2, condition) %>%
  summarise(
    x = unique(x),
    y_top = sum(perc),
    .groups = "drop"
  )

figS4C <- ggplot(
  df_cycle,
  aes(x = x, y = perc, fill = Phase)
) +
  geom_bar(
    stat = "identity",
    position = position_stack(reverse = TRUE),
    width = 0.3,
    color = "black"
  ) +
  geom_text(
    data = df_cycle_label,
    aes(x = x, y = y_top + 1.5, label = condition),
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_x_continuous(
    breaks = seq_along(levels(df_cycle$celltype2)),
    labels = levels(df_cycle$celltype2)
  ) +
  scale_fill_manual(
    values = color_cc,
    breaks = c("G1", "S", "G2M")
  ) +
  theme_classic() +
  labs(
    x = "Cell Type",
    y = "Percentage (%)",
    fill = "Phase"
  )

figS4C

ggsave(
  filename = "FigureS4C_CellCycle_Health_vs_EoE.pdf",
  plot = figS4C,
  device = "pdf",
  dpi = 300,
  width = 20,
  height = 15,
  units = "cm"
)


# =============================================================================
# Figure S4D
# IL13, IL33, and IL36B expression in healthy and EoE samples
# =============================================================================

figS4D_IL13 <- plot_umap4_facet2(
  seurat_obj = test.seu4,
  gene_name = "IL13",
  title_size = 7
)

figS4D_IL33 <- plot_umap4_facet2(
  seurat_obj = test.seu4,
  gene_name = "IL33",
  title_size = 7
)

figS4D_IL36B <- plot_umap4_facet2(
  seurat_obj = test.seu4,
  gene_name = "IL36B",
  title_size = 7
)

figS4D <- figS4D_IL13 / figS4D_IL33 / figS4D_IL36B

figS4D

ggsave(
  filename = "FigureS4D_IL13_IL33_IL36B_UMAP.pdf",
  plot = figS4D,
  device = "pdf",
  width = 12,
  height = 18,
  units = "cm"
)


# =============================================================================
# Figure S4E
# IL1A and IL1B expression in healthy and EoE samples
# =============================================================================

figS4E_IL1A <- plot_umap4_facet2(
  seurat_obj = test.seu4,
  gene_name = "IL1A",
  title_size = 7
)

figS4E_IL1B <- plot_umap4_facet2(
  seurat_obj = test.seu4,
  gene_name = "IL1B",
  title_size = 7
)

figS4E <- figS4E_IL1A / figS4E_IL1B

figS4E

ggsave(
  filename = "FigureS4E_IL1A_IL1B_UMAP.pdf",
  plot = figS4E,
  device = "pdf",
  width = 12,
  height = 12,
  units = "cm"
)


# =============================================================================
# Figure S4F
# IL1A expression in apical cells, shown separately for each donor
# =============================================================================

# Load the processed object used by the final IL1A plotting code in IL1a.R.
# The preprocessing code from IL1a.R is intentionally not included here.
test.seu_il1a <- qread("E:/EoE/[2026](05)EoE_afterQC_+Remi_reduce_harmony.qs")

#
# This section is taken from the plotting portion of IL1a.R.

apical <- subset(
  test.seu_il1a,
  subset = celltype == "Apical cell" &
    condition %in% c("Health", "EoE", "Remission")
)

apical_df <- FetchData(
  apical,
  vars = c("IL1A", "patient", "condition"),
  slot = "data"
) %>%
  rownames_to_column("cell") %>%
  filter(!is.na(patient), !is.na(condition))

apical_df$condition <- factor(
  apical_df$condition,
  levels = c("Health", "EoE", "Remission")
)

patient_order <- apical_df %>%
  distinct(patient, condition) %>%
  arrange(condition, patient) %>%
  pull(patient)

apical_df$patient <- factor(
  apical_df$patient,
  levels = patient_order
)

group_pos <- apical_df %>%
  distinct(patient, condition) %>%
  mutate(x = as.numeric(patient)) %>%
  group_by(condition) %>%
  summarise(
    xmin = min(x) - 0.45,
    xmax = max(x) + 0.45,
    xmid = mean(x),
    .groups = "drop"
  )

ymax <- max(apical_df$IL1A, na.rm = TRUE)
ylim_max <- ymax * 1.35

figS4F <- ggplot(
  apical_df,
  aes(x = patient, y = IL1A, fill = condition)
) +
  geom_violin(
    trim = FALSE,
    scale = "width",
    adjust = 0.75,
    width = 0.9,
    color = "grey50",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = cond_cols,
    breaks = c("Health", "EoE", "Remission"),
    labels = c("Healthy", "EoE", "Remission")
  ) +
  scale_y_continuous(
    limits = c(0, ylim_max),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = expression(italic("IL1A") ~ expression)
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text.x = element_text(
      angle = 50,
      hjust = 1,
      vjust = 1,
      color = "black",
      size = 9
    ),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.y = element_text(color = "black", size = 12),
    axis.title.x = element_blank(),
    plot.margin = margin(5, 5, 25, 5)
  ) +
  coord_cartesian(clip = "off") +
  geom_segment(
    data = group_pos,
    aes(x = xmin, xend = xmax, y = -0.35, yend = -0.35),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_text(
    data = group_pos,
    aes(x = xmid, y = -0.55, label = condition),
    inherit.aes = FALSE,
    size = 4
  )

figS4F

ggsave(
  filename = "FigureS4F_IL1A_Apical_per_patient.pdf",
  plot = figS4F,
  device = "pdf",
  width = 6.2,
  height = 3.2
)


# =============================================================================
# Figure S4G
# IL1 receptor-family expression across epithelial populations
# =============================================================================

genes <- c("IL1R1", "IL1R2", "IL1RAP", "IL1RN")

DefaultAssay(test.seu4) <- "RNA"

celltypes <- c(
  "Basal",
  "Cycling Suprabasal",
  "Suprabasal",
  "Apical"
)

conditions <- c("Health", "EoE")

group_levels <- as.vector(
  sapply(
    celltypes,
    function(ct) paste0(ct, "_", conditions)
  )
)

expr_mat <- GetAssayData(
  test.seu4,
  assay = "RNA",
  layer = "data"
)[genes, ]

expr_df <- as.data.frame(t(expr_mat))

meta <- test.seu4@meta.data[, c("celltype2", "condition")]

df_long <- expr_df %>%
  mutate(
    celltype2 = meta$celltype2,
    condition = meta$condition,
    cellcond = factor(
      paste0(celltype2, "_", condition),
      levels = group_levels
    )
  ) %>%
  pivot_longer(
    cols = all_of(genes),
    names_to = "gene",
    values_to = "expr"
  )

df_summary <- df_long %>%
  group_by(cellcond, gene) %>%
  summarise(
    avg_expr = mean(expr, na.rm = TRUE),
    pct_expr = mean(expr > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(gene) %>%
  mutate(
    norm_expr = scale(avg_expr)[, 1]
  ) %>%
  ungroup()

figS4G <- ggplot(
  df_summary,
  aes(
    x = cellcond,
    y = gene,
    color = norm_expr,
    size = pct_expr
  )
) +
  geom_point() +
  scale_color_gradient(
    low = "gray97",
    high = "#A90C38",
    name = "Scaled\nexpr.",
    limits = c(-1, 2),
    oob = squish
  ) +
  scale_size(
    range = c(2, 10),
    name = "% expr. cells",
    labels = percent
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

figS4G

ggsave(
  filename = "FigureS4G_IL1_receptors_DotPlot.pdf",
  plot = figS4G,
  device = "pdf",
  width = 6,
  height = 5
)
