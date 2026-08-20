# =============================================================================
# Figure 4 and Figure S4
# Intraepithelial IL1 signaling in EoE
# =============================================================================

.libPaths(c("~/SeuratV4", .libPaths()))

library(Seurat)
library(qs)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(ggpubr)
library(rstatix)
library(scales)
library(circlize)
library(multinichenetr)
library(nichenetr)
library(RColorBrewer)

setwd("E:\\2025Plot")
source("E:/2025Plot/EoE_source.R")

setwd("E:\\2025Plot\\IL1A paper")


# -----------------------------------------------------------------------------
# Single-cell objects used for Figure 4 / Figure S4
# -----------------------------------------------------------------------------

test.seu <- readRDS("E:/2025Plot/IL1A paper/[paper]EoE_all.RDS")

# Epithelial-only object used for the single-cell panels.
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


# =============================================================================
# Figure 4A
# UMAP of esophageal epithelial cells from healthy and EoE samples
# =============================================================================

fig4A <- plot_umap1(
  test.seu4,
  group.by = "celltype2",
  reduction = "umap",
  label = FALSE,
  cols = color_scheme
)

fig4A

ggsave(
  filename = "Figure4A_Epithelial_UMAP.pdf",
  plot = fig4A,
  device = "pdf",
  dpi = "retina",
  width = 14,
  height = 12,
  units = "cm"
)


# =============================================================================
# Figure 4E-F
# MultiNicheNet: top 20 ligand-receptor interactions directed toward basal cells
# =============================================================================
#
# The MultiNicheNet analysis itself is contained in multinichenetr.R.
# Run/source that script first so that `multinichenet_output` and
# `lr_network_all` are present in the R environment.
#
# source("E:/EoE/multinichenetr.R")


# Cell-type colors used for the circos plots.
cci_colors <- c(
  "#EFE2AA", "#FDB462", "#EF9A9A", "#E15759",
  "#CB95BB", "#FFCDD2", "#E78AC3", "#FFEB3B",
  "#DBC9B3", "#7EC0C3", "#BBDEFB", "#8ECFF8",
  "#1E88E5", "#83B4EF", "#1965B0", "#A5D6A7",
  "#AED0DF", "#64B5F6", "#89B780", "#EED0E0"
)

make_basal_circos <- function(group_name, output_file) {

  top_ids <- get_top_n_lr_pairs(
    multinichenet_output$prioritization_tables,
    top_n = 20,
    groups_oi = group_name,
    receivers_oi = "Basal"
  )$id

  tbl_all <- multinichenet_output$prioritization_tables$group_prioritization_tbl

  pri_tbl <- tbl_all %>%
    filter(id %in% top_ids) %>%
    group_by(group, sender, receiver) %>%
    summarise(
      prioritization_score = sum(prioritization_score),
      .groups = "drop"
    ) %>%
    mutate(
      ligand = sender,
      receptor = receiver,
      id = paste(group, sender, receiver, sep = "_")
    ) %>%
    filter(receiver == "Basal", group == group_name)

  celltypes <- sort(unique(c(pri_tbl$sender, pri_tbl$receiver)))

  if (length(celltypes) > length(cci_colors)) {
    plot_colors <- colorRampPalette(cci_colors)(length(celltypes))
  } else {
    plot_colors <- cci_colors[seq_along(celltypes)]
  }

  cols <- setNames(plot_colors, celltypes)

  pdf(output_file, width = 8, height = 8)
  circos.clear()

  make_circos_one_group(
    prioritized_tbl_oi = pri_tbl,
    colors_sender = cols,
    colors_receiver = cols
  )

  dev.off()
}


# Figure 4E: healthy esophagus -> basal-cell communication.
make_basal_circos(
  group_name = "Health",
  output_file = "Figure4E_MultiNiche_Health_to_Basal_top20.pdf"
)


# Figure 4F: EoE esophagus -> basal-cell communication.
make_basal_circos(
  group_name = "EoE",
  output_file = "Figure4F_MultiNiche_EoE_to_Basal_top20.pdf"
)


# =============================================================================
# Figure 4G
# MultiNicheNet bubble plot of the top increased ligand-receptor interactions
# directed toward basal cells in EoE
# =============================================================================

group_oi <- "EoE"

fig4G_lr <- get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables,
  top_n = 20,
  groups_oi = group_oi,
  receivers_oi = "Basal"
)

# Retain the OmniPath-annotated interactions, matching the original workflow.
fig4G_lr_omnipath <- fig4G_lr %>%
  inner_join(lr_network_all)

fig4G <- make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables,
  fig4G_lr_omnipath
)

fig4G

ggsave(
  filename = "Figure4G_MultiNiche_top20_to_Basal.pdf",
  plot = fig4G,
  device = "pdf",
  width = 40,
  height = 18,
  units = "cm"
)


# =============================================================================
# Figure 4H
# Protein-localization panel
# =============================================================================
#
# This panel is microscopy-based in the manuscript figure set.
# No R-generated image-processing panel corresponding to this microscopy image
# was identified in the supplied plotting scripts.
# The microscopy image itself should therefore be assembled from the original
# image-analysis/export workflow rather than recreated from the Seurat object.


# =============================================================================
# Figure 4I
# IL1A expression across healthy, active EoE, and remission epithelium
# =============================================================================
#
# This section is taken from the plotting part of IL1a.R.
# The preprocessing portion of IL1a.R is intentionally not reproduced here.

test.seu_il1a <- qread("E:/EoE/[2026](05)EoE_afterQC_+Remi_reduce_harmony.qs")

epithelial_types <- c(
  "Apical cell",
  "Quiescent basal cell",
  "Suprabasal",
  "Basal cell (cycling)"
)

epi.seu <- subset(
  test.seu_il1a,
  subset = celltype %in% epithelial_types
)

epi.seu$condition <- factor(
  epi.seu$condition,
  levels = c("Health", "EoE", "Remission")
)

cond_cols <- c(
  "Health" = "gray75",
  "EoE" = "#b40001",
  "Remission" = "#3D6B98"
)

cell_df <- FetchData(
  epi.seu,
  vars = c("IL1A", "condition", "patient", "celltype")
) %>%
  rownames_to_column("cell") %>%
  filter(!is.na(condition)) %>%
  mutate(
    condition = factor(
      condition,
      levels = c("Health", "EoE", "Remission")
    )
  )

my_comparisons <- list(
  c("Health", "EoE"),
  c("EoE", "Remission")
)

cell_stat <- cell_df %>%
  wilcox_test(
    IL1A ~ condition,
    comparisons = my_comparisons
  ) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  filter(p.adj <= 0.05)

if (nrow(cell_stat) > 0) {
  ymax1 <- max(cell_df$IL1A, na.rm = TRUE)

  cell_stat$y.position <- seq(
    from = ymax1 * 1.08,
    to = ymax1 * 1.20,
    length.out = nrow(cell_stat)
  )
}

fig4I <- ggplot(
  cell_df,
  aes(x = condition, y = IL1A, fill = condition)
) +
  geom_violin(
    trim = FALSE,
    scale = "width",
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(values = cond_cols) +
  labs(
    x = NULL,
    y = expression(italic("IL1A") ~ Expression)
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "none"
  )

if (nrow(cell_stat) > 0) {
  fig4I <- fig4I +
    stat_pvalue_manual(
      cell_stat,
      label = "p.adj.signif",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.4,
      size = 5,
      inherit.aes = FALSE
    )
}

fig4I

ggsave(
  filename = "Figure4I_IL1A_Epithelium_Health_EoE_Remission.pdf",
  plot = fig4I,
  device = "pdf",
  width = 5,
  height = 5
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
