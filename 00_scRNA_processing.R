setwd("E:\\EoE")

.libPaths(c("~/SeuratV4", .libPaths()))

library(Seurat)
packageVersion("Seurat")
library(DoubletFinder)
library(cowplot)
library(cols4all)
library(ggplot2)
library(dplyr)
library(patchwork)
library(harmony)
library(Matrix)


# =============================================================================
# Color settings used for QC and dimensional-reduction plots
# =============================================================================

color_scheme <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", "#FF7F00", "#FDB462",
                  "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02",
                  "#7570B3", "#BEAED4", "#666666", "#999999", "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3",
                  "#808000", "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00", "#F0F8FF", "#FAEBD7",
                  "#00FFFF", "#7FFFD4", "#F0FFFF", "#F5F5DC", "#FFE4C4", "#000000", "#FFEBCD", "#0000FF",
                  "#8A2BE2", "#A52A2A", "#DEB887", "#5F9EA0", "#7FFF00", "#D2691E", "#FF7F50", "#6495ED",
                  "#FFF8DC", "#DC143C", "#00FFFF", "#00008B", "#008B8B", "#B8860B", "#A9A9A9", "#006400",
                  "#A9A9A9", "#BDB76B", "#8B008B", "#556B2F")

color_cc <- c("#E15759", "#F0C420", "#1965B0")
color_ccm <- c("gray75", "#b40001")


# =============================================================================
# 1. Load the count matrix and create the Seurat object
# =============================================================================

# Load the decompressed barcode, feature, and matrix files.
barcodes <- read.table("E:/EoE/raw/raw/EoE_cell.tsv", header = FALSE)
features <- read.table("E:/EoE/raw/raw/EoE_gene.tsv", header = FALSE)
matrix <- readMM("E:/EoE/raw/raw/EoE.mtx")

# Assign gene and cell names to the sparse matrix.
rownames(matrix) <- features$V1
colnames(matrix) <- barcodes$V1

# Create the Seurat object.
test.seu <- CreateSeuratObject(counts = matrix)

head(test.seu@meta.data)
unique(test.seu$orig.ident)
print(ncol(test.seu))


# =============================================================================
# 2. Add clinical and reference metadata
# =============================================================================

# Extract the patient/sample identifier before ":" from orig.ident.
test.seu$patient <- sapply(strsplit(as.character(test.seu$orig.ident), ":"), `[`, 1)
head(test.seu$patient)

# Create the condition column.
test.seu$condition <- NA

# Define samples for each clinical condition.
health_samples <- c("E1904", "E1664", "E1611", "E1503", "E1343", "E1280", "E1178")
eoe_samples <- c("E1839", "E1712", "E1542", "E1218", "E1154", "E1100", "E1054", "E1036")
remission_samples <- c("E1939", "E1881", "E1729", "E1172", "E1158", "E1102", "E1021")

# Assign condition labels.
test.seu$condition[test.seu$patient %in% health_samples] <- "Health"
test.seu$condition[test.seu$patient %in% eoe_samples] <- "EoE"
test.seu$condition[test.seu$patient %in% remission_samples] <- "Remission"

head(test.seu@meta.data)


# Add sex and cell-type metadata.
metadata <- read.table("./raw/EoE_meta2.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(metadata)
print(nrow(metadata))
head(metadata$NAME)
head(rownames(test.seu@meta.data))

rownames(metadata) <- metadata$NAME
metadata <- metadata[, c("sex", "celltype")]
test.seu <- AddMetaData(object = test.seu, metadata = metadata)

head(test.seu@meta.data)


# Add the original 2D coordinates as metadata.
metadata <- read.table("./raw/EoE_coord_2d.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(metadata)
print(nrow(metadata))
head(metadata$NAME)
head(rownames(test.seu@meta.data))

rownames(metadata) <- metadata$NAME
metadata <- metadata[, c("X", "Y")]
test.seu <- AddMetaData(object = test.seu, metadata = metadata)

head(test.seu@meta.data)
unique(test.seu$celltype)
ncol(test.seu)
nrow(test.seu)


# =============================================================================
# 3. Downsample to 100,000 cells while retaining small cell-type groups
# =============================================================================

meta <- test.seu@meta.data
meta$group <- meta$celltype

target_total <- 100000

# Count cells in each group and process smaller groups first.
group_summary <- meta %>%
  group_by(group) %>%
  summarise(n_cells = n()) %>%
  arrange(n_cells)

sampled_cells <- c()
remaining_target <- target_total
remaining_groups <- nrow(group_summary)

# Retain all cells from small groups and randomly sample larger groups.
for (i in 1:nrow(group_summary)) {
  grp <- group_summary$group[i]
  cells_in_grp <- rownames(meta)[meta$group == grp]
  current_group_size <- length(cells_in_grp)

  target_per_group <- floor(remaining_target / remaining_groups)

  if (current_group_size <= target_per_group) {
    sampled_cells <- c(sampled_cells, cells_in_grp)
    remaining_target <- remaining_target - current_group_size
  } else {
    sampled <- sample(cells_in_grp, target_per_group)
    sampled_cells <- c(sampled_cells, sampled)
    remaining_target <- remaining_target - target_per_group
  }

  remaining_groups <- remaining_groups - 1
}

length(sampled_cells)

# Subset the Seurat object using the sampled cell barcodes.
test.seu.sub <- subset(test.seu, cells = sampled_cells)

head(test.seu.sub@meta.data)
unique(test.seu.sub$celltype)
table(test.seu.sub$celltype)
ncol(test.seu.sub)
nrow(test.seu.sub)

saveRDS(test.seu, file = './[2025](01)EoE_raw.RDS')
saveRDS(test.seu.sub, file = './[2025](02)EoE_raw_100000.RDS')


# =============================================================================
# 4. Quality control
# =============================================================================

test.seu <- readRDS('./[2025](01)EoE_raw.RDS')
test.seu.sub <- readRDS('./[2025](02)EoE_raw_100000.RDS')

test.seu <- test.seu.sub
rm(test.seu.sub)

# Calculate mitochondrial and ribosomal read percentages.
test.seu[["percent.mt"]] <- PercentageFeatureSet(test.seu, pattern = "^MT-")
test.seu[["percent.rp"]] <- PercentageFeatureSet(test.seu, pattern = "^RP[SL]")

# Visualize QC metrics by condition and cell type.
cor1 <- c("#D24744", "#C53341", "#B71D3E", "#4D7DA9", "#3D6B98", "#2E5A87")

plotQC <- VlnPlot(
  test.seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"),
  ncol = 2,
  group.by = "condition",
  pt.size = 0,
  cols = cor1
)
plotQC

plotQC2 <- VlnPlot(
  test.seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"),
  ncol = 2,
  group.by = "celltype",
  pt.size = 0,
  cols = color_scheme
)
plotQC2

# Apply the QC thresholds used in the original analysis.
test.seu2 <- subset(
  test.seu,
  subset = nFeature_RNA > 50 &
    percent.mt < 45
)

table(test.seu$celltype)
print(ncol(test.seu))
table(test.seu2$celltype)
print(ncol(test.seu2))

test.seu <- test.seu2
rm(test.seu2)

print(ncol(test.seu))

# Remove remission samples for the Health vs EoE analysis.
test.seu <- subset(test.seu, subset = condition != "Remission")

print(ncol(test.seu))
table(test.seu$celltype)

saveRDS(test.seu, file = './[2025](03)EoE_raw_100000_NoRemission.RDS')


# =============================================================================
# 5. Collapse the original annotations into broader cell-type categories
# =============================================================================

table(test.seu$celltype)
head(test.seu@meta.data)

test.seu@meta.data <- test.seu@meta.data %>%
  mutate(celltype2 = case_when(
    # Keep epithelial populations as separate groups.
    celltype == "Apical cell" ~ "Superficial",
    celltype == "Basal cell (cycling)" ~ "Cycling Basal",
    celltype == "Quiescent basal cell" ~ "Quiescent Basal ",
    celltype == "Suprabasal" ~ "Suprabasal",

    # T-cell populations.
    grepl("T cell|Tcm|Trm|Teff|Temra|Th|Treg|Gamma delta|MAIT", celltype) ~ "T cell",

    # B-cell and plasma-cell populations.
    grepl("B cell", celltype) ~ "B cell",
    grepl("Plasma cell", celltype) ~ "Plasma cell",

    # Myeloid populations.
    grepl("Monocyte", celltype) ~ "Monocyte",
    grepl("Macrophage", celltype) ~ "Macrophage",
    grepl("cDC|DC|pDC|Langerhans", celltype) ~ "Dendritic cell",

    # Other immune populations.
    grepl("NK", celltype) ~ "NK cell",
    grepl("Eosinophil", celltype) ~ "Eosinophil",
    grepl("Mast cell", celltype) ~ "Mast cell",
    grepl("ILC", celltype) ~ "Lymphoid Cell",

    # Stromal and vascular populations.
    grepl("Arterial|Capillary|BEC|Venous|LEC", celltype) ~ "Endothelial cell",
    grepl("Fibroblast", celltype) ~ "Fibroblast",
    grepl("Pericyte", celltype) ~ "Pericyte",

    # Other populations.
    grepl("Glia", celltype) ~ "Glia",
    grepl("Erythroid", celltype) ~ "Erythroid cell",
    TRUE ~ "Other"
  ))

head(test.seu@meta.data[, c("celltype", "celltype2")])
table(test.seu$celltype)
table(test.seu$celltype2)

saveRDS(test.seu, file = './[2025](03)EoE_raw_100000_NoRemission_REcelltype.RDS')


# =============================================================================
# 6. Doublet detection with DoubletFinder
# =============================================================================

test.seu <- readRDS('./[2025](03)EoE_raw_100000_NoRemission_REcelltype.RDS')

# Standard preprocessing before DoubletFinder.
test.seu <- NormalizeData(test.seu, normalization.method = "LogNormalize", scale.factor = 10000)
test.seu <- FindVariableFeatures(test.seu, selection.method = "vst", nfeatures = 2000)
test.seu <- ScaleData(test.seu, features = VariableFeatures(test.seu))
test.seu <- RunPCA(test.seu, features = VariableFeatures(test.seu), npcs = 50)
test.seu <- FindNeighbors(test.seu, dims = 1:20)
test.seu <- FindClusters(test.seu, resolution = 0.5)
test.seu <- RunUMAP(test.seu, dims = 1:20)
test.seu <- RunTSNE(test.seu, dims = 1:20)

# Identify the optimal pK value.
sweep.res.list <- paramSweep(test.seu, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pk_v <- as.numeric(as.character(bcmvn$pK))
pk_good <- pk_v[bcmvn$BCmetric == max(bcmvn$BCmetric)]

# Estimate the expected number of doublets as 10% of all cells.
nExp_poi <- round(0.1 * length(colnames(test.seu)))

test.seu <- doubletFinder(
  test.seu,
  PCs = 1:10,
  pN = 0.25,
  pK = pk_good,
  nExp = nExp_poi,
  reuse.pANN = FALSE,
  sct = FALSE
)

colnames(test.seu@meta.data)[ncol(test.seu@meta.data)] <- "DoubletFinder"

DF_df <- test.seu@meta.data[, "DoubletFinder", drop = FALSE]
write.table(
  DF_df,
  file = "DoubletFinder_result.txt",
  quote = FALSE,
  sep = '\t',
  row.names = TRUE,
  col.names = TRUE
)

DimPlot(test.seu, reduction = "tsne", pt.size = 1, group.by = "DoubletFinder")

saveRDS(test.seu, file = './[2025](03)EoE_raw_100000_NoRemission_REcelltype_doublet.RDS')

# Retain singlets only.
test.seu2 <- subset(test.seu, subset = (DoubletFinder == "Singlet"))

DimPlot(test.seu, reduction = "tsne", pt.size = 1, group.by = "DoubletFinder") +
  DimPlot(test.seu2, reduction = "tsne", pt.size = 1, group.by = "DoubletFinder")

saveRDS(test.seu2, file = './[2025](03)EoE_raw_100000_NoRemission_REcelltype_singlet.RDS')

test.seu <- test.seu2


# Save an additional singlet object used in the original workflow.
test.seu <- readRDS('./[2025](03)EoE_raw_100000_NoRemission_REcelltype_doublet.RDS')
test.seu2 <- subset(test.seu, subset = (DoubletFinder == "Singlet"))

table(test.seu2$celltype2)

test.seu$DoubletFinder <- factor(test.seu$DoubletFinder, levels = c("Singlet", "Doublet"))
unique(test.seu$condition)
test.seu$condition <- factor(test.seu$condition, levels = c("Health", "EoE"))

DimPlot(
  test.seu,
  reduction = "umap",
  pt.size = 0.5,
  cols = color_ccm,
  group.by = "DoubletFinder"
) +
  DimPlot(
    test.seu2,
    reduction = "umap",
    pt.size = 0.5,
    cols = color_ccm,
    group.by = "DoubletFinder"
  )

ggsave(
  filename = "umap_doublet.pdf",
  device = 'pdf',
  width = 50,
  height = 20,
  units = 'cm'
)

saveRDS(test.seu2, file = './[2025](04)EoE_raw_100000_NoRemission_REcelltype_singlet.RDS')

test.seu <- test.seu2
rm(test.seu2)


# =============================================================================
# 7. Normalization, PCA, Harmony integration, UMAP, and t-SNE
# =============================================================================

test.seu <- readRDS('./[2025](03)EoE_raw_100000_NoRemission_REcelltype_singlet.RDS')

test.seu <- NormalizeData(test.seu)
test.seu <- FindVariableFeatures(test.seu, selection.method = "vst", nfeatures = 2000)
test.seu <- ScaleData(test.seu, features = VariableFeatures(test.seu))
test.seu <- RunPCA(
  test.seu,
  features = VariableFeatures(object = test.seu),
  npcs = 50,
  verbose = TRUE
)

# Generate uncorrected embeddings.
test.seu <- RunUMAP(test.seu, dims = 1:30, reduction.name = "umap_naive")
test.seu <- RunTSNE(test.seu, dims = 1:30, reduction.name = "tsne_naive")

# Correct patient-associated batch effects with Harmony.
test.seu <- RunHarmony(
  test.seu,
  reduction = "pca",
  group.by.vars = "patient",
  reduction.save = "harmony"
)

test.seu <- RunUMAP(test.seu, reduction = "harmony", dims = 1:30, reduction.name = "umap")
test.seu <- RunTSNE(test.seu, reduction = "harmony", dims = 1:30, reduction.name = "tsne")

# Compare embeddings before and after Harmony correction.
p1 <- DimPlot(test.seu, reduction = "umap_naive", group.by = "condition", cols = color_cc)
p2 <- DimPlot(test.seu, reduction = "umap", group.by = "condition", cols = color_cc)
p4 <- DimPlot(test.seu, reduction = "tsne_naive", group.by = "condition", cols = color_cc)
p5 <- DimPlot(test.seu, reduction = "tsne", group.by = "condition", cols = color_cc)

p7 <- CombinePlots(plots = list(p1, p2, p4, p5), ncol = 2)
p7

ggsave(
  filename = "[2025]EoE_umap&tsne_harmony.pdf",
  plot = p7,
  device = 'pdf',
  width = 50,
  height = 40,
  units = 'cm'
)

saveRDS(test.seu, file = '[2025](05)EoE_afterQC_noRemi_reduce_harmony.RDS')


# =============================================================================
# 8. Cell-cycle scoring and clustering-parameter exploration
# =============================================================================

# Score cell-cycle phase using Seurat's built-in gene sets.
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

test.seu <- CellCycleScoring(
  test.seu,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = TRUE
)

# Explore different numbers of Harmony dimensions at resolution 0.7.
dims_list <- c(20, 22, 25, 27, 30)
resolutions_list <- c(0.7)

for (dims_val in dims_list) {
  for (resolution_val in resolutions_list) {

    test.seu <- FindNeighbors(
      object = test.seu,
      reduction = "harmony",
      dims = 1:dims_val,
      verbose = TRUE
    )

    test.seu <- FindClusters(
      object = test.seu,
      resolution = resolution_val,
      verbose = TRUE
    )

    test.seu <- RunUMAP(
      object = test.seu,
      reduction = "harmony",
      dims = 1:dims_val,
      verbose = TRUE,
      reduction.name = "umap"
    )

    gd2 <- DimPlot(
      object = test.seu,
      reduction = "umap",
      group.by = 'celltype2',
      label = TRUE,
      label.size = 3,
      cols = color_scheme
    )

    filename <- paste('[2025]2EoE_test_', dims_val, '_', resolution_val, '.pdf', sep = '')
    ggsave(filename = filename, plot = gd2, width = 12, height = 9)
  }
}


# Explore UMAP n.neighbors and min.dist using 25 Harmony dimensions.
dims_list <- c(25)
resolutions_list <- c(0.7)
nneighbors_list <- c(25, 30)
min_dist_list <- c(0.4, 0.25)

for (dims_val in dims_list) {
  for (resolution_val in resolutions_list) {
    for (n_neighbor in nneighbors_list) {
      for (min_dist in min_dist_list) {

        test.seu <- FindNeighbors(
          object = test.seu,
          reduction = "harmony",
          dims = 1:dims_val,
          verbose = TRUE
        )

        test.seu <- FindClusters(
          object = test.seu,
          resolution = resolution_val,
          verbose = TRUE
        )

        test.seu <- RunUMAP(
          object = test.seu,
          reduction = "harmony",
          dims = 1:dims_val,
          n.neighbors = n_neighbor,
          min.dist = min_dist,
          verbose = TRUE,
          reduction.name = "umap"
        )

        gd2 <- DimPlot(
          object = test.seu,
          reduction = "umap",
          group.by = 'celltype2',
          label = TRUE,
          label.size = 3,
          cols = color_scheme
        )

        filename <- paste(
          "[2025]EoE_test_",
          dims_val,
          "_",
          resolution_val,
          "{",
          n_neighbor,
          "_",
          min_dist,
          "}.pdf",
          sep = ""
        )

        ggsave(filename = filename, plot = gd2, width = 12, height = 9)
      }
    }
  }
}


# =============================================================================
# 9. Final clustering and dimensional reduction
# =============================================================================

# Final settings selected in the original analysis.
test.seu <- FindNeighbors(test.seu, dims = 1:25, reduction = "harmony")
test.seu <- FindClusters(test.seu, resolution = 0.7)

test.seu <- RunUMAP(
  object = test.seu,
  reduction = "harmony",
  dims = 1:25,
  reduction.name = "umap",
  n.neighbors = 30,
  min.dist = 0.3
)

test.seu <- RunTSNE(
  object = test.seu,
  reduction = "harmony",
  dims = 1:25,
  reduction.name = "tsne"
)

saveRDS(test.seu, file = './0701_EoE_afterQC_noRemi_reduce_harmony_25_0.7-20-0.2.RDS')


# =============================================================================
# 10. Final cell-type regrouping
# =============================================================================

table(test.seu$celltype2)

# Remove leading/trailing whitespace from the broad cell-type labels.
test.seu$celltype2 <- trimws(test.seu$celltype2)

# Remove glial cells.
test.seu <- subset(test.seu, subset = celltype2 != "Glia")

# Rename selected epithelial populations.
test.seu$celltype2[test.seu$celltype2 == "Cycling Basal"] <- "Cycling Suprabasal"
test.seu$celltype2[test.seu$celltype2 == "Quiescent Basal"] <- "Basal"
test.seu$celltype2[test.seu$celltype2 == "Superficial"] <- "Apical"

table(test.seu$celltype2)

# Define broad epithelial, immune, and stromal categories.
epithelial <- c("Cycling Suprabasal", "Basal", "Apical", "Suprabasal")
immune <- c(
  "B cell", "Dendritic cell", "Lymphoid Cell",
  "Macrophage", "Mast cell", "Monocyte",
  "NK cell", "Plasma cell", "T cell"
)
stromal <- c("Endothelial cell", "Pericyte", "Fibroblast")

test.seu$celltype3 <- ifelse(
  test.seu$celltype2 %in% epithelial,
  "Epithelial",
  ifelse(
    test.seu$celltype2 %in% immune,
    "Immune",
    ifelse(
      test.seu$celltype2 %in% stromal,
      "Stromal",
      ifelse(
        test.seu$celltype2 == "Eosinophil",
        "Eosinophil",
        ifelse(test.seu$celltype2 == "Erythroid cell", "Erythroid", NA)
      )
    )
  )
)

# Keep endothelial cells, pericytes, and fibroblasts as separate stromal groups.
test.seu$celltype4 <- ifelse(
  test.seu$celltype2 %in% epithelial,
  "Epithelial",
  ifelse(
    test.seu$celltype2 %in% immune,
    "Immune",
    ifelse(
      test.seu$celltype2 == "Endothelial cell",
      "Endothelial",
      ifelse(
        test.seu$celltype2 %in% c("Pericyte", "Fibroblast"),
        test.seu$celltype2,
        ifelse(
          test.seu$celltype2 == "Eosinophil",
          "Eosinophil",
          ifelse(test.seu$celltype2 == "Erythroid cell", "Erythroid", NA)
        )
      )
    )
  )
)

# Save the final processed Seurat object.
saveRDS(
  test.seu,
  file = './[2025](06.1)EoE_raw_100000_NoRemission_REcelltype_25-0.7-30-0.3.RDS'
)
