setwd("E:\\EoE")

.libPaths(c("~/SeuratV4", .libPaths()))
library(Seurat)
packageVersion("Seurat")

library(DoubletFinder)
library(cowplot)
library(cols4all)
library(ggplot2)
library(clustree)
library(RColorBrewer)
library(magrittr)
library(dplyr)
library(tibble)
library(patchwork)
library(ggmagnify)
library(ggforce)
library(harmony)
library(Matrix)
library(grDevices)


# =============================================================================
# Color settings
# =============================================================================

#color_scheme <- c4a("poly.wright25")
color_scheme <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

expanded_color_scheme <- colorRampPalette(color_scheme)(60)
color_scheme <- expanded_color_scheme

color_scheme <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", "#FF7F00", "#FDB462",
                  "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02",
                  "#7570B3", "#BEAED4", "#666666", "#999999", "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3",
                  "#808000", "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00", "#F0F8FF", "#FAEBD7",
                  "#00FFFF", "#7FFFD4", "#F0FFFF", "#F5F5DC", "#FFE4C4", "#000000", "#FFEBCD", "#0000FF",
                  "#8A2BE2", "#A52A2A", "#DEB887", "#5F9EA0", "#7FFF00", "#D2691E", "#FF7F50", "#6495ED",
                  "#FFF8DC", "#DC143C", "#00FFFF", "#00008B", "#008B8B", "#B8860B", "#A9A9A9", "#006400",
                  "#A9A9A9", "#BDB76B", "#8B008B", "#556B2F")

color_scheme_kelly <- c("#F14343", "#F68D41", "#F9C342", "#F9D740", "#B8B83B", "#78B44D",
                        "#5CBC54", "#33B260", "#2CB974", "#3AC0A0", "#4DCEC7", "#76D1CE",
                        "#94BFCF", "#9C94AC", "#D78EA8", "#D681B3", "#E1739B", "#DA4564",
                        "#A9263B", "#CA4659", "#CE8972", "#BFA76D", "#BFC56C", "#A1CC5B",
                        "#70D053", "#51CE60", "#6ACB8B", "#6CB8C2", "#3BA7B1", "#81A3AC")

color_scheme_circus <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6",
                         "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A",
                         "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02", "#7570B3", "#BEAED4",
                         "#666666", "#999999")

color_scheme_bear <- c("#D52B2A", "#E9492A", "#F0682A", "#F59535", "#F9BB45",
                       "#E1C64A", "#B5C44D", "#86C64E", "#60C45A", "#49C279",
                       "#3AB09E", "#31B4B4", "#2A94A1", "#688CA7", "#AA93AE",
                       "#D48FAD", "#DC75A6", "#DF578C", "#CB2E68", "#A72245",
                       "#C65A72", "#B96D6D", "#B69B6A", "#C0B14A", "#9CB840",
                       "#6AB65C", "#50BB7F", "#53A09B", "#5B92AA", "#7B8B9E")

color_scheme_kelly2 <- c("#3AC0A0", "#2CB974", "#B8B83B", "#D681B3", "#94BFCF", "#F9D740",
                         "#A9263B", "#6ACB8B", "#CE8972", "#CA4659", "#A1CC5B", "#70D053",
                         "#BFC56C", "#DA4564", "#9C94AC", "#4DCEC7", "#78B44D", "#81A3AC",
                         "#6CB8C2", "#3BA7B1", "#F9C342", "#5CBC54", "#F68D41", "#D78EA8",
                         "#E1739B", "#F14343", "#76D1CE", "#33B260", "#51CE60", "#BFA76D")

color_scheme_circus2 <- c("#999999", "#A6761D", "#B17BA6", "#33A02C", "#BEAED4", "#E6AB02",
                          "#E7298A", "#FDB462", "#E78AC3", "#DC050C", "#FF7F00", "#B2DF8A",
                          "#8DD3C7", "#7570B3", "#55A1B1", "#FB8072", "#882E72", "#7BAFDE",
                          "#1965B0", "#666666")

color_scheme_bear2 <- c("#B5C44D", "#E9492A", "#7B8B9E", "#53A09B", "#5B92AA", "#D52B2A",
                        "#50BB7F", "#3AB09E", "#C65A72", "#AA93AE", "#86C64E", "#C0B14A",
                        "#9CB840", "#6AB65C", "#B69B6A", "#F0682A", "#E1C64A", "#B96D6D",
                        "#F9BB45", "#A72245", "#2A94A1", "#DC75A6", "#49C279", "#688CA7",
                        "#D48FAD", "#DF578C", "#CB2E68", "#F59535", "#60C45A", "#31B4B4")

color_cc <- c("#E15759", "#F0C420", "#1965B0")
color_ccm <- c('gray75', '#b40001')
color_ccm2 <- c('#3D6B98','#F0F0F0','#C0203C')
color.gradient <- colorRampPalette(colors = c("#208421","#f0e9db","#b40001"))


# =============================================================================
# 1. Load raw data and create the Seurat object
# =============================================================================

# Load the decompressed barcode, gene, and count matrix files.
barcodes <- read.table("E:/EoE/raw/raw/EoE_cell.tsv", header = FALSE)
features <- read.table("E:/EoE/raw/raw/EoE_gene.tsv", header = FALSE)
matrix <- readMM("E:/EoE/raw/raw/EoE.mtx")

# Assign gene names and cell barcodes to the sparse matrix.
library(Matrix)
rownames(matrix) <- features$V1
colnames(matrix) <- barcodes$V1

# Create the Seurat object.
test.seu <- CreateSeuratObject(counts = matrix)

head(test.seu@meta.data)
unique(test.seu$orig.ident)
print(ncol(test.seu))


# =============================================================================
# 2. Add clinical metadata and reference coordinates
# =============================================================================

# Extract the patient identifier from orig.ident.
test.seu$patient <- sapply(strsplit(as.character(test.seu$orig.ident), ":"), `[`, 1)
head(test.seu$patient)

# Initialize the condition column.
test.seu$condition <- NA

# Define samples for each clinical condition.
health_samples <- c("E1904", "E1664", "E1611", "E1503", "E1343", "E1280", "E1178")
eoe_samples <- c("E1839", "E1712", "E1542", "E1218", "E1154", "E1100", "E1054", "E1036")
remission_samples <- c("E1939", "E1881", "E1729", "E1172", "E1158", "E1102", "E1021")

# Assign condition labels to each sample.
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


# Add reference 2D coordinates.
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
# 3. Downsample to 100,000 cells
# =============================================================================

library(dplyr)
library(Seurat)

# Use the original cell type as the grouping variable for balanced sampling.
meta <- test.seu@meta.data
meta$group <- meta$celltype

# Set the target number of cells.
target_total <- 100000

# Count cells in each group and sort groups from smallest to largest.
group_summary <- meta %>%
  group_by(group) %>%
  summarise(n_cells = n()) %>%
  arrange(n_cells)

sampled_cells <- c()
remaining_target <- target_total
remaining_groups <- nrow(group_summary)

# Keep all cells in smaller groups and sample larger groups as needed.
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

# Subset the Seurat object using the sampled cells.
test.seu.sub <- subset(test.seu, cells = sampled_cells)

head(test.seu.sub@meta.data)
unique(test.seu.sub$celltype)
table(test.seu.sub$celltype)

ncol(test.seu.sub)
nrow(test.seu.sub)

saveRDS(test.seu, file= './[2025](01)EoE_raw.RDS')
saveRDS(test.seu.sub, file= './[2025](02)EoE_raw_100000.RDS')


# =============================================================================
# 4. Quality control
# =============================================================================

test.seu <- readRDS('./[2025](01)EoE_raw.RDS')
test.seu.sub <- readRDS('./[2025](02)EoE_raw_100000.RDS')

test.seu <- test.seu.sub
rm(test.seu.sub)

# Calculate mitochondrial and ribosomal transcript percentages.
test.seu[["percent.mt"]] <- PercentageFeatureSet(test.seu, pattern = "^MT-")
test.seu[["percent.rp"]] <- PercentageFeatureSet(test.seu, pattern = "^RP[SL]")

# Visualize QC metrics by condition.
cor1 <- c("#D24744", "#C53341", "#B71D3E", "#4D7DA9", "#3D6B98", "#2E5A87")
plotQC <- VlnPlot(test.seu, features = c("nFeature_RNA", "nCount_RNA",
                                         "percent.mt", "percent.rp"),
                  ncol = 2,
                  group.by = "condition",
                  pt.size = 0,
                  cols = cor1)
plotQC

# Visualize QC metrics by original cell type.
plotQC2 <- VlnPlot(test.seu, features = c("nFeature_RNA", "nCount_RNA",
                                          "percent.mt", "percent.rp"),
                   ncol = 2,
                   group.by = "celltype",
                   pt.size = 0,
                   cols = color_scheme)
plotQC2

# This plot is retained from the original workflow.
plotQC3 <- VlnPlot(test.seu, features = c("nFeature_RNA", "nCount_RNA",
                                          "percent.mt", "percent.rp"),
                   ncol = 2,
                   group.by = "celltype2",
                   pt.size = 0,
                   cols = color_scheme)
plotQC3

# Visualize QC metrics for each patient.
plotQC4 <- VlnPlot(test.seu, features = c("nFeature_RNA", "nCount_RNA",
                                          "percent.mt", "percent.rp"),
                   ncol = 2,
                   group.by = "patient",
                   pt.size = 0,
                   cols = color_scheme_bear2)
plotQC4


# Apply the QC thresholds used in the analysis.
test.seu2 <- subset(test.seu,
                    subset = nFeature_RNA > 50 &
                      percent.mt < 45)
                      #percent.rp < 30 &
                      #nCount_RNA > 500)
#nCount_RNA < 100000)

table(test.seu$celltype)
print(ncol(test.seu))

table(test.seu2$celltype)
print(ncol(test.seu2))

test.seu <- test.seu2
rm(test.seu2)

print(ncol(test.seu))

#saveRDS(test.seu, file= './0701_EoE_afterQC.RDS')


# Plot the reference 2D coordinates colored by cell type.
plot_data <- data.frame(
  X = test.seu@meta.data$X,
  Y = test.seu@meta.data$Y,
  celltype = test.seu@meta.data$celltype
)

head(plot_data)

p <- ggplot(plot_data, aes(x = X, y = Y, color = celltype)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = color_scheme) +
  labs(title = "Cell Type Distribution", x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal()

p <- p +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave("celltype_sphere.pdf", plot = p, device = "pdf", width = 25, height = 9)


# Remove remission samples.
test.seu <- subset(test.seu, subset = condition != "Remission")

print(ncol(test.seu))
table(test.seu$celltype)

saveRDS(test.seu, file= './[2025](03)EoE_raw_100000_NoRemission.RDS')


# Plot the reference 2D coordinates after removing remission samples.
plot_data <- data.frame(
  X = test.seu@meta.data$X,
  Y = test.seu@meta.data$Y,
  celltype = test.seu@meta.data$celltype
)

p <- ggplot(plot_data, aes(x = X, y = Y, color = celltype)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = color_scheme) +
  labs(title = "Cell Type Distribution", x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal()

p <- p +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave("celltype_sphere_noRemi2.pdf", plot = p, device = "pdf", width = 25, height = 9)


# =============================================================================
# 5. Merge detailed cell-type annotations into broader categories
# =============================================================================

table(test.seu$celltype)
head(test.seu@meta.data)

library(dplyr)

test.seu@meta.data <- test.seu@meta.data %>%
  mutate(celltype2 = case_when(
    # Keep epithelial populations as separate categories.
    celltype == "Apical cell" ~ "Superficial",
    celltype == "Basal cell (cycling)" ~ "Cycling Basal",
    celltype == "Quiescent basal cell" ~ "Quiescent Basal ",
    celltype == "Suprabasal" ~ "Suprabasal",

    # T-cell populations.
    grepl("T cell|Tcm|Trm|Teff|Temra|Th|Treg|Gamma delta|MAIT", celltype) ~ "T cell",

    # B-cell populations.
    grepl("B cell", celltype) ~ "B cell",

    # Plasma cells.
    grepl("Plasma cell", celltype) ~ "Plasma cell",

    # Monocytes.
    grepl("Monocyte", celltype) ~ "Monocyte",

    # Macrophages.
    grepl("Macrophage", celltype) ~ "Macrophage",

    # Dendritic-cell populations, including cDC, DC, pDC, and Langerhans cells.
    grepl("cDC|DC|pDC|Langerhans", celltype) ~ "Dendritic cell",

    # NK cells.
    grepl("NK", celltype) ~ "NK cell",

    # Endothelial populations, including arterial, capillary, BEC, venous, and LEC cells.
    grepl("Arterial|Capillary|BEC|Venous|LEC", celltype) ~ "Endothelial cell",

    # Fibroblasts.
    grepl("Fibroblast", celltype) ~ "Fibroblast",

    # Pericytes.
    grepl("Pericyte", celltype) ~ "Pericyte",

    # Eosinophils.
    grepl("Eosinophil", celltype) ~ "Eosinophil",

    # Mast cells.
    grepl("Mast cell", celltype) ~ "Mast cell",

    # Innate lymphoid cells.
    grepl("ILC", celltype) ~ "Lymphoid Cell",

    # Glial cells.
    grepl("Glia", celltype) ~ "Glia",

    # Erythroid cells.
    grepl("Erythroid", celltype) ~ "Erythroid cell",

    # Assign all remaining populations to Other.
    TRUE ~ "Other"
  ))

head(test.seu@meta.data[, c("celltype", "celltype2")])

table(test.seu$celltype)
table(test.seu$celltype2)

saveRDS(test.seu, file= './[2025](03)EoE_raw_100000_NoRemission_REcelltype.RDS')


# =============================================================================
# 6. Doublet removal
# =============================================================================

test.seu <- readRDS('./[2025](03)EoE_raw_100000_NoRemission_REcelltype.RDS')

# Standard Seurat preprocessing before DoubletFinder.
test.seu <- NormalizeData(test.seu, normalization.method = "LogNormalize", scale.factor = 10000)
test.seu <- FindVariableFeatures(test.seu, selection.method = "vst", nfeatures = 2000)
test.seu <- ScaleData(test.seu, features = VariableFeatures(test.seu))
test.seu <- RunPCA(test.seu, features = VariableFeatures(test.seu),npcs = 50)
test.seu <- FindNeighbors(test.seu, dims = 1:20)
test.seu <- FindClusters(test.seu, resolution = 0.5)
test.seu <- RunUMAP(test.seu, dims = 1:20)
test.seu <- RunTSNE(test.seu, dims = 1:20)

color_scheme <- c4a("poly.wright25")
p1 <- DimPlot(test.seu, reduction = "tsne", cols = color_scheme) +
  ggtitle("t-SNE Plot raw")
p1

p3 <- DimPlot(test.seu, reduction = "umap", split.by  = 'condition',
              pt.size = 0.5, cols = color_scheme, ncol = 3)
p3

DimPlot(test.seu, reduction = "umap", cols = color_scheme)+
  ggtitle("uMAP Plot raw")


# Perform the DoubletFinder parameter sweep and select the optimal pK.
sweep.res.list <- paramSweep(test.seu, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_v <- as.numeric(as.character(bcmvn$pK))
pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]

# Assume an expected doublet rate of 10%.
nExp_poi <- round(0.1*length(colnames(test.seu)))

test.seu <- doubletFinder(test.seu, PCs = 1:10, pN = 0.25, pK = pk_good,
                          nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

colnames(test.seu@meta.data)[ncol(test.seu@meta.data)]="DoubletFinder"

DF_df <- test.seu@meta.data[,"DoubletFinder", drop = FALSE]
write.table(DF_df, file = "DoubletFinder_result.txt", quote = FALSE,
            sep = '\t', row.names = TRUE, col.names = TRUE)

DimPlot(test.seu, reduction = "tsne", pt.size = 1, group.by = "DoubletFinder")

saveRDS(test.seu, file= './[2025](03)EoE_raw_100000_NoRemission_REcelltype_doublet.RDS')

# Retain singlets only.
test.seu2 <- subset(test.seu, subset = (DoubletFinder == "Singlet"))

DimPlot(test.seu, reduction = "tsne", pt.size = 1, group.by = "DoubletFinder") +
  DimPlot(test.seu2, reduction = "tsne", pt.size = 1, group.by = "DoubletFinder")

saveRDS(test.seu2, file= './[2025](03)EoE_raw_100000_NoRemission_REcelltype_singlet.RDS')

test.seu <- test.seu2


# Re-load the doublet-classified object for visualization and save the singlet object.
test.seu <- readRDS('./[2025](03)EoE_raw_100000_NoRemission_REcelltype_doublet.RDS')

test.seu2 <- subset(test.seu, subset = (DoubletFinder == "Singlet"))

table(test.seu2$celltype2)

test.seu$DoubletFinder <- factor(test.seu$DoubletFinder, levels = c("Singlet", "Doublet"))

unique(test.seu$condition)
test.seu$condition <- factor(test.seu$condition, levels = c("Health", "EoE"))

DimPlot(test.seu, reduction = "umap", pt.size = 0.5,
        cols = color_ccm, group.by = "DoubletFinder") +
  DimPlot(test.seu2, reduction = "umap", pt.size = 0.5,
          cols = color_ccm, group.by = "DoubletFinder")

ggsave(filename = "umap_doublet.pdf", device = 'pdf',
       width = 50, height = 20, units = 'cm')

saveRDS(test.seu2, file= './[2025](04)EoE_raw_100000_NoRemission_REcelltype_singlet.RDS')

test.seu <- test.seu2
rm(test.seu2)


# =============================================================================
# 7. Normalization, dimensional reduction, and Harmony integration
# =============================================================================

test.seu <- readRDS('./[2025](03)EoE_raw_100000_NoRemission_REcelltype_singlet.RDS')

test.seu <- readRDS("./[2025](03)EoE_raw_100000_NoRemission_REcelltype_singlet.RDS",
                    refhook = function(x) {
                      message("Encountered an object of class: ", class(x))
                      x
                    })

#test.seu <- JoinLayers(test.seu)
#DefaultAssay(test.seu) <- "RNA"

test.seu <- NormalizeData(test.seu)
test.seu <- FindVariableFeatures(test.seu, selection.method = "vst", nfeatures = 2000)
test.seu <- ScaleData(test.seu, features = VariableFeatures(test.seu))
test.seu <- RunPCA(test.seu, features = VariableFeatures(object = test.seu), npcs = 50, verbose = T)

# Generate embeddings before batch correction.
test.seu <- RunUMAP(test.seu, dims = 1:30, reduction.name = "umap_naive")
test.seu <- RunTSNE(test.seu, dims = 1:30, reduction.name = "tsne_naive")

# Correct patient-associated batch effects using Harmony.
test.seu <- RunHarmony(test.seu,reduction = "pca",group.by.vars = "patient",reduction.save = "harmony")
test.seu <- RunUMAP(test.seu, reduction = "harmony", dims = 1:30,reduction.name = "umap")
test.seu <- RunTSNE(test.seu, reduction = "harmony", dims = 1:30,reduction.name = "tsne")

# Compare embeddings before and after Harmony correction.
p1 <- DimPlot(test.seu, reduction = "umap_naive",group.by = "condition", cols = color_cc)
p2 <- DimPlot(test.seu, reduction = "umap",group.by = "condition", cols = color_cc)
p3 <- p1+p2

#ggsave(filename = "umap_harmony.pdf", plot = p3, device = 'pdf', width = 50, height = 20, units = 'cm')

p4 <- DimPlot(test.seu, reduction = "tsne_naive",group.by = "condition", cols = color_cc)
p5 <- DimPlot(test.seu, reduction = "tsne",group.by = "condition", cols = color_cc)
p6 <- p1+p2

#ggsave(filename = "tsne_harmony.pdf", plot = p6, device = 'pdf', width = 50, height = 20, units = 'cm')

p7 <- CombinePlots(plots = list(p1,p2,p4,p5), ncol = 2)
p7

ggsave(filename = "[2025]EoE_umap&tsne_harmony.pdf", plot = p7,
       device = 'pdf', width = 50, height = 40, units = 'cm')

saveRDS(test.seu, file= '[2025](05)EoE_afterQC_noRemi_reduce_harmony.RDS')


# Optional Seurat v5-to-v4 conversion code retained from the original workflow.
## V5 to V4
#test.seu[["RNA"]] <- as(object = test.seu[["RNA"]], Class = "Assay")

## V5 to V4 - 2
#test.seu[["RNA5"]] <- as(object = test.seu[["RNA"]], Class = "Assay5")
#names(test.seu)
#DefaultAssay(test.seu)='RNA5'
#test.seu[['RNA']]=NULL


# =============================================================================
# 8. Cell-cycle scoring and clustering parameter exploration
# =============================================================================

# Calculate cell-cycle scores.
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
test.seu <- CellCycleScoring(test.seu, s.features = s.genes,
                             g2m.features = g2m.genes, set.ident = TRUE)


# Explore different numbers of Harmony dimensions.
dims_list <- c(20,22,25,27,30)
resolutions_list <- c(0.7)

for (dims_val in dims_list) {
  for (resolution_val in resolutions_list) {

    test.seu <- FindNeighbors(object = test.seu, reduction = "harmony",
                              dims = 1:dims_val, verbose = TRUE)
    test.seu <- FindClusters(object = test.seu,
                             resolution = resolution_val, verbose = TRUE)

    test.seu <- RunUMAP(object = test.seu, reduction = "harmony",
                        dims = 1:dims_val,
                        #n.neighbors = 25, min.dist = 0.4,
                        verbose = TRUE, reduction.name = "umap")

    gd2 <- DimPlot(object = test.seu, reduction = "umap",
                   group.by = 'celltype2',
                   label = TRUE, label.size = 3, cols = color_scheme)

    #gd3 <- DimPlot(test.seu, group.by = "Phase", reduction = "umap", pt.size = 0.2,
    #                cols = color_cc, label = TRUE, label.size = 5)

    #gd4 <- DimPlot(test.seu, group.by = "condition", reduction = "umap", pt.size = 0.2,
    #               cols = color_scheme, label = F, label.size = 5)

    #p1 <- plot_density(test.seu, c("IL1R1", "IL1R2", "IL1RAP", "IL1RN", "IL1A", "IL1B"),
    #                   joint = FALSE, reduction = "umap", combine = TRUE,
    #                   size = 0.2, pal = "inferno")

    #p2 <- plot_density(test.seu, c("COL17A1", "NGFR", "TP73", "KRT5", "KRT4", "CNFN"),
    #                   joint = FALSE, reduction = "umap", combine = TRUE,
    #                   size = 0.2, pal = "magma")

    filename <- paste('[2025]2EoE_test_', dims_val, '_', resolution_val, '.pdf', sep = '')
    ggsave(filename = filename, plot = gd2, width = 12, height = 9)

    # saveRDS(test.seu, file = paste('dc_UMAP_', dims_val, '_', resolution_val, '.rds', sep = ''))
  }
}


# Explore UMAP n.neighbors and min.dist values.
dims_list <- c(25)
resolutions_list <- c(0.7)
nneighbors_list <- c(25,30)
min_dist_list <- c(0.4,0.25)

for (dims_val in dims_list) {
  for (resolution_val in resolutions_list) {
    for (n_neighbor in nneighbors_list) {
      for (min_dist in min_dist_list) {

        test.seu <- FindNeighbors(object = test.seu, reduction = "harmony",
                                  dims = 1:dims_val, verbose = TRUE)
        test.seu <- FindClusters(object = test.seu,
                                 resolution = resolution_val, verbose = TRUE)

        test.seu <- RunUMAP(object = test.seu, reduction = "harmony",
                            dims = 1:dims_val,
                            n.neighbors = n_neighbor, min.dist = min_dist,
                            verbose = TRUE, reduction.name = "umap")

        gd2 <- DimPlot(object = test.seu, reduction = "umap",
                       group.by = 'celltype2',
                       label = TRUE, label.size = 3, cols = color_scheme)

        filename <- paste("[2025]EoE_test_", dims_val, "_", resolution_val,
                          "{", n_neighbor, "_", min_dist, "}.pdf", sep = "")

        ggsave(filename = filename, plot = gd2, width = 12, height = 9)

        # saveRDS(test.seu, file = paste("dc_UMAP_", dims_val, "_", resolution_val,
        #                               "{", n_neighbor, "_", min_dist, "}.rds", sep = ""))
      }
    }
  }
}


# =============================================================================
# 9. Final clustering and dimensional reduction
# =============================================================================

# Use 25 Harmony dimensions.
test.seu <- FindNeighbors(test.seu, dims = 1:25, reduction = "harmony")

# Use resolution 0.7.
test.seu <- FindClusters(test.seu, resolution = 0.7)

# Generate the final UMAP and t-SNE embeddings.
test.seu <- RunUMAP(object = test.seu, reduction = "harmony", dims = 1:25,
                    reduction.name = "umap",
                    n.neighbors = 30, min.dist = 0.3)

test.seu <- RunTSNE(object = test.seu, reduction = "harmony", dims = 1:25,
                    reduction.name = "tsne")

saveRDS(test.seu, file= './0701_EoE_afterQC_noRemi_reduce_harmony_25_0.7-20-0.2.RDS')


# Plot the final embeddings.
p1 <- DimPlot(test.seu, reduction = "tsne", pt.size=0.5,
              cols = color_scheme, group.by = 'celltype2',
              label = T, label.size = 2) +
  ggtitle('t-SNE plot')

p1

p2 <- DimPlot(test.seu, reduction = "umap", pt.size=0.5,
              cols = color_scheme, group.by = 'celltype',
              label = T, label.size = ) +
  ggtitle('UMAP plot')

p2

ggsave(filename = "[2025]EoE_tsne2.pdf", plot = p1, device = 'pdf',
       width = 40, height = 25, units = 'cm')

ggsave(filename = "[2025]EoE_umap2.pdf", plot = p2, device = 'pdf',
       width = 45, height = 19, units = 'cm')


p2 <- DimPlot(test.seu, reduction = "umap", pt.size=0.5,
              cols = color_scheme, group.by = 'celltype2',
              label = F) +
  ggtitle('UMAP plot')

p2

ggsave(filename = "[2025]EoE_umap3.pdf", plot = p2, device = 'pdf',
       width = 25, height = 19, units = 'cm')


p3 <- DimPlot(test.seu, reduction = "umap", split.by  = 'condition',
              group.by = 'celltype2', pt.size = 0.1, cols = color_scheme)

p3

ggsave(filename = "[2025]EoE_umap_conditions1.pdf", plot = p3,
       device = 'pdf', width = 30, height = 15, units = 'cm')


# Visualize QC metrics on the final UMAP.
p1 <- FeaturePlot(test.seu, features = "nFeature_RNA",
                  reduction = "umap", cols = color_ccm) +
  ggtitle("nFeature_RNA")

p2 <- FeaturePlot(test.seu, features = "nCount_RNA",
                  reduction = "umap", cols = color_ccm) +
  ggtitle("nCount_RNA")

p3 <- FeaturePlot(test.seu, features = "percent.mt",
                  reduction = "umap", cols = color_ccm) +
  ggtitle("percent.mt")

p4 <- FeaturePlot(test.seu, features = "percent.rp",
                  reduction = "umap", cols = color_ccm) +
  ggtitle("percent.rp")

combined_plot <- (p1 | p2) / (p3 | p4)
combined_plot

ggsave(filename = "[2025]EoE_umap_QC.pdf",
       device = 'pdf', width = 35, height = 28, units = 'cm')


# Visualize cell-cycle phase on the final UMAP.
gd3 <- DimPlot(test.seu, group.by = "Phase", reduction = "umap",
               pt.size = 0.2, cols = color_cc, label = F, label.size = 5)

gd3

ggsave(filename = "[2025]EoE_umap_Phase.pdf",
       device = 'pdf', width = 25, height = 19, units = 'cm')


# =============================================================================
# 10. Final cell-type regrouping
# =============================================================================

table(test.seu$celltype2)

# Remove leading and trailing whitespace.
test.seu$celltype2 <- trimws(test.seu$celltype2)

# Remove glial cells.
test.seu <- subset(test.seu, subset = celltype2 != "Glia")

# Rename selected epithelial populations.
test.seu$celltype2[test.seu$celltype2 == "Cycling Basal"] <- "Cycling Suprabasal"
test.seu$celltype2[test.seu$celltype2 == "Quiescent Basal"] <- "Basal"
test.seu$celltype2[test.seu$celltype2 == "Superficial"] <- "Apical"

table(test.seu$celltype2)


# Define broad cell classes.
epithelial <- c("Cycling Suprabasal", "Basal", "Apical", "Suprabasal")
immune <- c("B cell", "Dendritic cell", "Lymphoid Cell",
            "Macrophage", "Mast cell", "Monocyte", "NK cell", "Plasma cell", "T cell")
stromal <- c("Endothelial cell", "Pericyte", "Fibroblast")

# Create celltype3 using broad epithelial, immune, and stromal classes.
test.seu$celltype3 <- ifelse(test.seu$celltype2 %in% epithelial, "Epithelial",
                             ifelse(test.seu$celltype2 %in% immune, "Immune",
                                    ifelse(test.seu$celltype2 %in% stromal, "Stromal",
                                           ifelse(test.seu$celltype2 == "Eosinophil", "Eosinophil",
                                                  ifelse(test.seu$celltype2 == "Erythroid cell", "Erythroid", NA)))))

p2 <- DimPlot(test.seu, reduction = "umap", pt.size=0.5,
              cols = color_scheme, group.by = 'celltype3',
              label = F) +
  ggtitle('UMAP plot')

p2

ggsave(filename = "[2025]EoE_umap5.pdf", plot = p2, device = 'pdf',
       width = 25, height = 19, units = 'cm')


# Define a second broad grouping that keeps endothelial cells, pericytes,
# and fibroblasts as separate categories.
epithelial <- c("Cycling Suprabasal", "Basal", "Apical", "Suprabasal")
immune <- c("B cell", "Dendritic cell", "Lymphoid Cell",
            "Macrophage", "Mast cell", "Monocyte", "NK cell", "Plasma cell", "T cell")

test.seu$celltype4 <- ifelse(test.seu$celltype2 %in% epithelial, "Epithelial",
                             ifelse(test.seu$celltype2 %in% immune, "Immune",
                                    ifelse(test.seu$celltype2 == "Endothelial cell", "Endothelial",
                                           ifelse(test.seu$celltype2 %in% c("Pericyte", "Fibroblast"), test.seu$celltype2,
                                                  ifelse(test.seu$celltype2 == "Eosinophil", "Eosinophil",
                                                         ifelse(test.seu$celltype2 == "Erythroid cell", "Erythroid", NA))))))

p2 <- DimPlot(test.seu, reduction = "umap", pt.size=0.5,
              cols = color_scheme, group.by = 'celltype4',
              label = F) +
  ggtitle('UMAP plot')

p2

ggsave(filename = "[2025]EoE_umap6.pdf", plot = p2, device = 'pdf',
       width = 25, height = 19, units = 'cm')


saveRDS(test.seu, file= './[2025](06.1)EoE_raw_100000_NoRemission_REcelltype_25-0.7-30-0.3.RDS')


# =============================================================================
# 11. Generate the epithelial-only object
# =============================================================================

test.seu <- readRDS('./[2025](06.1)EoE_raw_100000_NoRemission_REcelltype_25-0.7-30-0.3.RDS')

p2 <- DimPlot(test.seu, reduction = "umap", pt.size=0.5,
              cols = color_scheme, group.by = 'celltype3',
              label = F) +
  ggtitle('UMAP plot')

p2

ggsave(filename = "[2025]EoE_umap5.pdf", plot = p2, device = 'pdf',
       width = 25, height = 19, units = 'cm')


p2 <- DimPlot(test.seu, reduction = "umap", pt.size=0.5,
              cols = color_scheme, group.by = 'celltype4',
              label = F) +
  ggtitle('UMAP plot')

p2

ggsave(filename = "[2025]EoE_umap6.pdf", plot = p2, device = 'pdf',
       width = 25, height = 19, units = 'cm')


p2 <- DimPlot(test.seu, reduction = "umap", pt.size=0.5,
              cols = color_scheme, group.by = 'celltype',
              label = T) +
  ggtitle('UMAP plot')

p2

ggsave(filename = "[2025]EoE_umap_all_celltype.pdf", plot = p2,
       device = 'pdf', width = 45, height = 19, units = 'cm')

table(test.seu$celltype)


# Subset epithelial cells only.
test.seu2 <- subset(test.seu, subset = celltype3 == "Epithelial")

p2 <- DimPlot(test.seu2, reduction = "umap", pt.size=0.5,
              cols = color_scheme, group.by = 'celltype2',
              label = T) +
  ggtitle('UMAP plot')

p2


# Remove cells outside the selected epithelial UMAP region.
test.seu2 <- subset(test.seu2, subset = UMAP_1 >= -7 & UMAP_2 <= 13.5 & UMAP_2 >= 0.5)

p2 <- DimPlot(test.seu2, reduction = "umap", pt.size=0.5,
              cols = color_scheme, group.by = 'celltype2',
              label = T) +
  ggtitle('UMAP plot')

p2

saveRDS(test.seu2, file= 'E:\\EoE\\[2025](06.2-Epi)EoE_raw_100000_NoRemission_REcelltype_25-0.7-30-0.3.RDS')
