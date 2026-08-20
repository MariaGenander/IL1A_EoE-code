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
library(clustree)
library(patchwork)
library(ggmagnify)
library(ggforce)
library(harmony)
library(nichenetr)
library(Nebulosa)
library(ggridges)
library(presto)
library(UCell)
library(clusterProfiler)
library(Matrix)
library(viridis)
library(readr)
library(dplyr)
library(stringr)
library(openxlsx)
library(ComplexHeatmap)
library(DESeq2)
library(circlize)
library(ComplexHeatmap)
library(monocle3)
library(Seurat)
library(ggplot2)
library(dplyr)
library(grDevices)

library(CellChat)
library(NMF)
library(ggalluvial)



#devtools::install_github("saeyslab/nichenetr")
#devtools::install_github("saeyslab/multinichenetr")

library(nichenetr)
library(multinichenetr)
packageVersion("multinichenetr")
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(muscat)


color_scheme <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

color_cc <- c("#E15759", "#F0C420", "#1965B0")

# Read Seurat object (this one is V5 so later we need to convert to v4)
sce <- readRDS("E:\\EoE\\0925_EoE_afterQC_noRemi_reduce_harmony_REcelltype.RDS")

head(sce)
unique(sce$condition)

# Load built-in S and G2M gene sets
cc.genes <- Seurat::cc.genes
sce <- CellCycleScoring(sce, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

gd3 <- DimPlot(sce, group.by = "Phase", reduction = "umap", 
               pt.size = 0.2, cols = color_cc, label = TRUE, label.size = 5)
gd3

gd2 <- DimPlot(object = sce, reduction = "umap", group.by = "simplified_celltype",
               label = TRUE, cols = color_scheme)
gd2


#### Subdividing Basal cells
# First, copy the simplified_celltype column
sce$celltype2 <- sce$simplified_celltype

# Modify based on conditions from Phase and simplified_celltype
sce$celltype2[sce$simplified_celltype == "Basal" & sce$Phase == "G1"] <- "Quiescent Basal"
sce$celltype2[sce$simplified_celltype == "Basal" & sce$Phase %in% c("S", "G2M")] <- "Cycling Basal"

# View results
table(sce$simplified_celltype)
table(sce$celltype2)




organism = "human"

options(timeout = 120)

if(organism == "human"){
  
  lr_network_all = 
    readRDS(url(
      "https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds"
    )) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(url(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"
  ))
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
} else if(organism == "mouse"){
  
  lr_network_all = readRDS(url(
    "https://zenodo.org/record/10229222/files/lr_network_mouse_allInfo_30112033.rds"
  )) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(url(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"
  ))
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
}



# SCE transform
options(future.globals.maxSize = 50 * 1024^3)  # Set the maximum size to 50 GiB

## Seurat V5 to V4 conversion
sce[["RNA3"]] <- as(object = sce[["RNA"]], Class = "Assay")
DefaultAssay(sce) <- 'RNA3'
## You need to remove the original RNA assay before renaming RNA3 to RNA, otherwise merging will fail later
sce[['RNA']]<-NULL
sce <- RenameAssays(sce,assay.name = 'RNA3',new.assay.name = 'RNA')


sce <- as.SingleCellExperiment(sce)

sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()


colData(sce)




###### START #######
# adjust celltype to fit R language
table(sce[[celltype_id]])
colData(sce)[[celltype_id]] <- make.names(colData(sce)[[celltype_id]])
table(sce[[celltype_id]])

table(sce[[celltype_id]], sce[[group_id]])
table(sce[[celltype_id]], sce[[sample_id]])

### define group
sample_id = "patient"
group_id = "condition"
celltype_id = "celltype2"

covariates = NA
batches = NA

### Define the contrasts of interest
contrasts_oi = c("'EoE-Health','Health-EoE'") 

contrast_tbl = tibble(
  contrast = c("EoE-Health","Health-EoE"), 
  group = c("EoE","Health")
)

#If you want to focus the analysis on specific cell types (e.g. because you know which cell types 
#reside in the same microenvironments based on spatial data), you can define this here. 
#If you have sufficient computational resources and no specific idea of cell-type colocalzations, 
#we recommend to consider all cell types as potential senders and receivers.

#Here we will consider all cell types in the data:

senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
]

#
conditions_keep = c("EoE", "Health")
sce = sce[, SummarizedExperiment::colData(sce)[,group_id] %in% 
            conditions_keep
]



#Parameters for step 1: Cell-type filtering:
min_cells = 10

#Parameters for step 2: Gene filtering
min_sample_prop = 0.50

#But how do we define which genes are expressed in a sample? 
#For this we will consider genes as expressed if they have non-zero expression 
#values in a fraction_cutoff fraction of cells of that cell type in that sample. By default, 
#we set fraction_cutoff = 0.05, 
#which means that genes should show non-zero expression values in at least 5% of cells in a sample.

fraction_cutoff = 0.05

#Parameters for step 4: DE analysis
empirical_pval = FALSE

#Parameters for step 5: Ligand activity prediction
logFC_threshold = 0.50
p_val_threshold = 0.05
p_val_adj = FALSE 

top_n_target = 250
n.cores = 1

#Parameters for step 6: Prioritization
scenario = "regular"
ligand_activity_down = FALSE



######---Running the MultiNicheNet wrapper function--######

#

multinichenet_output = multi_nichenet_analysis(
  sce = sce, 
  celltype_id = celltype_id, sample_id = sample_id, group_id = group_id, 
  batches = batches, covariates = covariates, 
  lr_network = lr_network, ligand_target_matrix = ligand_target_matrix, 
  contrasts_oi = contrasts_oi, contrast_tbl = contrast_tbl, 
  senders_oi = senders_oi, receivers_oi = receivers_oi,
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, 
  min_sample_prop = min_sample_prop,
  scenario = scenario, 
  ligand_activity_down = ligand_activity_down,
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj, 
  empirical_pval = empirical_pval, 
  top_n_target = top_n_target, 
  n.cores = n.cores, 
  verbose = TRUE
)












######Downstream analysis of the MultiNicheNet output######

#Normalized pseudobulk expression for each cell type - sample combination
multinichenet_output$celltype_info$pb_df %>% head()

multinichenet_output$celltype_info$pb_df_group %>% head()


#DE information for each cell type - contrast combination
multinichenet_output$celltype_de %>% head()


#Output of the NicheNet ligand activity analysis, and the NicheNet ligand-target inference
multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>% head()


#Tables with the final prioritization scores (results per group and per sample)
multinichenet_output$prioritization_tables$group_prioritization_tbl %>% head()





######Visualization of differential cell-cell interactions######
#
prioritized_tbl_oi_all = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  top_n = 50, 
  rank_per_group = FALSE
)

prioritized_tbl_oi = 
  multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% 
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

#colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
#colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)


colors_sender = RColorBrewer::brewer.pal(n = min(12, length(senders_receivers)), name = 'Set3') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Set3') %>% magrittr::set_names(senders_receivers)


circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)




#Interpretable bubble plots

group_oi = "EoE"

prioritized_tbl_oi_M_50 = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  top_n = 100, 
  groups_oi = group_oi)

plot_oi = make_sample_lr_prod_activity_plots(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_M_50)

plot_oi


#As a further help for further prioritization, we can assess the level of 
#curation of these LR pairs as defined by the Intercellular Communication 
#part of the Omnipath database

prioritized_tbl_oi_M_50_omnipath = prioritized_tbl_oi_M_50 %>% 
  inner_join(lr_network_all)

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_M_50_omnipath)

plot_oi


ggsave(filename = '1.pdf', plot = plot_oi, device = 'pdf', 
       width = 50, height = 80, units = 'cm')





#Quiescent.Basal as receiver:
prioritized_tbl_oi_M_50 = get_top_n_lr_pairs(
multinichenet_output$prioritization_tables, 
200, 
groups_oi = group_oi, 
receivers_oi = "Fibroblast")

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_M_50 %>% inner_join(lr_network_all))

#plot_oi

ggsave(filename = '(Top200)[Receiver]Fibroblast.pdf', plot = plot_oi, device = 'pdf', 
       width = 40, height = 60, units = 'cm')






#Fibro as receiver:
prioritized_tbl_oi_M_50 = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  50, 
  groups_oi = group_oi, 
  senders_oi = "Fibroblast")

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_M_50 %>% inner_join(lr_network_all))

#plot_oi

ggsave(filename = '[Sender]Fibroblast.pdf', plot = plot_oi, device = 'pdf', 
       width = 50, height = 80, units = 'cm')





#### Visualization of differential ligand-target links
# Without filtering of target genes based on LR-target expression correlation

group_oi = "EoE"
receiver_oi = "Quiescent.Basal"
prioritized_tbl_oi_M_10 = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  10, 
  groups_oi = group_oi, 
  receivers_oi = receiver_oi)

combined_plot = make_ligand_activity_target_plot(
  group_oi, 
  receiver_oi, 
  prioritized_tbl_oi_M_10,
  multinichenet_output$prioritization_tables, 
  multinichenet_output$ligand_activities_targets_DEgenes, contrast_tbl, 
  multinichenet_output$grouping_tbl, 
  multinichenet_output$celltype_info, 
  ligand_target_matrix, 
  plot_legend = FALSE)

combined_plot
combined_plot$combined_plot

design <- "
ABCCC
DDDDD
EEEEE"

p2 <- combined_plot$combined_plot + combined_plot$legends + plot_layout(design = design)
p2

ggsave(filename = '[Diff][Receive]Quiescent.Basal.pdf', 
       plot = p2, device = 'pdf', 
       width = 85, height = 35, units = 'cm')



##With filtering of target genes based on LR-target expression correlation

group_oi = "EoE"
receiver_oi = "Quiescent.Basal"
lr_target_prior_cor_filtered = multinichenet_output$lr_target_prior_cor %>%
  inner_join(
    multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>% 
      distinct(ligand, target, direction_regulation, contrast)
  ) %>% 
  inner_join(contrast_tbl) %>% filter(group == group_oi, receiver == receiver_oi)

lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>% 
  filter(direction_regulation == "up") %>% 
  filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>% 
  filter(direction_regulation == "down") %>% 
  filter( (rank_of_target < top_n_target) & (pearson < -0.50 | spearman < -0.50)) # downregulation -- negative correlation
lr_target_prior_cor_filtered = bind_rows(
  lr_target_prior_cor_filtered_up, 
  lr_target_prior_cor_filtered_down)

prioritized_tbl_oi = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  50, 
  groups_oi = group_oi, 
  receivers_oi = receiver_oi)

lr_target_correlation_plot = make_lr_target_correlation_plot(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi,  
  lr_target_prior_cor_filtered , 
  multinichenet_output$grouping_tbl, 
  multinichenet_output$celltype_info, 
  receiver_oi,
  plot_legend = FALSE)

lr_target_correlation_plot$combined_plot


design <- "
AABBDD
CCCCDD"

p3 <- lr_target_correlation_plot$combined_plot + lr_target_correlation_plot$legends + 
  plot_layout(design = design)
p3

ggsave(filename = '[Diff][Corre][Receive]Quiescent.Basal.pdf', 
       plot = p3, device = 'pdf', 
       width = 85, height = 35, units = 'cm')




## You can also visualize the expression correlation 
# in the following way for a selected LR pair and their targets:

ligand_oi = "IL1B"
receptor_oi = "IL1R2"
sender_oi = "Langerhans.cell"
receiver_oi = "Quiescent.Basal"
lr_target_scatter_plot = make_lr_target_scatter_plot(
  multinichenet_output$prioritization_tables, 
  ligand_oi, receptor_oi, sender_oi, receiver_oi, 
  multinichenet_output$celltype_info, 
  multinichenet_output$grouping_tbl, 
  lr_target_prior_cor_filtered)
lr_target_scatter_plot







###
##
#### Intercellular regulatory network inference and visualization
#


#
prioritized_tbl_oi = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  100, 
  rank_per_group = FALSE)

lr_target_prior_cor_filtered = 
  multinichenet_output$prioritization_tables$group_prioritization_tbl$group %>% unique() %>% 
  lapply(function(group_oi){
    lr_target_prior_cor_filtered = multinichenet_output$lr_target_prior_cor %>%
      inner_join(
        multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>%
          distinct(ligand, target, direction_regulation, contrast)
      ) %>% 
      inner_join(contrast_tbl) %>% filter(group == group_oi)
    
    lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>% 
      filter(direction_regulation == "up") %>% 
      filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
    
    lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>% 
      filter(direction_regulation == "down") %>% 
      filter( (rank_of_target < top_n_target) & (pearson < -0.50 | spearman < -0.50))
    lr_target_prior_cor_filtered = bind_rows(
      lr_target_prior_cor_filtered_up, 
      lr_target_prior_cor_filtered_down
    )
  }) %>% bind_rows()

lr_target_df = lr_target_prior_cor_filtered %>% 
  distinct(group, sender, receiver, ligand, receptor, id, target, direction_regulation) 


network = infer_intercellular_regulatory_network(lr_target_df, prioritized_tbl_oi)
network$links %>% head()

network$nodes %>% head()


#colors_sender["L_T_TIM3._CD38._HLADR."] = "pink" # the  original yellow background with white font is not very readable
network_graph = visualize_network(network, colors_sender)
network_graph$plot


###
#Interestingly, we can also use this network to further prioritize
#differential CCC interactions. Here we will assume that the most important 
#LR interactions are the ones that are involved
#in this intercellular regulatory network. We can get these interactions as follows:

network$prioritized_lr_interactions

prioritized_tbl_oi_network = prioritized_tbl_oi %>% inner_join(
  network$prioritized_lr_interactions)
prioritized_tbl_oi_network

group_oi = "EoE"

prioritized_tbl_oi_M = prioritized_tbl_oi_network %>% filter(group == group_oi)

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_M %>% inner_join(lr_network_all)
)
plot_oi

ggsave(filename = 'prioritized_lr.pdf', plot = plot_oi, device = 'pdf', 
       width = 40, height = 35, units = 'cm')

##






#

#