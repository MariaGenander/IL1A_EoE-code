# =============================================================================
# 12_BulkRNA_05_GSEA_GO.R
# GSEA and over-representation enrichment analysis
# =============================================================================

library(clusterProfiler)
library(GseaVis)
library(DOSE)
library(ggplot2)
library(dplyr)
library(cowplot)

setwd("D:\\Eve_Bulk\\R\\Pairwise")

# Gene-set database used in the original analysis.
kegg_gmt <- read.gmt("v2024.1.Hs.symbols.gmt")


# =============================================================================
# Helper functions
# =============================================================================

run_gsea <- function(result_file, output_csv) {

  res_df <- read.table(
    result_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  res_df <- res_df %>%
    filter(!is.na(log2FoldChange))

  gene_list <- res_df$log2FoldChange
  names(gene_list) <- res_df$geneNames

  gene_list <- sort(gene_list, decreasing = TRUE)

  gsea <- GSEA(
    gene_list,
    TERM2GENE = kegg_gmt,
    minGSSize = 5,
    pvalueCutoff = 0.2,
    pAdjustMethod = "fdr"
  )

  gsea_df <- as.data.frame(gsea)
  gsea_df <- gsea_df[order(gsea_df$enrichmentScore, decreasing = TRUE), ]

  write.csv(
    gsea_df,
    output_csv,
    quote = FALSE,
    row.names = FALSE
  )

  gsea
}


run_ora <- function(
    result_file,
    output_csv,
    lfc_cutoff = 0.585,
    pvalue_cutoff = 0.05
) {

  res_df <- read.table(
    result_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # The original ORA code used nominal P value < 0.05 together with
  # |log2FoldChange| > 0.585.
  res_sig <- res_df %>%
    filter(
      !is.na(pvalue),
      pvalue < pvalue_cutoff,
      abs(log2FoldChange) > lfc_cutoff
    )

  ego <- enricher(
    gene = res_sig$geneNames,
    TERM2GENE = kegg_gmt,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )

  out <- as.data.frame(ego)
  out <- out[order(out$p.adjust), ]

  write.csv(
    out,
    output_csv,
    quote = FALSE,
    row.names = FALSE
  )

  ego
}


plot_selected_gsea <- function(
    gsea_object,
    terms,
    output_file,
    ncol = 3,
    width = 40,
    height = 30
) {

  available <- intersect(
    terms,
    as.data.frame(gsea_object)$ID
  )

  if (length(available) == 0) {
    warning("None of the requested GSEA terms were found.")
    return(invisible(NULL))
  }

  plot_list <- lapply(
    available,
    function(term) {
      gseaNb(
        object = gsea_object,
        geneSetID = term,
        subPlot = 2,
        addPval = TRUE,
        termWidth = 70,
        pvalX = 0.90,
        pvalY = 0.88,
        pvalSize = 6,
        kegg = FALSE,
        geneCol = "#4d4d4d",
        add.geneExpHt = FALSE,
        rmSegment = TRUE,
        geneSize = 3,
        htHeight = 0.1,
        rmPrefix = FALSE
      )
    }
  )

  p <- cowplot::plot_grid(
    plotlist = plot_list,
    ncol = ncol,
    align = "hv"
  )

  ggsave(
    output_file,
    p,
    device = cairo_pdf,
    width = width,
    height = height,
    units = "cm"
  )

  invisible(p)
}


# =============================================================================
# 1. P0 NGFR+ vs NGFR-
# =============================================================================

gsea_P0_NGFR <- run_gsea(
  "{Pairwise}[DESeq2]P0_N+_vs_N-.txt",
  "[GSEA]P0 N+ vs N-.csv"
)

terms_P0_NGFR <- c(
  "PID_INTEGRIN1_PATHWAY",
  "REACTOME_LAMININ_INTERACTIONS",
  "GOBP_CELL_FATE_COMMITMENT",
  "GOBP_POSITIVE_REGULATION_OF_NOTCH_SIGNALING_PATHWAY",
  "GOBP_REGULATION_OF_FIBROBLAST_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY",
  "KEGG_HEDGEHOG_SIGNALING_PATHWAY",
  "KEGG_MEDICUS_REFERENCE_TIGHT_JUNCTION_ACTIN_SIGNALING_PATHWAY",
  "GOBP_POSITIVE_REGULATION_OF_RETINOIC_ACID_RECEPTOR_SIGNALING_PATHWAY",
  "REACTOME_NEGATIVE_REGULATION_OF_NOTCH4_SIGNALING",
  "GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY"
)

plot_selected_gsea(
  gsea_P0_NGFR,
  terms_P0_NGFR,
  "[GSEA_plot_multi] P0 N+ vs N- .pdf",
  ncol = 4,
  width = 75,
  height = 60
)

go_P0_NGFR <- run_ora(
  "{Pairwise}[DESeq2]P0_N+_vs_N-.txt",
  "[GO 0.585]{Pvalue}P0 N+ vs N-.csv"
)


# =============================================================================
# 2. P1 control NGFR+ vs NGFR-
# =============================================================================

gsea_P1_NGFR <- run_gsea(
  "{Pairwise}[DESeq2]P1_Ctrl N- vs N+.txt",
  "[GSEA]P1_Ctrl N- vs N+.csv"
)

terms_P1_NGFR <- c(
  "HALLMARK_E2F_TARGETS",
  "GOBP_CELL_CYCLE_CHECKPOINT_SIGNALING",
  "KEGG_HEDGEHOG_SIGNALING_PATHWAY",
  "GOBP_REGULATION_OF_FIBROBLAST_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_ACTIVIN_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY",
  "REACTOME_SIGNALING_BY_NOTCH",
  "GOBP_CELL_CELL_ADHESION",
  "GOBP_CELL_CELL_JUNCTION_ASSEMBLY",
  "GOBP_REGULATION_OF_EPITHELIAL_CELL_DIFFERENTIATION"
)

plot_selected_gsea(
  gsea_P1_NGFR,
  terms_P1_NGFR,
  "[GSEA_plot_multi] P1_Ctrl N- vs N+ .pdf",
  ncol = 4,
  width = 65,
  height = 50
)

go_P1_NGFR <- run_ora(
  "{Pairwise}[DESeq2]P1_Ctrl N- vs N+.txt",
  "[GO 0.585]{Pvalue}P1_Ctrl N- vs N+.csv"
)


# =============================================================================
# 3. NGFR+ control cells: P1 vs P0
# =============================================================================

gsea_Npos_P1P0 <- run_gsea(
  "{Pairwise}[DESeq2]N+ P1 vs P0.txt",
  "[GSEA][Ctrl] N+ P0 vs P1.csv"
)

terms_Npos_P1P0 <- c(
  "GAVISH_3CA_METAPROGRAM_EPITHELIAL_CELL_CYCLE",
  "REACTOME_CHOLESTEROL_BIOSYNTHESIS",
  "REACTOME_DECTIN_1_MEDIATED_NONCANONICAL_NF_KB_SIGNALING",
  "REACTOME_NEGATIVE_REGULATION_OF_NOTCH4_SIGNALING",
  "REACTOME_INTERLEUKIN_1_PROCESSING",
  "REACTOME_ERKS_ARE_INACTIVATED",
  "GOBP_CELL_SURFACE_RECEPTOR_PROTEIN_TYROSINE_KINASE_SIGNALING_PATHWAY",
  "GOBP_INFLAMMATORY_RESPONSE",
  "GOBP_NEGATIVE_REGULATION_OF_CELL_POPULATION_PROLIFERATION",
  "GOBP_NEGATIVE_REGULATION_OF_CELL_DIFFERENTIATION",
  "GOBP_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_MORPHOGENESIS_OF_AN_EPITHELIUM"
)

plot_selected_gsea(
  gsea_Npos_P1P0,
  terms_Npos_P1P0,
  "[GSEA_plot_multi] [Ctrl] N+ P0 vs P1 .pdf",
  ncol = 4,
  width = 65,
  height = 55
)

go_Npos_P1P0 <- run_ora(
  "{Pairwise}[DESeq2]N+ P1 vs P0.txt",
  "[GO 0.585]{Pvalue}[Ctrl] N+ P0 vs P1.csv"
)


# =============================================================================
# 4. NGFR- control cells: P1 vs P0
# =============================================================================

gsea_Nneg_P1P0 <- run_gsea(
  "{Pairwise}[DESeq2]N- P1 vs P0.txt",
  "[GSEA][Ctrl] N- P1 vs P0.csv"
)

terms_Nneg_P1P0 <- c(
  "REACTOME_INTERLEUKIN_1_PROCESSING",
  "REACTOME_CELL_EXTRACELLULAR_MATRIX_INTERACTIONS",
  "REACTOME_CHOLESTEROL_BIOSYNTHESIS",
  "KEGG_STEROID_BIOSYNTHESIS",
  "GOBP_EPITHELIAL_CELL_DIFFERENTIATION",
  "GOBP_REGULATION_OF_CELL_CELL_ADHESION",
  "GOBP_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE",
  "GOBP_POSITIVE_REGULATION_OF_CELL_DEVELOPMENT",
  "REACTOME_INTERFERON_SIGNALING"
)

plot_selected_gsea(
  gsea_Nneg_P1P0,
  terms_Nneg_P1P0,
  "[GSEA_plot_multi] [Ctrl] N- P1 vs P0 .pdf",
  ncol = 3,
  width = 55,
  height = 45
)

go_Nneg_P1P0 <- run_ora(
  "{Pairwise}[DESeq2]N- P1 vs P0.txt",
  "[GO 0.585]{Pvalue}[Ctrl] N- P0 vs P1.csv"
)


# =============================================================================
# 5. P1 IL1A vs control, NGFR+ and NGFR- combined
# =============================================================================

gsea_P1_IL1A <- run_gsea(
  "{Pairwise}[DESeq2]P1 IL1A vs Ctr.txt",
  "[GSEA]P1 IL1A vs Ctr.csv"
)

# This pathway was explicitly plotted in the original workflow.
if ("HALLMARK_TNFA_SIGNALING_VIA_NFKB" %in% as.data.frame(gsea_P1_IL1A)$ID) {

  p_nfkb <- gseaNb(
    object = gsea_P1_IL1A,
    geneSetID = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    subPlot = 2,
    addPval = TRUE,
    termWidth = 70,
    pvalX = 0.95,
    pvalY = 0.9,
    pvalSize = 6,
    addGene = TRUE,
    markTopgene = TRUE,
    topGeneN = 10,
    kegg = FALSE,
    geneCol = "#4d4d4d",
    rmSegment = TRUE,
    rmPrefix = FALSE
  )

  ggsave(
    "[GSEA All IL1A]HALLMARK_TNFA_SIGNALING_VIA_NFKB.pdf",
    p_nfkb,
    device = cairo_pdf,
    width = 15,
    height = 18,
    units = "cm"
  )
}

go_P1_IL1A <- run_ora(
  "{Pairwise}[DESeq2]P1 IL1A vs Ctr.txt",
  "[GO 0.585]{Pvalue}P1 IL1A vs Ctr.csv"
)


# =============================================================================
# 6. P1 NGFR+ IL1A vs control
# =============================================================================

gsea_Npos_IL1A <- run_gsea(
  "{Pairwise}[DESeq2]P1 N+ IL1A vs Ctr.txt",
  "[GSEA]P1 N+ IL1A vs Ctr.csv"
)

terms_Npos_IL1A <- c(
  "WP_INTEGRINMEDIATED_CELL_ADHESION",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "GOBP_CELLULAR_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN",
  "GOBP_NEUTROPHIL_CHEMOTAXIS"
)

plot_selected_gsea(
  gsea_Npos_IL1A,
  terms_Npos_IL1A,
  "[GSEA_plot_multi]  P1 N+ IL1A vs Ctr .pdf",
  ncol = 2,
  width = 35,
  height = 36
)

go_Npos_IL1A <- run_ora(
  "{Pairwise}[DESeq2]P1 N+ IL1A vs Ctr.txt",
  "[GO 0.585]{Pvalue}P1 N+ IL1A vs Ctr.csv"
)


# =============================================================================
# 7. P1 NGFR- IL1A vs control
# =============================================================================

gsea_Nneg_IL1A <- run_gsea(
  "{Pairwise}[DESeq2]P1 N- IL1A vs Ctr.txt",
  "[GSEA]P1 N- IL1A vs Ctr.csv"
)

terms_Nneg_IL1A <- c(
  "GOBP_REGULATION_OF_CELL_ADHESION",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "GOBP_KERATINIZATION",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)

plot_selected_gsea(
  gsea_Nneg_IL1A,
  terms_Nneg_IL1A,
  "[GSEA_plot_multi] P1 N- IL1A vs Ctr .pdf",
  ncol = 2,
  width = 40,
  height = 36
)

go_Nneg_IL1A <- run_ora(
  "{Pairwise}[DESeq2]P1 N- IL1A vs Ctr.txt",
  "[GO 0.585]{Pvalue}P1 N- IL1A vs Ctr.csv"
)
