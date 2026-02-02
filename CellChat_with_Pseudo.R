rm(list = ls()); gc()

# Load required libraries - ADD MISSING LIBRARIES
suppressPackageStartupMessages({
  library(Seurat); library(qs); library(dplyr); library(ggplot2)
  library(patchwork)  # ADDED - provides wrap_plots()
  library(ggpubr)     # ADDED - provides stat_compare_means()
  library(slingshot); library(SingleCellExperiment); library(tradeSeq)
  library(CellChat); library(ComplexHeatmap); library(circlize)
  library(ggrepel); library(viridis); library(RColorBrewer)
  library(igraph); library(pheatmap); library(corrplot); library(reshape2)
})

# ================================
# SETUP & DATA LOADING
# ================================

project_name <- "deep_pseudotime_cellchat_analysis"
base_dir <- file.path(getwd(), project_name)

# Create detailed directory structure
dirs <- c("data/processed", "results/plots/pseudotime", "results/plots/cellchat", 
          "results/plots/integration", "results/tables", "results/publication_figures")
sapply(dirs, function(d) dir.create(file.path(base_dir, d), recursive = TRUE, showWarnings = FALSE))

# Load data
seurat_obj <- qread("手动注释.qs")
DefaultAssay(seurat_obj) <- "RNA"

cat("📊 Dataset: ", ncol(seurat_obj), "cells,", nrow(seurat_obj), "genes\n")
cat("Conditions:", paste(unique(seurat_obj$condition), collapse = ", "), "\n")

# Split by condition for CellChat analysis
seurat_split_list <- SplitObject(seurat_obj, split.by = "condition")
cat("✅ Created condition-specific objects:", names(seurat_split_list), "\n")

# ================================
# PART 1: DEEP PSEUDOTIME ANALYSIS
# ================================

cat("\n🔬 DEEP PSEUDOTIME ANALYSIS\n")

# Focus on epithelial malignant progression
epithelial_cells <- c("AT2", "Malignant AT2", "Club Cells")
epithelial_subset <- subset(seurat_obj, celltype %in% epithelial_cells)

cat("Epithelial cells for trajectory:", ncol(epithelial_subset), "\n")
print(table(epithelial_subset$celltype, epithelial_subset$condition))

# Prepare SCE object for slingshot
DefaultAssay(epithelial_subset) <- "SCT"
counts_matrix <- GetAssayData(epithelial_subset, layer = "counts", assay = "SCT")
data_matrix <- GetAssayData(epithelial_subset, layer = "data", assay = "SCT")

sce <- SingleCellExperiment(
  assays = list(counts = counts_matrix, logcounts = data_matrix),
  colData = epithelial_subset@meta.data,
  reducedDims = list(
    PCA = Embeddings(epithelial_subset, "pca")[, 1:30],
    UMAP = Embeddings(epithelial_subset, "umap")
  )
)

# Run slingshot on UMAP coordinates
sce <- slingshot(sce, 
                 clusterLabels = 'celltype',
                 reducedDim = 'UMAP',  
                 start.clus = "AT2",
                 approx_points = 200)

# Add pseudotime back to Seurat object
pseudotime_matrix <- slingPseudotime(sce)
epithelial_subset@meta.data$pseudotime <- pseudotime_matrix[, 1]

cat("✅ Slingshot completed on UMAP coordinates\n")
cat("Lineages found:", length(slingLineages(sce)), "\n")
for(i in seq_along(slingLineages(sce))) {
  lineage <- slingLineages(sce)[[i]]
  cat("Lineage", i, ":", paste(lineage, collapse = " → "), "\n")
}

# Fix the smooth trajectory curve creation
create_smooth_trajectory <- function(plot_data_umap) {
  ordered_data <- plot_data_umap[order(plot_data_umap$pseudotime), ]
  loess_fit_x <- loess(UMAP1 ~ pseudotime, data = ordered_data, span = 0.3)
  loess_fit_y <- loess(UMAP2 ~ pseudotime, data = ordered_data, span = 0.3)
  pseudotime_seq <- seq(min(ordered_data$pseudotime), max(ordered_data$pseudotime), length.out = 150)
  smooth_umap1 <- predict(loess_fit_x, newdata = pseudotime_seq)
  smooth_umap2 <- predict(loess_fit_y, newdata = pseudotime_seq)
  valid_idx <- !is.na(smooth_umap1) & !is.na(smooth_umap2)
  curve_df <- data.frame(UMAP1 = smooth_umap1[valid_idx], UMAP2 = smooth_umap2[valid_idx])
  return(curve_df)
}

# ================================
# PUBLICATION-QUALITY PLOTS
# ================================

# Plot 1: Cell types on UMAP
p1_celltype <- DimPlot(epithelial_subset, reduction = "umap", group.by = "celltype",
                       label = TRUE, label.size = 6, repel = TRUE, pt.size = 1.5) +
  scale_color_manual(values = c("AT2" = "#1f77b4", "Club Cells" = "#ff7f0e", 
                                "Malignant AT2" = "#d62728")) +
  ggtitle("Epithelial Cell Types in Lung Adenocarcinoma") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    panel.grid = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  labs(color = "Cell Type")

ggsave(file.path(base_dir, "results/publication_figures/Fig1A_epithelial_celltypes_UMAP.pdf"), 
       p1_celltype, width = 12, height = 9, dpi = 300)

# Plot 2: Pseudotime trajectory on UMAP
umap_coords <- Embeddings(epithelial_subset, "umap")
plot_data_umap <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  celltype = epithelial_subset$celltype,
  pseudotime = epithelial_subset$pseudotime,
  condition = epithelial_subset$condition
) %>% filter(!is.na(pseudotime))

curve_df <- create_smooth_trajectory(plot_data_umap)

p2_trajectory <- ggplot(plot_data_umap, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = pseudotime), size = 2, alpha = 0.8) +
  scale_color_viridis_c(name = "Pseudotime", option = "plasma") +
  geom_path(data = curve_df, aes(x = UMAP1, y = UMAP2), 
            color = "#FF0000", linewidth = 4, alpha = 0.9, inherit.aes = FALSE) +
  geom_path(data = curve_df, aes(x = UMAP1, y = UMAP2), 
            color = "#FFFFFF", linewidth = 2, alpha = 0.8, inherit.aes = FALSE) +
  ggtitle("Malignant Progression Trajectory in UMAP Space") +
  theme_void() +
  theme(
    title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    legend.key.size = unit(1, "cm")
  ) +
  guides(color = guide_colorbar(barwidth = 2, barheight = 10))

ggsave(file.path(base_dir, "results/publication_figures/Fig1B_trajectory_UMAP.pdf"), 
       p2_trajectory, width = 12, height = 9, dpi = 300)

# ================================
# ENHANCED TRAJECTORY ANALYSIS FOR 2 LINEAGES
# ================================

cat("\n🔬 ANALYZING 2 TRAJECTORY LINEAGES\n")
cat("Lineage 1: AT2 → Malignant AT2 (Direct malignant transformation)\n")
cat("Lineage 2: AT2 → Club Cells (Transdifferentiation pathway)\n")

# Get both pseudotime lineages
pseudotime_lin1 <- slingPseudotime(sce)[, 1]  
pseudotime_lin2 <- slingPseudotime(sce)[, 2]  

epithelial_subset@meta.data$pseudotime_lineage1 <- pseudotime_lin1
epithelial_subset@meta.data$pseudotime_lineage2 <- pseudotime_lin2

# ================================
# DEEP TRAJECTORY GENE ANALYSIS
# ================================

cat("\n🧬 DEEP TRAJECTORY GENE ANALYSIS FOR BOTH LINEAGES\n")

# Prepare for tradeSeq
counts_filtered <- GetAssayData(epithelial_subset, layer = "counts", assay = "SCT")
keep_genes <- names(which(Matrix::rowSums(counts_filtered > 0) >= ncol(counts_filtered) * 0.1))

var_genes <- intersect(VariableFeatures(seurat_obj), keep_genes)
if(length(var_genes) < 1000) {
  mean_expr <- Matrix::rowMeans(counts_filtered[keep_genes, ])
  additional_genes <- names(sort(mean_expr, decreasing = TRUE))[1:(1000 - length(var_genes))]
  trajectory_genes <- c(var_genes, setdiff(additional_genes, var_genes))
} else {
  trajectory_genes <- var_genes[1:1000]
}

cat("Genes for trajectory analysis:", length(trajectory_genes), "\n")

# Fit GAM models
pseudotime_both <- slingPseudotime(sce, na = FALSE)
cellweights_both <- slingCurveWeights(sce)

sce_filtered <- fitGAM(
  counts = counts_filtered[trajectory_genes, ],
  pseudotime = pseudotime_both,
  cellWeights = cellweights_both,
  nknots = 6,
  verbose = TRUE
)

# Gene association tests
association_test_lin1 <- associationTest(sce_filtered, lineages = TRUE, l2fc = 0.5)
association_test_lin1$gene <- rownames(association_test_lin1)
association_test_lin1 <- association_test_lin1[!is.na(association_test_lin1$pvalue_1), ]

sig_genes_lin1 <- association_test_lin1[association_test_lin1$pvalue_1 < 0.05, ]
sig_genes_lin1 <- sig_genes_lin1[order(sig_genes_lin1$pvalue_1), ]

sig_genes_lin2 <- association_test_lin1[association_test_lin1$pvalue_2 < 0.05 & 
                                          !is.na(association_test_lin1$pvalue_2), ]
sig_genes_lin2 <- sig_genes_lin2[order(sig_genes_lin2$pvalue_2), ]

pattern_test <- patternTest(sce_filtered, global = TRUE, pairwise = TRUE)
pattern_test$gene <- rownames(pattern_test)
pattern_test <- pattern_test[!is.na(pattern_test$pvalue), ]
lineage_diff_genes <- pattern_test[pattern_test$pvalue < 0.05, ]

cat("Lineage 1 genes (AT2 → Malignant AT2):", nrow(sig_genes_lin1), "\n")
cat("Lineage 2 genes (AT2 → Club Cells):", nrow(sig_genes_lin2), "\n")
cat("Differential genes between lineages:", nrow(lineage_diff_genes), "\n")

# ================================
# KIF20A SPECIFIC ANALYSIS SECTION
# ================================

cat("\n🎯 SPECIFIC GENE ANALYSIS: KIF20A\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

# Function for comprehensive single gene analysis
analyze_specific_gene <- function(gene_name, epithelial_subset, sig_genes_lin1, sig_genes_lin2, lineage_diff_genes) {
  
  cat("\n🧬 Analyzing", gene_name, "\n")
  cat(paste(rep("-", 30), collapse = ""), "\n")
  
  # Check if gene exists
  if (!gene_name %in% rownames(epithelial_subset)) {
    cat("❌ Gene", gene_name, "not found in dataset\n")
    return(NULL)
  }
  
  # Get expression data
  gene_expr <- GetAssayData(epithelial_subset, layer = "data", assay = "SCT")[gene_name, ]
  
  # Basic statistics
  n_expressing <- sum(gene_expr > 0)
  pct_expressing <- round(n_expressing / length(gene_expr) * 100, 2)
  mean_expr <- round(mean(gene_expr), 4)
  max_expr <- round(max(gene_expr), 4)
  
  cat("📊 Expression Summary:\n")
  cat("  - Expressing cells:", n_expressing, "/", length(gene_expr), "(", pct_expressing, "%)\n")
  cat("  - Mean expression:", mean_expr, "\n")
  cat("  - Max expression:", max_expr, "\n")
  
  # Expression by cell type
  expr_by_celltype <- data.frame(
    celltype = epithelial_subset$celltype,
    expression = gene_expr,
    condition = epithelial_subset$condition,
    pseudotime_lin1 = epithelial_subset$pseudotime_lineage1,
    pseudotime_lin2 = epithelial_subset$pseudotime_lineage2
  )
  
  celltype_summary <- expr_by_celltype %>%
    group_by(celltype) %>%
    summarise(
      n_cells = n(),
      n_expressing = sum(expression > 0),
      pct_expressing = round(n_expressing / n_cells * 100, 2),
      mean_expr = round(mean(expression), 4),
      median_expr = round(median(expression), 4),
      max_expr = round(max(expression), 4),
      .groups = 'drop'
    )
  
  cat("\n📈 Expression by Cell Type:\n")
  print(celltype_summary)
  
  # Check trajectory status
  trajectory_status <- list(
    in_lineage1 = gene_name %in% sig_genes_lin1$gene,
    in_lineage2 = gene_name %in% sig_genes_lin2$gene,
    is_differential = gene_name %in% lineage_diff_genes$gene
  )
  
  cat("\n🎯 Trajectory Status:\n")
  cat("  - Significant in Lineage 1 (Direct malignant):", trajectory_status$in_lineage1, "\n")
  cat("  - Significant in Lineage 2 (Transdifferentiation):", trajectory_status$in_lineage2, "\n")
  cat("  - Differential between lineages:", trajectory_status$is_differential, "\n")
  
  return(list(
    gene_name = gene_name,
    expression_data = expr_by_celltype,
    summary = celltype_summary,
    trajectory_status = trajectory_status,
    n_expressing = n_expressing,
    pct_expressing = pct_expressing
  ))
}

# Function to create comprehensive gene plots
create_gene_plots <- function(gene_analysis_result, epithelial_subset, save_prefix) {
  
  gene_name <- gene_analysis_result$gene_name
  expr_data <- gene_analysis_result$expression_data
  
  cat("🎨 Creating plots for", gene_name, "\n")
  
  # Plot 1: Expression on UMAP
  umap_coords <- Embeddings(epithelial_subset, "umap")
  umap_data <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    expression = expr_data$expression,
    celltype = expr_data$celltype,
    condition = expr_data$condition
  )
  
  p1_umap <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = expression), size = 2, alpha = 0.8) +
    scale_color_gradient(low = "lightgray", high = "red", name = "Expression") +
    ggtitle(paste(gene_name, "Expression on UMAP")) +
    theme_void() +
    theme(
      title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold")
    )
  
  # Plot 2: Expression by cell type
  p2_celltype <- ggplot(expr_data, aes(x = celltype, y = expression)) +
    geom_violin(aes(fill = celltype), alpha = 0.7) +
    geom_boxplot(width = 0.3, fill = "white", alpha = 0.8) +
    scale_fill_manual(values = c("AT2" = "#1f77b4", "Club Cells" = "#ff7f0e", 
                                 "Malignant AT2" = "#d62728")) +
    ggtitle(paste(gene_name, "Expression by Cell Type")) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(size = 14, face = "bold"),
      title = element_text(size = 16, face = "bold"),
      legend.position = "none"
    ) +
    labs(x = "Cell Type", y = "Expression") +
    stat_summary(fun = median, geom = "point", size = 3, color = "red")
  
  # Plot 3: Expression along both lineages
  expr_both_lineages <- data.frame(
    pseudotime = c(expr_data$pseudotime_lin1, expr_data$pseudotime_lin2),
    expression = rep(expr_data$expression, 2),
    lineage = rep(c("Direct Malignant", "Transdifferentiation"), each = nrow(expr_data)),
    celltype = rep(expr_data$celltype, 2)
  ) %>% filter(!is.na(pseudotime))
  
  p3_lineages <- ggplot(expr_both_lineages, aes(x = pseudotime, y = expression)) +
    geom_point(aes(color = lineage), size = 1.5, alpha = 0.6) +
    geom_smooth(aes(color = lineage), method = "gam", linewidth = 3, se = TRUE, alpha = 0.3) +
    scale_color_manual(values = c("Direct Malignant" = "#d62728", 
                                  "Transdifferentiation" = "#2ca02c")) +
    ggtitle(paste(gene_name, "Expression Along Both Trajectories")) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold")
    ) +
    labs(x = "Pseudotime", y = "Expression", color = "Trajectory")
  
  # Plot 4: Binary expression analysis
  expr_data$expressing <- expr_data$expression > 0
  binary_summary <- expr_data %>%
    group_by(celltype, condition) %>%
    summarise(
      n_total = n(),
      n_expressing = sum(expressing),
      pct_expressing = n_expressing / n_total * 100,
      .groups = 'drop'
    )
  
  p4_binary <- ggplot(binary_summary, aes(x = celltype, y = pct_expressing)) +
    geom_col(aes(fill = condition), position = position_dodge(0.8), alpha = 0.8) +
    scale_fill_viridis_d(name = "Disease Stage") +
    ggtitle(paste("Percentage of Cells Expressing", gene_name)) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(size = 14, face = "bold"),
      title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 12)
    ) +
    labs(x = "Cell Type", y = "% Expressing Cells") +
    geom_text(aes(label = paste0(n_expressing, "/", n_total)), 
              position = position_dodge(0.8), vjust = -0.5, size = 3)
  
  # Save individual plots
  ggsave(file.path(base_dir, "results/publication_figures", 
                   paste0(save_prefix, "_", gene_name, "_UMAP.pdf")), 
         p1_umap, width = 12, height = 9, dpi = 300)
  
  ggsave(file.path(base_dir, "results/publication_figures", 
                   paste0(save_prefix, "_", gene_name, "_celltype.pdf")), 
         p2_celltype, width = 10, height = 8, dpi = 300)
  
  ggsave(file.path(base_dir, "results/publication_figures", 
                   paste0(save_prefix, "_", gene_name, "_trajectories.pdf")), 
         p3_lineages, width = 12, height = 8, dpi = 300)
  
  ggsave(file.path(base_dir, "results/publication_figures", 
                   paste0(save_prefix, "_", gene_name, "_binary.pdf")), 
         p4_binary, width = 12, height = 8, dpi = 300)
  
  # Create comprehensive figure
  comprehensive_plot <- (p1_umap + p2_celltype) / (p3_lineages + p4_binary)
  
  ggsave(file.path(base_dir, "results/publication_figures", 
                   paste0(save_prefix, "_", gene_name, "_comprehensive.pdf")), 
         comprehensive_plot, width = 20, height = 14, dpi = 300)
  
  return(list(umap = p1_umap, celltype = p2_celltype, trajectories = p3_lineages, 
              binary = p4_binary, comprehensive = comprehensive_plot))
}

# ================================
# ANALYZE KIF20A AND OTHER GENES
# ================================

# Define genes of interest
genes_of_interest <- c("KIF20A", "MKI67", "PCNA", "TP53", "EGFR", "KRAS", 
                       "SCGB1A1", "SFTPC", "CDKN2A", "CCND1")

# Filter to available genes
available_genes <- intersect(genes_of_interest, rownames(epithelial_subset))
cat("Genes available for analysis:", paste(available_genes, collapse = ", "), "\n")

# Analyze each gene
gene_analysis_results <- list()
for(gene in available_genes) {
  gene_analysis_results[[gene]] <- analyze_specific_gene(
    gene, epithelial_subset, sig_genes_lin1, sig_genes_lin2, lineage_diff_genes
  )
  
  if(!is.null(gene_analysis_results[[gene]])) {
    create_gene_plots(gene_analysis_results[[gene]], epithelial_subset, "SpecificGene")
  }
}

# ================================
# SPECIAL KIF20A ANALYSIS
# ================================

if("KIF20A" %in% available_genes) {
  
  cat("\n🎯 SPECIAL FOCUS ON KIF20A\n")
  cat(paste(rep("=", 40), collapse = ""), "\n")
  
  kif20a_result <- gene_analysis_results[["KIF20A"]]
  
  cat("📚 KIF20A BIOLOGICAL BACKGROUND:\n")
  cat("  • Full name: Kinesin Family Member 20A\n")
  cat("  • Function: Mitotic kinesin essential for cytokinesis\n")
  cat("  • Role: Required for proper chromosome segregation\n")
  cat("  • Cancer relevance: Overexpressed in many cancer types\n")
  cat("  • Clinical significance: Associated with poor prognosis\n")
  cat("  • Therapeutic potential: Target for cancer therapy\n\n")
  
  # Enhanced KIF20A analysis
  kif20a_expr <- kif20a_result$expression_data$expression
  
  # Clinical interpretation
  cat("🏥 CLINICAL INTERPRETATION FOR KIF20A IN YOUR DATA:\n")
  if(kif20a_result$pct_expressing < 5) {
    cat("  ⚠️  VERY LOW EXPRESSION: Only", kif20a_result$pct_expressing, "% of cells express KIF20A\n")
    cat("  💡 Possible interpretations:\n")
    cat("     • Cells are not actively dividing (quiescent state)\n")
    cat("     • KIF20A regulation occurs post-transcriptionally\n")
    cat("     • Expression might be cell cycle-dependent\n")
    cat("     • May not be a major driver in this specific cohort\n")
  } else {
    cat("  ✅ MODERATE EXPRESSION:", kif20a_result$pct_expressing, "% of cells express KIF20A\n")
    cat("  💡 This suggests active cell division in subset of cells\n")
  }
  
  # Find co-expressed genes with KIF20A
  if(sum(kif20a_expr > 0) > 50) {
    cat("\n🔍 FINDING KIF20A CO-EXPRESSED GENES:\n")
    
    all_expr <- GetAssayData(epithelial_subset, layer = "data", assay = "SCT")
    expressing_genes <- names(which(Matrix::rowSums(all_expr > 0) >= ncol(all_expr) * 0.05))
    
    correlations <- sapply(expressing_genes, function(g) {
      if(g == "KIF20A") return(1)
      cor(kif20a_expr, all_expr[g, ], use = "complete.obs")
    })
    
    correlations <- correlations[!is.na(correlations)]
    correlations <- sort(correlations, decreasing = TRUE)
    
    top_positive_corr <- head(correlations[correlations > 0.3], 10)
    
    cat("Top positively correlated genes (r > 0.3):\n")
    for(i in 1:length(top_positive_corr)) {
      cat(sprintf("  %d. %s (r = %.3f)\n", i, names(top_positive_corr)[i], top_positive_corr[i]))
    }
    
    # Create co-expression plot
    if(length(top_positive_corr) > 3) {
      coexpr_genes <- names(head(top_positive_corr, 6))
      coexpr_data <- data.frame(
        Gene = rep(coexpr_genes, each = ncol(epithelial_subset)),
        Expression = as.numeric(t(all_expr[coexpr_genes, ])),
        CellType = rep(epithelial_subset$celltype, length(coexpr_genes)),
        Pseudotime_Lin1 = rep(epithelial_subset$pseudotime_lineage1, length(coexpr_genes))
      ) %>% filter(!is.na(Pseudotime_Lin1))
      
      p_coexpr <- ggplot(coexpr_data, aes(x = Pseudotime_Lin1, y = Expression)) +
        geom_point(aes(color = CellType), size = 1, alpha = 0.6) +
        geom_smooth(method = "gam", color = "red", se = TRUE) +
        scale_color_manual(values = c("AT2" = "#1f77b4", "Club Cells" = "#ff7f0e", 
                                      "Malignant AT2" = "#d62728")) +
        facet_wrap(~Gene, scales = "free_y", ncol = 3) +
        ggtitle("KIF20A Co-expressed Genes Along Direct Malignant Trajectory") +
        theme_minimal() +
        theme(
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          title = element_text(size = 14, face = "bold"),
          strip.text = element_text(size = 11, face = "bold")
        ) +
        labs(x = "Pseudotime (Direct Malignant)", y = "Expression")
      
      ggsave(file.path(base_dir, "results/publication_figures/SpecialFocus_KIF20A_coexpression.pdf"), 
             p_coexpr, width = 18, height = 12, dpi = 300)
    }
  }
}

# ================================
# PART 2: ADVANCED CELLCHAT ANALYSIS
# ================================

cat("\n📞 ADVANCED CELLCHAT ANALYSIS WITH SOPHISTICATED FEATURES\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# ================================
# HELPER FUNCTIONS FOR PSEUDOTIME STORAGE
# ================================

# Function to store pseudotime analysis
store_pseudotime_analysis <- function(cellchat_obj, pseudotime_cellchats) {
  attr(cellchat_obj, "pseudotime_analysis") <- pseudotime_cellchats
  return(cellchat_obj)
}

# Function to retrieve pseudotime analysis
get_pseudotime_analysis <- function(cellchat_obj) {
  return(attr(cellchat_obj, "pseudotime_analysis"))
}

# ================================
# ADVANCED CELLCHAT ANALYSIS FUNCTION
# ================================

run_advanced_analysis_safely <- function(condition, pseudotime_data = NULL) {
  
  cat("🔬 Advanced CellChat analysis for", condition, "\n")
  
  # Get condition data
  condition_data <- seurat_split_list[[condition]]
  if (is.null(condition_data)) {
    cat("❌ No data found for condition:", condition, "\n")
    return(NULL)
  }
  
  # Prepare data
  data.input <- GetAssayData(condition_data, assay = "RNA", slot = "data")
  meta <- data.frame(
    celltype = Idents(condition_data), 
    samples = condition_data$orig.ident,
    row.names = names(Idents(condition_data))
  )
  
  # Convert to factors and clean levels
  meta$celltype <- droplevels(as.factor(meta$celltype))
  meta$samples <- as.factor(meta$samples)
  
  cat("  Cell types included:", paste(levels(meta$celltype), collapse = ", "), "\n")
  cat("  Total cells:", nrow(meta), "\n")
  
  tryCatch({
    # Create main CellChat object
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
    cellchat <- addMeta(cellchat, meta = meta)
    cellchat <- setIdent(cellchat, ident.use = "celltype")
    
    # Set up database
    cellchat@DB <- CellChatDB.human
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    # Project data with error handling
    tryCatch({
      cellchat <- projectData(cellchat, PPI.human)
    }, error = function(e) {
      cat("  ⚠️ Warning: PPI projection failed, continuing without it\n")
    })
    
    cat("  Significant interactions found:", nrow(cellchat@LR$LRsig), "\n")
    
    # Run inference
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    # Compute centrality
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    
    # Try pseudotime analysis if data available
    pseudotime_cellchats <- list()
    if (!is.null(pseudotime_data)) {
      cat("  Integrating pseudotime information...\n")
      
      # Get common cells
      common_cells <- intersect(colnames(data.input), names(pseudotime_data))
      
      if (length(common_cells) > 100) {
        # Create pseudotime bins
        pseudotime_values <- pseudotime_data[common_cells]
        pseudotime_bins <- cut(pseudotime_values, breaks = 3, 
                               labels = c("Early", "Middle", "Late"))
        
        for (stage in c("Early", "Middle", "Late")) {
          stage_cells <- common_cells[which(pseudotime_bins == stage & !is.na(pseudotime_bins))]
          
          if (length(stage_cells) >= 50) {
            cat("    Processing pseudotime stage:", stage, "with", length(stage_cells), "cells\n")
            
            stage_data <- data.input[, stage_cells]
            stage_meta <- meta[stage_cells, ]
            stage_meta$celltype <- droplevels(stage_meta$celltype)
            
            if (nlevels(stage_meta$celltype) >= 2) {
              tryCatch({
                stage_cellchat <- createCellChat(object = stage_data, meta = stage_meta, 
                                                 group.by = "celltype")
                stage_cellchat <- addMeta(stage_cellchat, meta = stage_meta)
                stage_cellchat <- setIdent(stage_cellchat, ident.use = "celltype")
                stage_cellchat@DB <- CellChatDB.human
                
                stage_cellchat <- subsetData(stage_cellchat)
                stage_cellchat <- identifyOverExpressedGenes(stage_cellchat)
                stage_cellchat <- identifyOverExpressedInteractions(stage_cellchat)
                stage_cellchat <- computeCommunProb(stage_cellchat)
                stage_cellchat <- filterCommunication(stage_cellchat, min.cells = 5)
                stage_cellchat <- computeCommunProbPathway(stage_cellchat)
                stage_cellchat <- aggregateNet(stage_cellchat)
                
                pseudotime_cellchats[[stage]] <- stage_cellchat
                cat("      ✅ Stage", stage, "completed\n")
              }, error = function(e) {
                cat("      ❌ Error in stage", stage, ":", e$message, "\n")
              })
            }
          }
        }
      }
    }
    
    # Store pseudotime analysis if available
    if (length(pseudotime_cellchats) > 0) {
      cellchat <- store_pseudotime_analysis(cellchat, pseudotime_cellchats)
      cat("  ✅ Pseudotime integration completed with", length(pseudotime_cellchats), "stages\n")
    }
    
    cat("  ✅ Advanced analysis completed successfully\n")
    return(cellchat)
    
  }, error = function(e) {
    cat("❌ Error in advanced analysis for", condition, ":", e$message, "\n")
    cat("🔄 Falling back to standard analysis...\n")
    
    # Fallback to basic analysis
    tryCatch({
      cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
      cellchat <- addMeta(cellchat, meta = meta)
      cellchat <- setIdent(cellchat, ident.use = "celltype")
      cellchat@DB <- CellChatDB.human
      cellchat <- subsetData(cellchat)
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat)
      cellchat <- computeCommunProb(cellchat)
      cellchat <- filterCommunication(cellchat, min.cells = 10)
      cellchat <- computeCommunProbPathway(cellchat)
      cellchat <- aggregateNet(cellchat)
      
      # Try to compute centrality for fallback too
      tryCatch({
        cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
      }, error = function(e2) {
        cat("  ⚠️ Could not compute centrality in fallback\n")
      })
      
      return(cellchat)
    }, error = function(e2) {
      cat("❌ Fallback analysis also failed:", e2$message, "\n")
      return(NULL)
    })
  })
}

# ================================
# PREPARE PSEUDOTIME DATA AND RUN ANALYSIS
# ================================

# Prepare pseudotime data for integration
epithelial_pseudotime <- NULL
if (exists("epithelial_subset")) {
  epithelial_pseudotime <- epithelial_subset@meta.data$pseudotime_lineage1
  names(epithelial_pseudotime) <- colnames(epithelial_subset)
  cat("✅ Prepared pseudotime data for", length(epithelial_pseudotime), "cells\n")
}

# Run analysis for all conditions
cellchat_ais_advanced <- run_advanced_analysis_safely("AIS", epithelial_pseudotime)
cellchat_iac_advanced <- run_advanced_analysis_safely("IAC", epithelial_pseudotime)  
cellchat_mia_advanced <- run_advanced_analysis_safely("MIA", epithelial_pseudotime)

# ================================
# VISUALIZATION FUNCTIONS
# ================================

# Function 1: Advanced Network Analysis
create_advanced_network_plots <- function(cellchat_obj, condition_name) {
  
  cat("🎨 Creating advanced network plots for", condition_name, "\n")
  
  if (is.null(cellchat_obj)) {
    cat("  ❌ No CellChat object provided for", condition_name, "\n")
    return(FALSE)
  }
  
  plot_dir <- file.path(base_dir, "results/publication_figures")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  # 1. Network centrality analysis
  tryCatch({
    pdf(file.path(plot_dir, paste0("Advanced_", condition_name, "_network_centrality.pdf")), 
        width = 16, height = 12)
    
    netAnalysis_signalingRole_scatter(cellchat_obj)
    
    if (!is.null(cellchat_obj@netP$centr)) {
      netAnalysis_signalingRole_heatmap(cellchat_obj, pattern = "outgoing")
      netAnalysis_signalingRole_heatmap(cellchat_obj, pattern = "incoming")
    }
    
    dev.off()
    cat("  ✅ Network centrality plots saved\n")
  }, error = function(e) {
    cat("  ❌ Error creating centrality plots:", e$message, "\n")
    if (dev.cur() != 1) dev.off()
  })
  
  # 2. Circle plots
  tryCatch({
    pdf(file.path(plot_dir, paste0("Advanced_", condition_name, "_circle_plots.pdf")), 
        width = 16, height = 8)
    
    par(mfrow = c(1, 2))
    groupSize <- as.numeric(table(cellchat_obj@idents))
    
    netVisual_circle(cellchat_obj@net$count, vertex.weight = groupSize, 
                     weight.scale = TRUE, label.edge = FALSE, vertex.label.cex = 1)
    
    netVisual_circle(cellchat_obj@net$weight, vertex.weight = groupSize, 
                     weight.scale = TRUE, label.edge = FALSE, vertex.label.cex = 1)
    
    dev.off()
    cat("  ✅ Circle plots saved\n")
  }, error = function(e) {
    cat("  ❌ Error creating circle plots:", e$message, "\n")
    if (dev.cur() != 1) dev.off()
  })
  
  # 3. Pathway analysis
  pathways <- cellchat_obj@netP$pathways
  if (length(pathways) > 0) {
    tryCatch({
      pdf(file.path(plot_dir, paste0("Advanced_", condition_name, "_pathways.pdf")), 
          width = 14, height = 10)
      
      top_pathways <- head(pathways, min(6, length(pathways)))
      
      par(mfrow = c(2, 3))
      for (pathway in top_pathways) {
        netVisual_aggregate(cellchat_obj, signaling = pathway, layout = "circle")
        title(paste(pathway, "in", condition_name))
      }
      
      dev.off()
      cat("  ✅ Pathway plots saved for", length(top_pathways), "pathways\n")
    }, error = function(e) {
      cat("  ❌ Error creating pathway plots:", e$message, "\n")
      if (dev.cur() != 1) dev.off()
    })
  }
  
  return(TRUE)
}

# Function 2: Comparative Analysis
perform_sophisticated_comparison <- function(cellchat_list) {
  
  cat("🔬 Performing sophisticated comparative analysis\n")
  
  valid_objects <- Filter(function(x) !is.null(x) && inherits(x, "CellChat"), cellchat_list)
  
  if (length(valid_objects) < 2) {
    cat("  ❌ Need at least 2 valid CellChat objects for comparison\n")
    return(NULL)
  }
  
  cat("  Found", length(valid_objects), "valid objects:", names(valid_objects), "\n")
  
  tryCatch({
    cellchat_merged <- mergeCellChat(valid_objects, add.names = names(valid_objects))
    
    pdf(file.path(base_dir, "results/publication_figures/Advanced_comparative_analysis.pdf"), 
        width = 20, height = 16)
    
    gg1 <- compareInteractions(cellchat_merged, show.legend = TRUE, 
                               group = seq_len(length(valid_objects)))
    gg2 <- compareInteractions(cellchat_merged, show.legend = TRUE, 
                               group = seq_len(length(valid_objects)), 
                               measure = "weight")
    
    print(gg1 + gg2 + plot_layout(ncol = 2))
    
    if (length(valid_objects) == 2) {
      par(mfrow = c(1, 2))
      netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "count")
      netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "weight")
    }
    
    tryCatch({
      gg3 <- rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
      gg4 <- rankNet(cellchat_merged, mode = "comparison", stacked = FALSE, do.stat = TRUE)
      print(gg3)
      print(gg4)
    }, error = function(e) {
      cat("  ⚠️ Warning: Could not create pathway ranking plots:", e$message, "\n")
    })
    
    dev.off()
    cat("  ✅ Comparative analysis completed\n")
    
    return(cellchat_merged)
    
  }, error = function(e) {
    cat("  ❌ Error in comparative analysis:", e$message, "\n")
    if (dev.cur() != 1) dev.off()
    return(NULL)
  })
}

# Function 3: Pseudotime Communication Analysis
analyze_pseudotime_communication <- function(cellchat_obj, condition_name) {
  
  pseudotime_analysis <- get_pseudotime_analysis(cellchat_obj)
  
  if (is.null(pseudotime_analysis)) {
    cat("⚠️  No pseudotime analysis available for", condition_name, "\n")
    return(NULL)
  }
  
  cat("🔄 Analyzing communication changes along pseudotime for", condition_name, "\n")
  
  if (length(pseudotime_analysis) < 2) {
    cat("⚠️  Not enough pseudotime stages for comparison (", length(pseudotime_analysis), ")\n")
    return(NULL)
  }
  
  tryCatch({
    pt_merged <- mergeCellChat(pseudotime_analysis, add.names = names(pseudotime_analysis))
    
    pdf(file.path(base_dir, "results/publication_figures", 
                  paste0("Advanced_", condition_name, "_pseudotime_communication.pdf")), 
        width = 18, height = 12)
    
    gg1 <- compareInteractions(pt_merged, show.legend = TRUE)
    print(gg1)
    
    if (length(pseudotime_analysis) == 2) {
      par(mfrow = c(1, 2))
      netVisual_diffInteraction(pt_merged, weight.scale = TRUE, measure = "count")
      netVisual_diffInteraction(pt_merged, weight.scale = TRUE, measure = "weight")
    }
    
    tryCatch({
      gg2 <- rankNet(pt_merged, mode = "comparison", stacked = FALSE, do.stat = TRUE)
      print(gg2)
    }, error = function(e) {
      cat("    ⚠️ Warning: Could not create pseudotime pathway ranking:", e$message, "\n")
    })
    
    dev.off()
    
    cat("  ✅ Pseudotime communication analysis completed\n")
    return(pt_merged)
    
  }, error = function(e) {
    cat("  ❌ Error in pseudotime analysis:", e$message, "\n")
    if (dev.cur() != 1) dev.off()
    return(NULL)
  })
}

# ================================
# RUN SOPHISTICATED ANALYSES
# ================================

# Store advanced objects
cellchat_advanced_list <- list(
  AIS = cellchat_ais_advanced,
  IAC = cellchat_iac_advanced,
  MIA = cellchat_mia_advanced
)

# Save advanced objects
qsave(cellchat_ais_advanced, file.path(base_dir, "data/processed/cellchat_AIS_advanced.qs"))
qsave(cellchat_iac_advanced, file.path(base_dir, "data/processed/cellchat_IAC_advanced.qs"))
qsave(cellchat_mia_advanced, file.path(base_dir, "data/processed/cellchat_MIA_advanced.qs"))

# Create advanced network plots
for (cond in names(cellchat_advanced_list)) {
  if (!is.null(cellchat_advanced_list[[cond]])) {
    create_advanced_network_plots(cellchat_advanced_list[[cond]], cond)
  }
}

# Sophisticated comparative analysis
cellchat_sophisticated_merged <- perform_sophisticated_comparison(cellchat_advanced_list)

# Pseudotime-communication integration analysis
for (cond in names(cellchat_advanced_list)) {
  if (!is.null(cellchat_advanced_list[[cond]])) {
    analyze_pseudotime_communication(cellchat_advanced_list[[cond]], cond)
  }
}

# ================================
# NETWORK METRICS COMPUTATION
# ================================

cat("\n📊 COMPUTING ADVANCED NETWORK METRICS\n")

compute_network_metrics_safe <- function(cellchat_obj, condition_name) {
  
  if (is.null(cellchat_obj)) {
    cat("⚠️  No CellChat object for", condition_name, "\n")
    return(NULL)
  }
  
  tryCatch({
    net_matrix <- cellchat_obj@net$weight
    
    if (is.null(net_matrix) || all(net_matrix == 0) || any(is.na(net_matrix)) || any(is.infinite(net_matrix))) {
      cat("⚠️  Invalid network matrix for", condition_name, "\n")
      return(NULL)
    }
    
    g <- graph_from_adjacency_matrix(net_matrix, mode = "directed", weighted = TRUE, diag = FALSE)
    g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
    
    if (vcount(g) == 0 || ecount(g) == 0) {
      cat("⚠️  Empty graph for", condition_name, "\n")
      return(NULL)
    }
    
    metrics <- list(
      condition = condition_name,
      n_nodes = vcount(g),
      n_edges = ecount(g)
    )
    
    metrics$density <- tryCatch(edge_density(g), error = function(e) NA)
    metrics$reciprocity <- tryCatch(reciprocity(g), error = function(e) NA)
    metrics$transitivity <- tryCatch(transitivity(g, type = "global"), error = function(e) NA)
    
    if (is.connected(g)) {
      metrics$diameter <- tryCatch(diameter(g, weights = NA), error = function(e) NA)
      metrics$avg_path_length <- tryCatch(mean_distance(g, weights = NA), error = function(e) NA)
    } else {
      metrics$diameter <- NA
      metrics$avg_path_length <- NA
    }
    
    metrics$centralization_degree <- tryCatch(centr_degree(g)$centralization, error = function(e) NA)
    metrics$centralization_betweenness <- tryCatch(centr_betw(g)$centralization, error = function(e) NA)
    
    metrics$modularity <- tryCatch({
      if (ecount(g) > 0) {
        wc <- cluster_walktrap(g)
        modularity(wc)
      } else {
        NA
      }
    }, error = function(e) NA)
    
    cat("  ✅ Metrics computed for", condition_name, "\n")
    return(metrics)
    
  }, error = function(e) {
    cat("❌ Error computing metrics for", condition_name, ":", e$message, "\n")
    return(NULL)
  })
}

# Compute metrics for all conditions
network_metrics <- list()
for (cond in names(cellchat_advanced_list)) {
  if (!is.null(cellchat_advanced_list[[cond]])) {
    metrics <- compute_network_metrics_safe(cellchat_advanced_list[[cond]], cond)
    if (!is.null(metrics)) {
      network_metrics[[cond]] <- metrics
    }
  }
}

# Convert to dataframe and save
if (length(network_metrics) > 0) {
  metrics_df <- do.call(rbind, lapply(network_metrics, function(x) {
    data.frame(x[sapply(x, function(y) length(y) == 1 && !is.null(y))])
  }))
  
  write.csv(metrics_df, file.path(base_dir, "results/tables/advanced_network_metrics.csv"), 
            row.names = FALSE)
  
  if (nrow(metrics_df) > 0 && ncol(metrics_df) > 2) {
    metrics_long <- reshape2::melt(metrics_df, id.vars = "condition")
    metrics_long <- metrics_long[!is.na(metrics_long$value), ]
    
    if (nrow(metrics_long) > 0) {
      p_metrics <- ggplot(metrics_long, aes(x = condition, y = value, fill = condition)) +
        geom_bar(stat = "identity", alpha = 0.8) +
        facet_wrap(~variable, scales = "free_y") +
        scale_fill_viridis_d() +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("Advanced Network Metrics Across Disease Stages")
      
      ggsave(file.path(base_dir, "results/publication_figures/Advanced_network_metrics.pdf"), 
             p_metrics, width = 16, height = 12, dpi = 300)
      
      cat("✅ Network metrics analysis completed\n")
    }
  }
}

# ================================
# SAVE COMPREHENSIVE RESULTS
# ================================

# Save all results
final_results <- list(
  epithelial_subset = epithelial_subset,
  sce = sce,
  lineage1_genes = sig_genes_lin1,
  lineage2_genes = sig_genes_lin2,
  differential_genes = lineage_diff_genes,
  specific_genes_analysis = gene_analysis_results,
  cellchat_advanced = cellchat_advanced_list,
  cellchat_sophisticated_merged = cellchat_sophisticated_merged,
  network_metrics = if(exists("metrics_df")) metrics_df else NULL
)

qsave(final_results, file.path(base_dir, "data/processed/complete_analysis_integrated.qs"))

# Save gene tables
if(exists("sig_genes_lin1")) {
  write.csv(sig_genes_lin1, file.path(base_dir, "results/tables/lineage1_genes.csv"), row.names = FALSE)
}
if(exists("sig_genes_lin2")) {
  write.csv(sig_genes_lin2, file.path(base_dir, "results/tables/lineage2_genes.csv"), row.names = FALSE)
}
if(exists("lineage_diff_genes")) {
  write.csv(lineage_diff_genes, file.path(base_dir, "results/tables/differential_genes.csv"), row.names = FALSE)
}

# Gene comparison table
if(exists("gene_analysis_results") && length(gene_analysis_results) > 1) {
  comparison_summary <- data.frame(
    Gene = names(gene_analysis_results),
    Pct_Expressing = sapply(gene_analysis_results, function(x) x$pct_expressing),
    Mean_Expression = sapply(gene_analysis_results, function(x) mean(x$expression_data$expression)),
    In_Lineage1 = sapply(gene_analysis_results, function(x) x$trajectory_status$in_lineage1),
    In_Lineage2 = sapply(gene_analysis_results, function(x) x$trajectory_status$in_lineage2),
    Is_Differential = sapply(gene_analysis_results, function(x) x$trajectory_status$is_differential)
  )
  
  write.csv(comparison_summary, file.path(base_dir, "results/tables/specific_genes_comparison.csv"), 
            row.names = FALSE)
}

# ================================
# FINAL SUMMARY
# ================================

cat("\n🎉 INTEGRATED ANALYSIS COMPLETED!\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("🔬 Features Applied:\n")
cat("  ✅ Deep pseudotime analysis with slingshot and tradeSeq\n")
cat("  ✅ Advanced CellChat analysis with pseudotime integration\n") 
cat("  ✅ Comprehensive gene-specific analysis (including KIF20A)\n")
cat("  ✅ Network topology metrics with igraph\n")
cat("  ✅ Publication-quality visualizations\n")
cat("  ✅ Sophisticated comparative analysis\n")

cat("\n📁 Results Generated:\n")
cat("  - Pseudotime trajectory plots\n")
cat("  - Gene-specific comprehensive analysis\n")
cat("  - Network centrality analysis\n")
cat("  - Comparative pathway analysis\n")
cat("  - Advanced network metrics\n")

cat("\n📊 Analysis Summary:\n")
cat("  - Lineage 1 genes (Direct malignant):", nrow(sig_genes_lin1), "\n")
cat("  - Lineage 2 genes (Transdifferentiation):", nrow(sig_genes_lin2), "\n")
cat("  - Differential genes between lineages:", nrow(lineage_diff_genes), "\n")
cat("  - Genes analyzed:", length(gene_analysis_results), "\n")

for (cond in names(cellchat_advanced_list)) {
  if (!is.null(cellchat_advanced_list[[cond]])) {
    n_pathways <- length(cellchat_advanced_list[[cond]]@netP$pathways)
    n_interactions <- nrow(cellchat_advanced_list[[cond]]@net$count)
    pseudotime_stages <- length(get_pseudotime_analysis(cellchat_advanced_list[[cond]]))
    
    cat("  ", cond, ":")
    cat(" Pathways:", n_pathways)
    cat(" | Cell types:", n_interactions)
    cat(" | Pseudotime stages:", pseudotime_stages, "\n")
  }
}

if (exists("metrics_df")) {
  cat("  Network metrics computed for", nrow(metrics_df), "conditions\n")
}

cat("\n💾 All results saved to:", base_dir, "\n")
