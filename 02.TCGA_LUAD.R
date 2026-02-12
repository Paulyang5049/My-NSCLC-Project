#!/usr/bin/env R

# =============================================================================
# CORRECTED TCGA-LUAD Analysis with Enhanced Visualizations
# Fixed scoping issues and improved error handling
# =============================================================================

# Clear environment and load libraries
rm(list = ls())

# Load required libraries
required_packages <- c(
  "TCGAbiolinks", "DESeq2", "org.Hs.eg.db", "AnnotationDbi", 
  "dplyr", "ggplot2", "pheatmap", "scatterplot3d", "factoextra",
  "FactoMineR", "openxlsx", "qs", "RColorBrewer", "limma",
  "edgeR", "biomaRt"
)

# Install and load packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Set working directory and create output folders
setwd("/Volumes/Haoran's SSD/Scientific research/My research/LUAD_我的第一篇文章/TCGA_LUAD_Analysis")
dir.create("results", showWarnings = FALSE)
dir.create("results/plots", showWarnings = FALSE)
dir.create("results/data", showWarnings = FALSE)

# =============================================================================
# 1. DATA DOWNLOAD AND PREPARATION (CORRECTED)
# =============================================================================

cat("=== Starting TCGA-LUAD Data Download and Processing ===\n")

# Download TCGA-LUAD data
download_tcga_data <- function() {
  cat("Downloading TCGA-LUAD data...\n")
  
  query <- GDCquery(
    project = "TCGA-LUAD",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor", "Solid Tissue Normal")
  )
  
  GDCdownload(query)
  data <- GDCprepare(query)
  
  return(data)
}

# Prepare expression matrix and sample information
prepare_data <- function(data) {
  cat("Preparing expression matrix and sample information...\n")
  
  # Extract count matrix and sample info
  count_matrix <- assay(data)
  sample_info <- colData(data)
  
  # Clean sample information
  sample_info$patient_id <- substr(sample_info$submitter_id, 1, 12)
  sample_info$sample_type <- factor(sample_info$sample_type)
  
  # Remove samples with missing sample type
  valid_samples <- !is.na(sample_info$sample_type)
  count_matrix <- count_matrix[, valid_samples]
  sample_info <- sample_info[valid_samples, ]
  
  # Create Group variable matching your style
  Group_tcga <- ifelse(sample_info$sample_type == "Primary Tumor", "Tumor", "Control")
  Group_tcga <- factor(Group_tcga, levels = c("Control", "Tumor"))
  names(Group_tcga) <- colnames(count_matrix)  # Important: add names!
  
  # Gene ID conversion
  cat("Converting gene IDs...\n")
  rownames(count_matrix) <- gsub("\\..*", "", rownames(count_matrix))
  
  gene_symbols <- mapIds(org.Hs.eg.db,
                         keys = rownames(count_matrix),
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
  
  # Handle missing symbols
  gene_symbols[is.na(gene_symbols)] <- rownames(count_matrix)[is.na(gene_symbols)]
  
  # Handle duplicated gene symbols
  duplicated_symbols <- duplicated(gene_symbols) | duplicated(gene_symbols, fromLast = TRUE)
  gene_symbols[duplicated_symbols] <- paste0(gene_symbols[duplicated_symbols], 
                                             "_", rownames(count_matrix)[duplicated_symbols])
  
  rownames(count_matrix) <- gene_symbols
  
  # Filter low count genes
  min_samples <- ceiling(0.1 * ncol(count_matrix))
  keep_genes <- rowSums(count_matrix >= 10) >= min_samples
  mrna_expr_count <- count_matrix[keep_genes, ]
  
  cat("Data preparation complete:\n")
  cat("Samples:", ncol(mrna_expr_count), "\n")
  cat("Genes:", nrow(mrna_expr_count), "\n")
  cat("Tumor samples:", sum(Group_tcga == "Tumor"), "\n")
  cat("Control samples:", sum(Group_tcga == "Control"), "\n")
  
  return(list(
    mrna_expr_count = mrna_expr_count,
    Group_tcga = Group_tcga,
    sample_info = sample_info,
    clinical_indexed = sample_info
  ))
}

# =============================================================================
# 2. DIFFERENTIAL EXPRESSION ANALYSIS (FIXED)
# =============================================================================

# Fixed differential analysis function
diff_analysis <- function(exprset, Group_tcga, project = "TCGA-LUAD", save = FALSE) {
  cat("Performing differential expression analysis...\n")
  
  # Ensure Group_tcga matches column names
  if (!all(names(Group_tcga) == colnames(exprset))) {
    if (length(Group_tcga) == ncol(exprset)) {
      names(Group_tcga) <- colnames(exprset)
    } else {
      stop("Group_tcga length doesn't match expression matrix columns")
    }
  }
  
  # Create sample data frame
  sample_df <- data.frame(
    condition = Group_tcga,
    row.names = colnames(exprset)
  )
  
  # DESeq2 analysis
  dds <- DESeqDataSetFromMatrix(
    countData = exprset,
    colData = sample_df,
    design = ~ condition
  )
  
  dds$condition <- relevel(dds$condition, ref = "Control")
  dds <- DESeq(dds)
  
  # Get results
  deseq2_results <- results(dds, contrast = c("condition", "Tumor", "Control"))
  
  # Convert to limma-style results for compatibility
  limma_result <- data.frame(
    logFC = deseq2_results$log2FoldChange,
    P.Value = deseq2_results$pvalue,
    adj.P.Val = deseq2_results$padj,
    AveExpr = deseq2_results$baseMean,
    gene_symbol = rownames(deseq2_results),
    stringsAsFactors = FALSE
  )
  
  # Remove rows with NA values
  limma_result <- limma_result[complete.cases(limma_result), ]
  
  return(list(
    dds = dds,
    limma_result = limma_result,
    deseq2_results = deseq2_results
  ))
}

# Add change column function
add_change_column <- function(deg, logFC_t = 1, P.Value_t = 0.05) {
  k1 <- (deg$P.Value < P.Value_t) & (deg$logFC < -logFC_t)
  k2 <- (deg$P.Value < P.Value_t) & (deg$logFC > logFC_t)
  deg <- mutate(deg, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
  return(deg)
}

# =============================================================================
# 3. VISUALIZATION FUNCTIONS (Your Preferred Style)
# =============================================================================

# Enhanced PCA Analysis with 3D plot
perform_pca_analysis <- function(mrna_expr_count, Group_tcga) {
  cat("Performing PCA analysis...\n")
  
  # Define colors matching your style
  color.bin <- c("#00599F", "#D01910")
  
  # Perform PCA (transpose for samples as rows)
  res.pca.comp <- prcomp(t(mrna_expr_count), scale = FALSE)
  
  # Prepare plot data
  plot.data <- as.data.frame(res.pca.comp$x[, 1:min(10, ncol(res.pca.comp$x))])
  plot.data <- plot.data %>% 
    mutate(ID = rownames(plot.data),
           Type = Group_tcga[rownames(plot.data)],
           TypeColor = color.bin[as.numeric(as.factor(Group_tcga[rownames(plot.data)]))])
  
  # 3D PCA Plot
  pdf("results/plots/PCA3d.pdf", width = 7, height = 7)
  scatterplot3d(x = plot.data$PC2, 
                y = plot.data$PC1, 
                z = plot.data$PC3,
                color = plot.data$TypeColor,
                pch = 16, cex.symbols = 1,
                scale.y = 0.7, angle = 45,
                xlab = "PC2", ylab = "PC1", zlab = "PC3",
                main = "3D Scatter Plot of PCA",
                col.axis = "#444444", col.grid = "#CCCCCC")
  legend("bottom", legend = levels(as.factor(Group_tcga)),
         col = color.bin, pch = 16,
         inset = -0.3, xpd = TRUE, horiz = TRUE)
  dev.off()
  
  # Save PCA data
  write.xlsx(plot.data, "results/data/PCA3d.xlsx", overwrite = TRUE)
  
  # Variance explained plot
  p <- fviz_eig(res.pca.comp)
  var_explained <- get_eig(res.pca.comp)
  
  ggsave("results/plots/PCA_percent.pdf", p, width = 5, height = 5)
  write.xlsx(var_explained, "results/data/PCA_percent.xlsx", rowNames = TRUE, overwrite = TRUE)
  
  return(list(
    pca_result = res.pca.comp,
    plot_data = plot.data,
    var_explained = var_explained
  ))
}

# Enhanced Heatmap Function
generate_heatmap <- function(exp, Group) {
  cat("Generating heatmap...\n")
  
  # Define custom color palette matching your style
  custom_palette <- colorRampPalette(c("#599CB4", "#92B5CA", "#AECFD4", "#CCE4EF",
                                       "#F5DFDB", "#EDB8B0", "#E69191", "#C25759"))(100)
  
  # Calculate standard deviation for each row and select top 1000 genes
  cg <- names(tail(sort(apply(exp, 1, sd)), 1000))
  n <- exp[cg, ]
  
  # Order columns by group to cluster Control and Tumor separately
  ordered_indices <- order(Group)
  n <- n[, ordered_indices]
  Group_ordered <- Group[ordered_indices]
  
  # Create annotation data frame
  annotation_col <- data.frame(group = Group_ordered)
  rownames(annotation_col) <- colnames(n)
  
  # Generate the heatmap
  pdf("results/plots/heatmap.pdf", width = 12, height = 10)
  heatmap <- pheatmap(n,
                      color = custom_palette,
                      show_colnames = FALSE,
                      show_rownames = FALSE,
                      annotation_col = annotation_col,
                      scale = "row",
                      breaks = seq(-3, 3, length.out = 100),
                      fontsize = 10,
                      fontsize_row = 6,
                      fontsize_col = 6,
                      main = "TCGA-LUAD Heatmap",
                      clustering_distance_rows = "euclidean",
                      clustering_method = "complete",
                      cluster_cols = FALSE)
  dev.off()
  
  return(heatmap)
}

# Enhanced Volcano Plot Function
draw_volcano <- function(data, logFC_cutoff = 1, p_value_cutoff = 0.05, title = "Volcano Plot") {
  cat("Drawing volcano plot...\n")
  
  data$change <- "stable"
  data$change[data$logFC > logFC_cutoff & data$P.Value < p_value_cutoff] <- "up"
  data$change[data$logFC < -logFC_cutoff & data$P.Value < p_value_cutoff] <- "down"
  
  volcano <- ggplot(data, aes(x = logFC, y = -log10(P.Value), color = change)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("down" = "#2B6688", "up" = "#F1A93B", "stable" = "#A8ACB9")) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "gray90")
    ) +
    labs(x = "log2 Fold Change", y = "-log10 P Value") +
    ggtitle(title) +
    geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(p_value_cutoff), linetype = "dashed", color = "black")
  
  ggsave("results/plots/volcano_plot.pdf", volcano, width = 8, height = 6)
  
  return(volcano)
}

# =============================================================================
# 4. MAIN ANALYSIS WORKFLOW (CORRECTED)
# =============================================================================

main_analysis <- function() {
  cat("=== Starting Main Analysis Workflow ===\n")
  
  # Step 1: Download and prepare data
  tcga_data <- download_tcga_data()
  data_list <- prepare_data(tcga_data)
  
  # Extract components
  mrna_expr_count <- data_list$mrna_expr_count
  Group_tcga <- data_list$Group_tcga
  sample_info <- data_list$sample_info
  clinical_indexed <- data_list$clinical_indexed
  
  # Step 2: Differential expression analysis (pass Group_tcga as parameter)
  luad_deg_result <- diff_analysis(exprset = mrna_expr_count, 
                                   Group_tcga = Group_tcga, 
                                   project = "TCGA-LUAD", 
                                   save = FALSE)
  limma_result <- luad_deg_result$limma_result
  
  # Add change column
  LUAD_deg <- add_change_column(limma_result, logFC_t = 1, P.Value_t = 0.05)
  
  # Step 3: Generate visualizations
  # PCA Analysis
  pca_results <- perform_pca_analysis(mrna_expr_count, Group_tcga)
  
  # Heatmap
  heatmap <- generate_heatmap(mrna_expr_count, Group_tcga)
  
  # Volcano Plot
  volcano_plot <- draw_volcano(LUAD_deg, logFC_cutoff = 1, p_value_cutoff = 0.05, title = "TCGA-LUAD Volcano Plot")
  
  # Step 4: Summary statistics
  cat("\n=== Analysis Summary ===\n")
  cat("Total genes analyzed:", nrow(LUAD_deg), "\n")
  cat("Significantly upregulated genes:", sum(LUAD_deg$change == "up"), "\n")
  cat("Significantly downregulated genes:", sum(LUAD_deg$change == "down"), "\n")
  cat("Stable genes:", sum(LUAD_deg$change == "stable"), "\n")
  
  # Step 5: Save all results using qs
  cat("\n=== Saving Results ===\n")
  
  # Save main results using qs (fast save/load)
  qsave(mrna_expr_count, "results/data/mrna_expr_count.qs")
  qsave(Group_tcga, "results/data/Group_tcga.qs")
  qsave(clinical_indexed, "results/data/clinical_indexed.qs")
  qsave(LUAD_deg, "results/data/LUAD_deg.qs")
  qsave(luad_deg_result, "results/data/luad_deg_result.qs")
  qsave(pca_results, "results/data/pca_results.qs")
  
  # Save complete workspace in your preferred format
  save(mrna_expr_count, Group_tcga, clinical_indexed, LUAD_deg, luad_deg_result, 
       pca_results, file = "results/data/LUAD.Rdata")
  
  # Save as single qs file (fastest option)
  all_results <- list(
    mrna_expr_count = mrna_expr_count,
    Group_tcga = Group_tcga,
    clinical_indexed = clinical_indexed,
    LUAD_deg = LUAD_deg,
    luad_deg_result = luad_deg_result,
    pca_results = pca_results,
    sample_info = sample_info
  )
  qsave(all_results, "results/data/TCGA_LUAD_complete_analysis.qs")
  
  # Save differential expression results as Excel
  write.xlsx(LUAD_deg, "results/data/LUAD_differential_expression.xlsx", overwrite = TRUE)
  
  # Save significant genes separately
  sig_up <- LUAD_deg[LUAD_deg$change == "up", ]
  sig_down <- LUAD_deg[LUAD_deg$change == "down", ]
  
  write.xlsx(list(
    "All_Results" = LUAD_deg,
    "Upregulated" = sig_up,
    "Downregulated" = sig_down
  ), "results/data/LUAD_significant_genes.xlsx", overwrite = TRUE)
  
  cat("Analysis complete! All results saved in 'results' directory.\n")
  cat("Main data files:\n")
  cat("- TCGA_LUAD_complete_analysis.qs (fastest loading)\n")
  cat("- LUAD.Rdata (R workspace format)\n")
  cat("- Individual .qs files for each component\n")
  
  return(all_results)
}

# =============================================================================
# 5. QUICK LOAD FUNCTION
# =============================================================================

# Function to quickly load all results
load_tcga_results <- function(method = "qs") {
  if (method == "qs") {
    # Fastest method
    results <- qread("results/data/TCGA_LUAD_complete_analysis.qs")
    list2env(results, envir = .GlobalEnv)
    cat("Data loaded successfully using qs format!\n")
  } else if (method == "rdata") {
    # Traditional method
    load("results/data/LUAD.Rdata", envir = .GlobalEnv)
    cat("Data loaded successfully using Rdata format!\n")
  }
}

# =============================================================================
# 6. RUN ANALYSIS
# =============================================================================

# Execute the complete analysis
cat("Starting TCGA-LUAD integrated analysis...\n")
cat("This will take some time for initial download.\n")
cat("Data seems to be already downloaded, so this should be faster.\n\n")

# Run the analysis
results <- main_analysis()

cat("\n=== Usage Examples ===\n")
cat("To reload data quickly: load_tcga_results('qs')\n")
cat("To reload with Rdata: load_tcga_results('rdata')\n")
cat("Individual components: mrna_expr_count <- qread('results/data/mrna_expr_count.qs')\n")

