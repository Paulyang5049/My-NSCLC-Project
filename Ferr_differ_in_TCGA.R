# Load required libraries
library(openxlsx)
library(dplyr)
library(qs)
library(ggplot2)
library(ggrepel)

# Set working directory
setwd("/Volumes/Haoran's SSD/Scientific research/My research/LUAD_我的第一篇文章/TCGA_LUAD_Analysis")

# Define gene groups
driver_genes <- c(
  "AOC1", "METTL1", "TGFB2", "SMAD7", "CRYAB", "EPAS1", "SOCS2", "KLF2", 
  "APOL3", "CDO1", "COX4I2", "CD36", "SNCA", "CLTRN", "COX7A1", "SCARA5", 
  "CAVIN1", "LINC00472", "SYNPR-AS1", "PVT1", "ARGLU1-DT", "H3C2"
)

suppressor_genes <- c(
  "PDK4", "LTF", "CYP24A1", "CP", "ENO1", "FBLN1", "TFAP2C", "RARRES2", 
  "SLC35F2", "TULP1", "AADAC", "NR5A2", "HELLS", "PRDX4", "HSPA12B", "KL", 
  "ETS1", "IL33", "TFAP2A", "BRDT", "FGF2", "GATA6", "RXRG", "HYOU1", 
  "PROM2", "PDIA4", "PROK2", "CCT3", "AIM2", "SPDYA", "KLF15", "CST1", 
  "AQP4", "ETV4", "B3GNT3", "PYCR1", "LRRK2", "C8orf76", "TCF4", "AGER", 
  "GPX3", "EMP2", "VAMP2", "ITGB2-AS1", "LINC01134", "SELENOP"
)

# Verify all genes are accounted for
all_genes <- c(driver_genes, suppressor_genes)
cat("Total driver genes:", length(driver_genes), "\n")
cat("Total suppressor genes:", length(suppressor_genes), "\n")
cat("Total genes:", length(all_genes), "\n\n")

# Function to analyze gene groups separately
analyze_gene_groups <- function(driver_genes, suppressor_genes) {
  
  cat("=== TCGA-LUAD Analysis: Driver vs Suppressor Genes ===\n\n")
  
  # Load differential expression results
  if (file.exists("results/data/LUAD_differential_expression.xlsx")) {
    cat("Loading differential expression results from Excel...\n")
    deg_results <- read.xlsx("results/data/LUAD_differential_expression.xlsx")
  } else if (file.exists("results/data/LUAD_deg.qs")) {
    cat("Loading differential expression results from qs file...\n")
    deg_results <- qread("results/data/LUAD_deg.qs")
  } else {
    cat("Differential expression results not found! Please run the analysis first.\n")
    return(NULL)
  }
  
  cat("Total genes in differential expression results:", nrow(deg_results), "\n\n")
  
  # Find driver genes
  driver_matches <- deg_results[deg_results$gene_symbol %in% driver_genes, ]
  driver_found <- driver_matches$gene_symbol
  driver_missing <- driver_genes[!driver_genes %in% driver_found]
  
  # Find suppressor genes
  suppressor_matches <- deg_results[deg_results$gene_symbol %in% suppressor_genes, ]
  suppressor_found <- suppressor_matches$gene_symbol
  suppressor_missing <- suppressor_genes[!suppressor_genes %in% suppressor_found]
  
  # Add gene type column
  driver_matches$gene_type <- "Driver"
  suppressor_matches$gene_type <- "Suppressor"
  
  # Summary statistics
  cat("=== SEARCH RESULTS SUMMARY ===\n")
  cat("Driver genes searched:", length(driver_genes), "\n")
  cat("Driver genes found:", nrow(driver_matches), "\n")
  cat("Driver genes missing:", length(driver_missing), "\n")
  cat("Suppressor genes searched:", length(suppressor_genes), "\n")
  cat("Suppressor genes found:", nrow(suppressor_matches), "\n")
  cat("Suppressor genes missing:", length(suppressor_missing), "\n\n")
  
  # Analysis of Driver genes
  if (nrow(driver_matches) > 0) {
    cat("=== DRIVER GENES ANALYSIS ===\n")
    
    # Sort by significance
    driver_sorted <- driver_matches[order(driver_matches$P.Value), ]
    print(driver_sorted[, c("gene_symbol", "logFC", "P.Value", "adj.P.Val", "change")])
    
    cat("\nDriver genes regulation summary:\n")
    driver_reg_summary <- table(driver_matches$change)
    print(driver_reg_summary)
    
    cat("\nDriver genes statistics:\n")
    cat("Upregulated (logFC > 1, p < 0.05):", 
        sum(driver_matches$logFC > 1 & driver_matches$P.Value < 0.05), "\n")
    cat("Downregulated (logFC < -1, p < 0.05):", 
        sum(driver_matches$logFC < -1 & driver_matches$P.Value < 0.05), "\n")
    cat("Significant (p < 0.05):", sum(driver_matches$P.Value < 0.05), "\n")
    cat("Mean log2FC:", round(mean(driver_matches$logFC, na.rm = TRUE), 3), "\n")
    cat("Median log2FC:", round(median(driver_matches$logFC, na.rm = TRUE), 3), "\n\n")
    
    if (length(driver_missing) > 0) {
      cat("Missing driver genes:", paste(driver_missing, collapse = ", "), "\n\n")
    }
  }
  
  # Analysis of Suppressor genes
  if (nrow(suppressor_matches) > 0) {
    cat("=== SUPPRESSOR GENES ANALYSIS ===\n")
    
    # Sort by significance
    suppressor_sorted <- suppressor_matches[order(suppressor_matches$P.Value), ]
    print(suppressor_sorted[, c("gene_symbol", "logFC", "P.Value", "adj.P.Val", "change")])
    
    cat("\nSuppressor genes regulation summary:\n")
    suppressor_reg_summary <- table(suppressor_matches$change)
    print(suppressor_reg_summary)
    
    cat("\nSuppressor genes statistics:\n")
    cat("Upregulated (logFC > 1, p < 0.05):", 
        sum(suppressor_matches$logFC > 1 & suppressor_matches$P.Value < 0.05), "\n")
    cat("Downregulated (logFC < -1, p < 0.05):", 
        sum(suppressor_matches$logFC < -1 & suppressor_matches$P.Value < 0.05), "\n")
    cat("Significant (p < 0.05):", sum(suppressor_matches$P.Value < 0.05), "\n")
    cat("Mean log2FC:", round(mean(suppressor_matches$logFC, na.rm = TRUE), 3), "\n")
    cat("Median log2FC:", round(median(suppressor_matches$logFC, na.rm = TRUE), 3), "\n\n")
    
    if (length(suppressor_missing) > 0) {
      cat("Missing suppressor genes:", paste(suppressor_missing, collapse = ", "), "\n\n")
    }
  }
  
  # Comparative analysis
  cat("=== COMPARATIVE ANALYSIS: DRIVERS vs SUPPRESSORS ===\n")
  
  if (nrow(driver_matches) > 0 && nrow(suppressor_matches) > 0) {
    # Compare mean expression changes
    driver_mean_fc <- mean(driver_matches$logFC, na.rm = TRUE)
    suppressor_mean_fc <- mean(suppressor_matches$logFC, na.rm = TRUE)
    
    cat("Mean log2FC - Drivers:", round(driver_mean_fc, 3), "\n")
    cat("Mean log2FC - Suppressors:", round(suppressor_mean_fc, 3), "\n")
    cat("Difference (Driver - Suppressor):", round(driver_mean_fc - suppressor_mean_fc, 3), "\n\n")
    
    # Statistical test
    fc_test <- t.test(driver_matches$logFC, suppressor_matches$logFC)
    cat("T-test comparing log2FC between groups:\n")
    cat("P-value:", format(fc_test$p.value, scientific = TRUE, digits = 3), "\n")
    cat("Interpretation:", ifelse(fc_test$p.value < 0.05, "Significant difference", "No significant difference"), "\n\n")
    
    # Significance rates
    driver_sig_rate <- sum(driver_matches$P.Value < 0.05) / nrow(driver_matches) * 100
    suppressor_sig_rate <- sum(suppressor_matches$P.Value < 0.05) / nrow(suppressor_matches) * 100
    
    cat("Significance rates (p < 0.05):\n")
    cat("Drivers:", round(driver_sig_rate, 1), "%\n")
    cat("Suppressors:", round(suppressor_sig_rate, 1), "%\n\n")
  }
  
  # Save results
  cat("=== SAVING RESULTS ===\n")
  
  # Combine all results
  all_found_genes <- rbind(driver_matches, suppressor_matches)
  
  # Create comprehensive Excel file
  gene_analysis_results <- list(
    "Summary" = data.frame(
      Gene_Type = c("Driver", "Suppressor", "Total"),
      Searched = c(length(driver_genes), length(suppressor_genes), length(c(driver_genes, suppressor_genes))),
      Found = c(nrow(driver_matches), nrow(suppressor_matches), nrow(all_found_genes)),
      Missing = c(length(driver_missing), length(suppressor_missing), length(c(driver_missing, suppressor_missing))),
      Upregulated = c(sum(driver_matches$change == "up"), sum(suppressor_matches$change == "up"), 
                      sum(all_found_genes$change == "up")),
      Downregulated = c(sum(driver_matches$change == "down"), sum(suppressor_matches$change == "down"),
                        sum(all_found_genes$change == "down")),
      Significant = c(sum(driver_matches$P.Value < 0.05), sum(suppressor_matches$P.Value < 0.05),
                      sum(all_found_genes$P.Value < 0.05))
    ),
    "All_Found_Genes" = all_found_genes,
    "Driver_Genes" = driver_matches,
    "Suppressor_Genes" = suppressor_matches
  )
  
  # Add missing genes if any
  if (length(driver_missing) > 0 || length(suppressor_missing) > 0) {
    missing_df <- data.frame(
      Gene_Symbol = c(driver_missing, suppressor_missing),
      Gene_Type = c(rep("Driver", length(driver_missing)), rep("Suppressor", length(suppressor_missing)))
    )
    gene_analysis_results[["Missing_Genes"]] <- missing_df
  }
  
  # Add statistical comparison
  if (nrow(driver_matches) > 0 && nrow(suppressor_matches) > 0) {
    gene_analysis_results[["Statistical_Comparison"]] <- data.frame(
      Metric = c("Mean_log2FC_Drivers", "Mean_log2FC_Suppressors", "T_test_p_value", 
                 "Sig_rate_Drivers_percent", "Sig_rate_Suppressors_percent"),
      Value = c(round(driver_mean_fc, 3), round(suppressor_mean_fc, 3), 
                format(fc_test$p.value, scientific = TRUE, digits = 3),
                round(driver_sig_rate, 1), round(suppressor_sig_rate, 1))
    )
  }
  
  write.xlsx(gene_analysis_results, "results/data/Driver_vs_Suppressor_Analysis.xlsx", overwrite = TRUE)
  
  # Save individual files
  qsave(all_found_genes, "results/data/all_driver_suppressor_genes.qs")
  qsave(driver_matches, "results/data/driver_genes_results.qs")
  qsave(suppressor_matches, "results/data/suppressor_genes_results.qs")
  
  cat("Results saved to:\n")
  cat("- results/data/Driver_vs_Suppressor_Analysis.xlsx\n")
  cat("- results/data/all_driver_suppressor_genes.qs\n")
  cat("- results/data/driver_genes_results.qs\n")
  cat("- results/data/suppressor_genes_results.qs\n\n")
  
  return(list(
    driver_genes = driver_matches,
    suppressor_genes = suppressor_matches,
    all_genes = all_found_genes,
    driver_missing = driver_missing,
    suppressor_missing = suppressor_missing,
    statistics = if(exists("fc_test")) list(
      driver_mean_fc = driver_mean_fc,
      suppressor_mean_fc = suppressor_mean_fc,
      t_test = fc_test,
      driver_sig_rate = driver_sig_rate,
      suppressor_sig_rate = suppressor_sig_rate
    ) else NULL
  ))
}

# Enhanced visualization for grouped analysis
create_grouped_visualizations <- function(analysis_results) {
  
  if (is.null(analysis_results)) {
    cat("No results to visualize.\n")
    return(NULL)
  }
  
  driver_genes <- analysis_results$driver_genes
  suppressor_genes <- analysis_results$suppressor_genes
  all_genes <- analysis_results$all_genes
  
  if (nrow(all_genes) == 0) {
    cat("No genes found to visualize.\n")
    return(NULL)
  }
  
  # 1. Combined volcano plot with gene types
  volcano_plot <- ggplot(all_genes, aes(x = logFC, y = -log10(P.Value), 
                                        color = gene_type, shape = change)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text_repel(aes(label = gene_symbol), size = 3, 
                    box.padding = 0.5, max.overlaps = Inf,
                    segment.color = "gray50") +
    scale_color_manual(values = c("Driver" = "#E74C3C", "Suppressor" = "#3498DB"),
                       name = "Gene Type") +
    scale_shape_manual(values = c("up" = 17, "down" = 25, "stable" = 16),
                       name = "Regulation") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.position = "bottom",
      legend.box = "horizontal"
    ) +
    labs(x = "log2 Fold Change", y = "-log10(P Value)", 
         title = "Driver vs Suppressor Genes in TCGA-LUAD",
         subtitle = paste("Drivers:", nrow(driver_genes), "| Suppressors:", nrow(suppressor_genes))) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5, color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5, color = "gray50") +
    guides(color = guide_legend(override.aes = list(size = 4)),
           shape = guide_legend(override.aes = list(size = 4)))
  
  ggsave("results/plots/driver_vs_suppressor_volcano.pdf", volcano_plot, width = 16, height = 12)
  
  # 2. Side-by-side comparison plots
  if (nrow(driver_genes) > 0 && nrow(suppressor_genes) > 0) {
    
    # Box plot comparing log2FC distributions
    boxplot_data <- rbind(
      data.frame(Gene_Type = "Driver", logFC = driver_genes$logFC),
      data.frame(Gene_Type = "Suppressor", logFC = suppressor_genes$logFC)
    )
    
    boxplot_comparison <- ggplot(boxplot_data, aes(x = Gene_Type, y = logFC, fill = Gene_Type)) +
      geom_boxplot(alpha = 0.7, outlier.size = 2) +
      geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
      scale_fill_manual(values = c("Driver" = "#E74C3C", "Suppressor" = "#3498DB")) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none"
      ) +
      labs(x = "Gene Type", y = "log2 Fold Change", 
           title = "Expression Change Distribution: Drivers vs Suppressors") +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)
    
    ggsave("results/plots/driver_vs_suppressor_boxplot.pdf", boxplot_comparison, width = 8, height = 6)
    
    # Regulation status comparison
    reg_summary <- rbind(
      data.frame(Gene_Type = "Driver", Regulation = names(table(driver_genes$change)), 
                 Count = as.numeric(table(driver_genes$change))),
      data.frame(Gene_Type = "Suppressor", Regulation = names(table(suppressor_genes$change)), 
                 Count = as.numeric(table(suppressor_genes$change)))
    )
    
    reg_comparison <- ggplot(reg_summary, aes(x = Gene_Type, y = Count, fill = Regulation)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
      scale_fill_manual(values = c("down" = "#2B6688", "up" = "#F1A93B", "stable" = "#A8ACB9")) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      ) +
      labs(x = "Gene Type", y = "Number of Genes", 
           title = "Regulation Status: Drivers vs Suppressors") +
      geom_text(aes(label = Count), position = position_dodge(width = 0.9), 
                vjust = -0.5, fontface = "bold")
    
    ggsave("results/plots/driver_vs_suppressor_regulation.pdf", reg_comparison, width = 10, height = 6)
  }
  
  # 3. Individual heatmaps for each group
  if (nrow(driver_genes) > 0) {
    driver_heatmap_data <- driver_genes[order(driver_genes$logFC, decreasing = TRUE), ]
    driver_heatmap_data$Gene_Order <- factor(driver_heatmap_data$gene_symbol, 
                                             levels = driver_heatmap_data$gene_symbol)
    
    driver_heatmap <- ggplot(driver_heatmap_data, aes(x = 1, y = Gene_Order, fill = logFC)) +
      geom_tile(color = "white", size = 0.5) +
      scale_fill_gradient2(low = "#2B6688", mid = "white", high = "#F1A93B",
                           midpoint = 0, name = "log2FC") +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.y = element_text(size = 9)
      ) +
      labs(title = "Driver Genes Expression Heatmap",
           subtitle = "Ordered by log2 Fold Change",
           y = "Gene Symbol")
    
    ggsave("results/plots/driver_genes_heatmap.pdf", driver_heatmap, 
           width = 6, height = max(6, nrow(driver_genes) * 0.3))
  }
  
  if (nrow(suppressor_genes) > 0) {
    suppressor_heatmap_data <- suppressor_genes[order(suppressor_genes$logFC, decreasing = TRUE), ]
    suppressor_heatmap_data$Gene_Order <- factor(suppressor_heatmap_data$gene_symbol, 
                                                 levels = suppressor_heatmap_data$gene_symbol)
    
    suppressor_heatmap <- ggplot(suppressor_heatmap_data, aes(x = 1, y = Gene_Order, fill = logFC)) +
      geom_tile(color = "white", size = 0.5) +
      scale_fill_gradient2(low = "#2B6688", mid = "white", high = "#F1A93B",
                           midpoint = 0, name = "log2FC") +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.y = element_text(size = 9)
      ) +
      labs(title = "Suppressor Genes Expression Heatmap",
           subtitle = "Ordered by log2 Fold Change",
           y = "Gene Symbol")
    
    ggsave("results/plots/suppressor_genes_heatmap.pdf", suppressor_heatmap, 
           width = 6, height = max(6, nrow(suppressor_genes) * 0.3))
  }
  
  cat("Group-specific visualizations saved:\n")
  cat("- results/plots/driver_vs_suppressor_volcano.pdf\n")
  cat("- results/plots/driver_vs_suppressor_boxplot.pdf\n")
  cat("- results/plots/driver_vs_suppressor_regulation.pdf\n")
  cat("- results/plots/driver_genes_heatmap.pdf\n")
  cat("- results/plots/suppressor_genes_heatmap.pdf\n\n")
  
  return(list(volcano = volcano_plot, boxplot = boxplot_comparison, 
              regulation = reg_comparison, driver_heatmap = driver_heatmap, 
              suppressor_heatmap = suppressor_heatmap))
}

# Run the grouped analysis
cat("Starting Driver vs Suppressor gene analysis...\n\n")
grouped_results <- analyze_gene_groups(driver_genes, suppressor_genes)

# Create visualizations
if (!is.null(grouped_results)) {
  cat("Creating grouped visualizations...\n")
  grouped_plots <- create_grouped_visualizations(grouped_results)
  
  # Final summary
  cat("=== FINAL GROUPED ANALYSIS SUMMARY ===\n")
  cat("Driver genes found:", nrow(grouped_results$driver_genes), "/", length(driver_genes), "\n")
  cat("Suppressor genes found:", nrow(grouped_results$suppressor_genes), "/", length(suppressor_genes), "\n")
  
  if (!is.null(grouped_results$statistics)) {
    cat("Mean log2FC - Drivers:", round(grouped_results$statistics$driver_mean_fc, 3), "\n")
    cat("Mean log2FC - Suppressors:", round(grouped_results$statistics$suppressor_mean_fc, 3), "\n")
    cat("Statistical difference (p-value):", format(grouped_results$statistics$t_test$p.value, scientific = TRUE, digits = 3), "\n")
  }
  
  cat("\nAnalysis complete! Check the Excel file and plots for detailed results.\n")
}
