getwd()
setwd("/Users/yangpaul/Desktop/")

# Load required libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(scales)

# Updated theme for publication-quality plots
theme_pub <- theme_minimal() +
    theme(
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        strip.text = element_text(size = 12, face = "bold")
    )

# Read CC data
read_go_cc <- function() {
    data <- data.frame(
        GO_ID = c("GO:0009925", "GO:0005902", "GO:1990204"),
        Term = c("basal plasma membrane", "microvillus", "oxidoreductase complex"),
        GeneRatio = c("5/48", "4/48", "4/48"),
        BgRatio = c("298/19960", "100/19960", "150/19960"),
        RichFactor = c(0.016778523, 0.04, 0.026666667),
        FoldEnrichment = c(6.977069351, 16.63333333, 11.08888889),
        zScore = c(5.104052667, 7.694750311, 6.089467218),
        pvalue = c(0.000726265, 9.74385E-05, 0.000461081),
        p.adjust = c(0.028104741, 0.017149173, 0.027050069),
        qvalue = c(0.020338957, 0.012410586, 0.019575707),
        geneID = c("AQP4/PROM2/AURKA/CD44/CA9", "PROM2/LRRK2/CD44/CA9", "COX7A1/DUOX1/COX4I2/CYBB"),
        Count = c(5, 4, 4),
        stringsAsFactors = FALSE
    )
    
    data$neg_log10_padj <- -log10(data$p.adjust)
    data$GO_Type <- "CC"
    
    return(data)
}

# Read MF data
read_go_mf <- function() {
    data <- data.frame(
        GO_ID = c("GO:0015631", "GO:0001216"),
        Term = c("tubulin binding", "DNA-binding transcription activator activity"),
        GeneRatio = c("7/48", "6/48"),
        BgRatio = c("392/18737", "487/18737"),
        RichFactor = c(0.017857143, 0.012320329),
        FoldEnrichment = c(6.970610119, 4.809291581),
        zScore = c(6.054369376, 4.316614732),
        pvalue = c(5.84445E-05, 0.00145417),
        p.adjust = c(0.013442234, 0.041989249),
        qvalue = c(0.010089366, 0.031515958),
        geneID = c("CRYAB/SNCA/S100A8/GABARAPL1/GJA1/KIF20A/LRRK2", "EPAS1/ATF3/EGR1/ETV4/TFAP2A/JUN"),
        Count = c(7, 6),
        stringsAsFactors = FALSE
    )
    
    data$neg_log10_padj <- -log10(data$p.adjust)
    data$GO_Type <- "MF"
    
    return(data)
}

# Read BP data with unique terms
read_go_bp <- function() {
    data <- data.frame(
        Term = c("response to oxidative stress", 
                 "cellular response to chemical stimulus", 
                 "response to oxygen-containing compound", 
                 "cellular response to organic substance",
                 "regulation of cell death", 
                 "cellular response to endogenous stimulus",
                 "regulation of apoptotic process", 
                 "cellular response to stimulus",
                 "positive regulation of cell death", 
                 "regulation of programmed cell death",
                 "cellular response to oxygen compound", 
                 "cellular response to stress",
                 "response to endogenous stimulus", 
                 "response to organic substance",
                 "positive regulation of apoptotic process", 
                 "regulation of metabolic process"),
        BgRatio = c("242/18986", "367/18986", "405/18986", "467/18986", "344/18986", 
                    "374/18986", "446/18986", "486/18986", "237/18986", "325/18986",
                    "331/18986", "403/18986", "438/18986", "490/18986", "252/18986", "296/18986"),
        RichFactor = c(0.037190083, 0.024523161, 0.022222222, 0.019271949, 0.023255814,
                       0.021390374, 0.01793722, 0.016460905, 0.029535865, 0.021538462,
                       0.021148036, 0.017369727, 0.015981735, 0.014285714, 0.023809524, 0.016891892),
        FoldEnrichment = c(14.71022727, 9.69993188, 8.789814815, 7.622858672, 9.198643411,
                           8.460784314, 7.094917788, 6.510973937, 11.68266526, 8.519358974,
                           8.364929507, 6.870450786, 6.321442161, 5.650595238, 9.41765873, 6.681447072),
        pvalue = c(8.32692E-09, 2.94605E-07, 6.72237E-07, 2.18288E-06, 2.14821E-06,
                   3.98872E-06, 1.44111E-05, 2.66978E-05, 2.05856E-06, 1.62739E-05,
                   1.83091E-05, 6.40505E-05, 0.000107798, 0.00021523, 3.96323E-05, 0.000880951),
        p.adjust = c(1.67704E-05, 0.000296668, 0.000451295, 0.000662825, 0.000662825,
                     0.000953504, 0.00229869, 0.002688464, 0.000662825, 0.00229869,
                     0.00229869, 0.004717053, 0.006846689, 0.009524038, 0.003336748, 0.02245868),
        Count = c(9, 9, 9, 9, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 6, 5),
        stringsAsFactors = FALSE
    )
    
    data$neg_log10_padj <- -log10(data$p.adjust)
    data$GO_Type <- "BP"
    
    return(data)
}

# Fixed function to create publication-quality GO plot
create_go_plot <- function(data, title, top_n = 10) {
    # Select top terms by significance and ensure uniqueness
    plot_data <- data %>%
        arrange(p.adjust) %>%
        head(top_n) %>%
        mutate(
            # Create unique identifiers first
            Term_ID = paste0(Term, "_", row_number()),
            # Wrap text
            Term_wrapped = str_wrap(Term, width = 40),
            # Create display labels (without ID suffix)
            Display_Term = Term_wrapped
        ) %>%
        # Ensure no duplicates in display terms
        mutate(
            Display_Term = ifelse(duplicated(Display_Term), 
                                  paste0(Display_Term, " (", row_number(), ")"), 
                                  Display_Term)
        ) %>%
        mutate(
            Display_Term = factor(Display_Term, levels = rev(Display_Term))
        )
    
    # Create color palette
    colors <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#F4A582", "#D6604D", "#B2182B"))(6)
    
    p <- ggplot(plot_data, aes(x = Display_Term, y = neg_log10_padj)) +
        geom_col(aes(fill = neg_log10_padj), width = 0.7, color = "black", linewidth = 0.3) +
        scale_fill_gradientn(
            colors = colors,
            name = expression(-log[10](italic(P)[adj])),
            guide = guide_colorbar(
                title.position = "top",
                title.hjust = 0.5,
                barwidth = 0.8,
                barheight = 8
            )
        ) +
        scale_y_continuous(
            expand = expansion(mult = c(0, 0.05)),
            breaks = pretty_breaks(n = 6)
        ) +
        coord_flip() +
        labs(
            title = title,
            x = "GO Terms",
            y = expression(-log[10](italic(P)[adjusted])),
            caption = paste("Top", nrow(plot_data), "most significant terms")
        ) +
        theme_pub +
        theme(
            legend.position = "right",
            plot.caption = element_text(size = 10, color = "grey50", hjust = 1)
        )
    
    # Add significance threshold line
    if (max(plot_data$neg_log10_padj, na.rm = TRUE) > 1.301) {
        p <- p + geom_hline(yintercept = 1.301, linetype = "dashed", 
                            color = "red", alpha = 0.7, linewidth = 0.8)
    }
    
    return(p)
}

# Read data
go_cc <- read_go_cc()
go_mf <- read_go_mf()
go_bp <- read_go_bp()

# Create individual plots
p_cc <- create_go_plot(go_cc, "Cellular Component", top_n = min(10, nrow(go_cc)))
p_mf <- create_go_plot(go_mf, "Molecular Function", top_n = min(10, nrow(go_mf)))
p_bp <- create_go_plot(go_bp, "Biological Process", top_n = min(15, nrow(go_bp)))

# Display plots
print(p_cc)
print(p_mf)
print(p_bp)

# Save high-resolution plots for publication
ggsave("GO_Cellular_Component.png", p_cc, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("GO_Molecular_Function.png", p_mf, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("GO_Biological_Process.png", p_bp, width = 12, height = 8, dpi = 300, bg = "white")

# Save as PDF for vector graphics (recommended for publications)
ggsave("GO_Cellular_Component.pdf", p_cc, width = 10, height = 6, device = "pdf")
ggsave("GO_Molecular_Function.pdf", p_mf, width = 10, height = 6, device = "pdf")
ggsave("GO_Biological_Process.pdf", p_bp, width = 12, height = 8, device = "pdf")

# Optional: Create a combined plot
library(cowplot)
combined_plot <- plot_grid(p_cc, p_mf, p_bp, 
                           labels = c("A", "B", "C"), 
                           label_size = 16,
                           ncol = 1, 
                           rel_heights = c(0.8, 0.8, 1.2))

ggsave("GO_Combined_Analysis.pdf", combined_plot, width = 12, height = 16, device = "pdf")

print("All plots created successfully!")

# Print summary
cat("\nSummary of GO Analysis Results:\n")
cat("Cellular Component terms:", nrow(go_cc), "\n")
cat("Molecular Function terms:", nrow(go_mf), "\n")
cat("Biological Process terms:", nrow(go_bp), "\n")





# Set theme for publication-quality plots
theme_pub <- theme_minimal() +
    theme(
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "right"
    )

# Read KEGG data
read_kegg_data <- function(file_path) {
    data <- read.csv(file_path, stringsAsFactors = FALSE)
    
    # Clean and prepare data
    data <- data %>%
        mutate(
            # Calculate -log10(p.adjust) for color scale
            neg_log10_padj = -log10(p.adjust),
            # Extract gene ratio numerator for bubble size
            gene_count = as.numeric(str_extract(GeneRatio, "^\\d+")),
            # Create pathway labels (remove pathway ID for cleaner look)
            Pathway = str_wrap(Description, width = 50),
            # Create category-subcategory labels
            Category_full = paste(category, subcategory, sep = " - ")
        ) %>%
        # Order by significance
        arrange(p.adjust) %>%
        # Add row number for ordering in plot
        mutate(
            Pathway = factor(Pathway, levels = rev(Pathway))
        )
    
    return(data)
}

# Read Reactome data
read_reactome_data <- function(file_path) {
    data <- read.csv(file_path, stringsAsFactors = FALSE)
    
    # Clean and prepare data
    data <- data %>%
        mutate(
            # Calculate -log10(p.adjust) for color scale
            neg_log10_padj = -log10(p.adjust),
            # Extract gene ratio numerator for bubble size
            gene_count = as.numeric(str_extract(GeneRatio, "^\\d+")),
            # Create pathway labels
            Pathway = str_wrap(Description, width = 50)
        ) %>%
        # Order by significance
        arrange(p.adjust) %>%
        # Add row number for ordering in plot
        mutate(
            Pathway = factor(Pathway, levels = rev(Pathway))
        )
    
    return(data)
}

# Function to create bubble plot
create_bubble_plot <- function(data, title, color_var = "neg_log10_padj", 
                               size_var = "gene_count", x_var = "RichFactor",
                               top_n = 15) {
    
    # Select top pathways by significance
    plot_data <- data %>%
        head(top_n)
    
    # Create color palette
    colors <- colorRampPalette(c("#3182BD", "#6BAED6", "#9ECAE1", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#D94701"))(7)
    
    p <- ggplot(plot_data, aes_string(x = x_var, y = "Pathway")) +
        geom_point(aes_string(size = size_var, color = color_var), alpha = 0.8) +
        scale_color_gradientn(
            colors = colors,
            name = expression(-log[10](italic(P)[adj])),
            guide = guide_colorbar(
                title.position = "top",
                title.hjust = 0.5,
                barwidth = 1,
                barheight = 8
            )
        ) +
        scale_size_continuous(
            name = "Gene Count",
            range = c(3, 12),
            breaks = pretty_breaks(n = 4),
            guide = guide_legend(
                title.position = "top",
                title.hjust = 0.5,
                override.aes = list(alpha = 0.8)
            )
        ) +
        scale_x_continuous(
            breaks = pretty_breaks(n = 6),
            expand = expansion(mult = c(0.1, 0.1))
        ) +
        labs(
            title = title,
            x = "Rich Factor",
            y = "Pathways",
            caption = paste("Top", nrow(plot_data), "most significant pathways")
        ) +
        theme_pub +
        theme(
            axis.text.y = element_text(size = 10),
            plot.caption = element_text(size = 10, color = "grey50", hjust = 1),
            legend.box = "horizontal",
            legend.box.just = "top"
        )
    
    # Add significance threshold line if applicable
    if (max(plot_data[[color_var]], na.rm = TRUE) > 1.301) {
        # Add a subtle background color for significant pathways
        sig_pathways <- plot_data[plot_data$p.adjust < 0.05, ]
        if (nrow(sig_pathways) > 0) {
            p <- p + 
                geom_rect(data = sig_pathways, 
                          aes(xmin = -Inf, xmax = Inf, 
                              ymin = as.numeric(Pathway) - 0.4, 
                              ymax = as.numeric(Pathway) + 0.4),
                          fill = "lightblue", alpha = 0.1, inherit.aes = FALSE)
        }
    }
    
    return(p)
}

# Function to create enhanced bubble plot with fold enrichment
create_enhanced_bubble_plot <- function(data, title, top_n = 15) {
    
    # Select top pathways by significance
    plot_data <- data %>%
        head(top_n)
    
    # Create color palette for fold enrichment
    colors <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#F4A582", "#D6604D", "#B2182B"))(6)
    
    p <- ggplot(plot_data, aes(x = RichFactor, y = Pathway)) +
        geom_point(aes(size = gene_count, color = FoldEnrichment), alpha = 0.8) +
        scale_color_gradientn(
            colors = colors,
            name = "Fold\nEnrichment",
            trans = "log10",
            breaks = c(2, 5, 10, 20, 50),
            labels = c("2", "5", "10", "20", "50"),
            guide = guide_colorbar(
                title.position = "top",
                title.hjust = 0.5,
                barwidth = 1,
                barheight = 8
            )
        ) +
        scale_size_continuous(
            name = "Gene\nCount",
            range = c(3, 12),
            breaks = pretty_breaks(n = 4),
            guide = guide_legend(
                title.position = "top",
                title.hjust = 0.5,
                override.aes = list(alpha = 0.8)
            )
        ) +
        scale_x_continuous(
            breaks = pretty_breaks(n = 6),
            expand = expansion(mult = c(0.1, 0.1))
        ) +
        labs(
            title = title,
            x = "Rich Factor",
            y = "Pathways",
            caption = paste("Top", nrow(plot_data), "most significant pathways")
        ) +
        theme_pub +
        theme(
            axis.text.y = element_text(size = 10),
            plot.caption = element_text(size = 10, color = "grey50", hjust = 1),
            legend.box = "horizontal",
            legend.box.just = "top"
        )
    
    return(p)
}

# Read the data
kegg_data <- read_kegg_data("kegg.csv")
reactome_data <- read_reactome_data("reactome.csv")

# Create bubble plots - Version 1: Colored by -log10(p.adjust)
p_kegg_padj <- create_bubble_plot(kegg_data, "KEGG Pathway Enrichment Analysis", 
                                  top_n = min(15, nrow(kegg_data)))

p_reactome_padj <- create_bubble_plot(reactome_data, "Reactome Pathway Enrichment Analysis", 
                                      top_n = min(15, nrow(reactome_data)))

# Create bubble plots - Version 2: Colored by Fold Enrichment
p_kegg_fold <- create_enhanced_bubble_plot(kegg_data, "KEGG Pathway Enrichment Analysis", 
                                           top_n = min(15, nrow(kegg_data)))

p_reactome_fold <- create_enhanced_bubble_plot(reactome_data, "Reactome Pathway Enrichment Analysis", 
                                               top_n = min(15, nrow(reactome_data)))

# Display plots
print("=== KEGG Bubble Plot (colored by p.adjust) ===")
print(p_kegg_padj)

print("=== Reactome Bubble Plot (colored by p.adjust) ===")
print(p_reactome_padj)

print("=== KEGG Bubble Plot (colored by Fold Enrichment) ===")
print(p_kegg_fold)

print("=== Reactome Bubble Plot (colored by Fold Enrichment) ===")
print(p_reactome_fold)

# Save high-resolution plots for publication
# Version 1: p.adjust coloring
ggsave("KEGG_Bubble_Plot_pAdj.png", p_kegg_padj, width = 12, height = 8, dpi = 300, bg = "white")
ggsave("Reactome_Bubble_Plot_pAdj.png", p_reactome_padj, width = 12, height = 8, dpi = 300, bg = "white")

ggsave("KEGG_Bubble_Plot_pAdj.pdf", p_kegg_padj, width = 12, height = 8, device = "pdf")
ggsave("Reactome_Bubble_Plot_pAdj.pdf", p_reactome_padj, width = 12, height = 8, device = "pdf")

# Version 2: Fold enrichment coloring
ggsave("KEGG_Bubble_Plot_FoldEnrich.png", p_kegg_fold, width = 12, height = 8, dpi = 300, bg = "white")
ggsave("Reactome_Bubble_Plot_FoldEnrich.png", p_reactome_fold, width = 12, height = 8, dpi = 300, bg = "white")

ggsave("KEGG_Bubble_Plot_FoldEnrich.pdf", p_kegg_fold, width = 12, height = 8, device = "pdf")
ggsave("Reactome_Bubble_Plot_FoldEnrich.pdf", p_reactome_fold, width = 12, height = 8, device = "pdf")

# Create combined plots using cowplot
library(cowplot)

# Combined plot with p.adjust coloring
combined_padj <- plot_grid(p_kegg_padj, p_reactome_padj, 
                           labels = c("A", "B"), 
                           label_size = 16,
                           ncol = 1,
                           rel_heights = c(1, 1))

# Combined plot with fold enrichment coloring
combined_fold <- plot_grid(p_kegg_fold, p_reactome_fold, 
                           labels = c("A", "B"), 
                           label_size = 16,
                           ncol = 1,
                           rel_heights = c(1, 1))

# Save combined plots
ggsave("Combined_Pathway_Bubble_pAdj.pdf", combined_padj, width = 12, height = 14, device = "pdf")
ggsave("Combined_Pathway_Bubble_FoldEnrich.pdf", combined_fold, width = 12, height = 14, device = "pdf")

print("All bubble plots saved successfully!")

# Print summary statistics
cat("\nSummary of Pathway Analysis Results:\n")
cat("KEGG pathways analyzed:", nrow(kegg_data), "\n")
cat("Reactome pathways analyzed:", nrow(reactome_data), "\n")
cat("KEGG significant pathways (p.adj < 0.05):", sum(kegg_data$p.adjust < 0.05), "\n")
cat("Reactome significant pathways (p.adj < 0.05):", sum(reactome_data$p.adjust < 0.05), "\n")

# Create a summary table of top pathways
cat("\nTop 5 KEGG Pathways:\n")
print(kegg_data[1:min(5, nrow(kegg_data)), c("Description", "gene_count", "p.adjust", "FoldEnrichment")])

cat("\nTop 5 Reactome Pathways:\n")
print(reactome_data[1:min(5, nrow(reactome_data)), c("Description", "gene_count", "p.adjust", "FoldEnrichment")])


