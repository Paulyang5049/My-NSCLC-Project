getwd()
setwd("/Users/yangpaul/Desktop/")

# Load required libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(scales)
library(cowplot)

# Updated theme for publication-quality plots
theme_pub <- theme_minimal() +
    theme(
        axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11, face = "bold"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        strip.text = element_text(size = 11, face = "bold"),
        plot.margin = margin(10, 10, 10, 10)
    )

# Define unified color palette for all plots
unified_colors <- c("#2166AC", "#4393C3", "#92C5DE", "#F7F7F7", "#F4A582", "#D6604D", "#B2182B")

# Read CC data - standardized columns
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
    data$GO_Type <- "Cellular Component"
    
    return(data)
}

# Read MF data - standardized columns
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
    data$GO_Type <- "Molecular Function"
    
    return(data)
}

# Read BP data - standardized columns (added missing columns)
read_go_bp <- function() {
    data <- data.frame(
        GO_ID = paste0("GO:BP_", 1:16), # Placeholder GO IDs for BP terms
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
        GeneRatio = c("9/48", "9/48", "9/48", "9/48", "8/48", "8/48", "8/48", "8/48", 
                      "7/48", "7/48", "7/48", "7/48", "7/48", "7/48", "6/48", "5/48"),
        BgRatio = c("242/18986", "367/18986", "405/18986", "467/18986", "344/18986", 
                    "374/18986", "446/18986", "486/18986", "237/18986", "325/18986",
                    "331/18986", "403/18986", "438/18986", "490/18986", "252/18986", "296/18986"),
        RichFactor = c(0.037190083, 0.024523161, 0.022222222, 0.019271949, 0.023255814,
                       0.021390374, 0.01793722, 0.016460905, 0.029535865, 0.021538462,
                       0.021148036, 0.017369727, 0.015981735, 0.014285714, 0.023809524, 0.016891892),
        FoldEnrichment = c(14.71022727, 9.69993188, 8.789814815, 7.622858672, 9.198643411,
                           8.460784314, 7.094917788, 6.510973937, 11.68266526, 8.519358974,
                           8.364929507, 6.870450786, 6.321442161, 5.650595238, 9.41765873, 6.681447072),
        zScore = rep(NA, 16), # Placeholder values
        pvalue = c(8.32692E-09, 2.94605E-07, 6.72237E-07, 2.18288E-06, 2.14821E-06,
                   3.98872E-06, 1.44111E-05, 2.66978E-05, 2.05856E-06, 1.62739E-05,
                   1.83091E-05, 6.40505E-05, 0.000107798, 0.00021523, 3.96323E-05, 0.000880951),
        p.adjust = c(1.67704E-05, 0.000296668, 0.000451295, 0.000662825, 0.000662825,
                     0.000953504, 0.00229869, 0.002688464, 0.000662825, 0.00229869,
                     0.00229869, 0.004717053, 0.006846689, 0.009524038, 0.003336748, 0.02245868),
        qvalue = rep(NA, 16), # Placeholder values
        geneID = rep(NA, 16), # Placeholder values
        Count = c(9, 9, 9, 9, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 6, 5),
        stringsAsFactors = FALSE
    )
    
    data$neg_log10_padj <- -log10(data$p.adjust)
    data$GO_Type <- "Biological Process"
    
    return(data)
}

# Function to create unified GO plot with all terms on same axis
create_unified_go_plot <- function(go_cc, go_mf, go_bp, cc_top = 3, mf_top = 2, bp_top = 10) {
    
    # Select top terms from each category
    cc_data <- go_cc %>% 
        arrange(p.adjust) %>% 
        head(cc_top) %>%
        mutate(
            Term_wrapped = str_wrap(Term, width = 40),
            Category_Term = paste0(GO_Type, ": ", Term_wrapped)
        )
    
    mf_data <- go_mf %>% 
        arrange(p.adjust) %>% 
        head(mf_top) %>%
        mutate(
            Term_wrapped = str_wrap(Term, width = 40),
            Category_Term = paste0(GO_Type, ": ", Term_wrapped)
        )
    
    bp_data <- go_bp %>% 
        arrange(p.adjust) %>% 
        head(bp_top) %>%
        mutate(
            Term_wrapped = str_wrap(Term, width = 40),
            Category_Term = paste0(GO_Type, ": ", Term_wrapped)
        )
    
    # Combine all data
    combined_data <- rbind(cc_data, mf_data, bp_data)
    
    # Order by GO type and significance within type
    combined_data <- combined_data %>%
        arrange(GO_Type, p.adjust) %>%
        mutate(
            # Create ordered factor for plotting
            GO_Type = factor(GO_Type, levels = c("Biological Process", "Molecular Function", "Cellular Component")),
            Category_Term = factor(Category_Term, levels = rev(Category_Term)),
            # Create a grouping variable for visual separation
            Group_order = as.numeric(GO_Type)
        )
    
    # Calculate max for consistent scaling
    max_log10_padj <- max(combined_data$neg_log10_padj, na.rm = TRUE)
    
    # Create the plot
    p <- ggplot(combined_data, aes(x = neg_log10_padj, y = Category_Term)) +
        geom_col(aes(fill = neg_log10_padj), width = 0.7, color = "black", linewidth = 0.2) +
        scale_fill_gradientn(
            colors = unified_colors,
            name = expression(-log[10](italic(P)[adj])),
            limits = c(0, max_log10_padj),
            guide = guide_colorbar(
                title.position = "top",
                title.hjust = 0.5,
                barwidth = 1,
                barheight = 8
            )
        ) +
        scale_x_continuous(
            limits = c(0, max_log10_padj * 1.05),
            expand = c(0, 0),
            breaks = pretty_breaks(n = 6)
        ) +
        # Add category separators
        geom_hline(yintercept = c(bp_top + 0.5, bp_top + mf_top + 0.5), 
                   linetype = "solid", color = "grey60", alpha = 0.5, linewidth = 0.8) +
        # Add significance threshold line
        {if (max_log10_padj > 1.301) {
            geom_vline(xintercept = 1.301, linetype = "dashed", 
                       color = "red", alpha = 0.7, linewidth = 0.8)
        }} +
        labs(
            title = "Gene Ontology Enrichment Analysis",
            x = expression(-log[10](italic(P)[adjusted])),
            y = NULL
        ) +
        theme_pub +
        theme(
            legend.position = "right",
            axis.text.y = element_text(size = 10, hjust = 1),
            panel.grid.major.y = element_blank(),
            plot.title = element_text(size = 16)
        ) +
        # Add category labels on the right
        annotate("text", x = max_log10_padj * 1.02, 
                 y = c(bp_top/2, bp_top + mf_top/2, bp_top + mf_top + cc_top/2),
                 label = c("BP", "MF", "CC"),
                 angle = 90, vjust = 0.5, hjust = 0.5,
                 size = 3.5, fontface = "bold", color = "grey30")
    
    return(p)
}

# Read pathway data (keeping the existing functions)
read_kegg_data <- function(file_path) {
    data <- read.csv(file_path, stringsAsFactors = FALSE)
    data <- data %>%
        mutate(
            neg_log10_padj = -log10(p.adjust),
            gene_count = as.numeric(str_extract(GeneRatio, "^\\d+")),
            Pathway = str_wrap(Description, width = 40)
        ) %>%
        arrange(p.adjust) %>%
        mutate(Pathway = factor(Pathway, levels = rev(Pathway)))
    return(data)
}

read_reactome_data <- function(file_path) {
    data <- read.csv(file_path, stringsAsFactors = FALSE)
    data <- data %>%
        mutate(
            neg_log10_padj = -log10(p.adjust),
            gene_count = as.numeric(str_extract(GeneRatio, "^\\d+")),
            Pathway = str_wrap(Description, width = 40)
        ) %>%
        arrange(p.adjust) %>%
        mutate(Pathway = factor(Pathway, levels = rev(Pathway)))
    return(data)
}

# Function to create unified pathway bubble plot
create_unified_bubble_plot <- function(kegg_data, reactome_data, max_vals = NULL) {
    
    # Prepare data
    kegg_top <- kegg_data %>% head(10) %>% mutate(Database = "KEGG")
    reactome_top <- reactome_data %>% head(10) %>% mutate(Database = "Reactome")
    
    # Combine data
    combined_data <- rbind(kegg_top, reactome_top) %>%
        mutate(
            Database = factor(Database, levels = c("KEGG", "Reactome")),
            Plot_order = row_number(),
            Pathway_short = str_wrap(Description, width = 45)
        ) %>%
        arrange(Database, p.adjust) %>%
        mutate(
            Pathway_short = factor(Pathway_short, levels = rev(Pathway_short))
        )
    
    # Set unified scales
    if(is.null(max_vals)) {
        max_gene_count <- max(combined_data$gene_count, na.rm = TRUE)
        max_fold_enrich <- max(combined_data$FoldEnrichment, na.rm = TRUE)
        max_rich_factor <- max(combined_data$RichFactor, na.rm = TRUE)
    } else {
        max_gene_count <- max_vals$gene_count
        max_fold_enrich <- max_vals$fold_enrich
        max_rich_factor <- max_vals$rich_factor
    }
    
    p <- ggplot(combined_data, aes(x = RichFactor, y = Pathway_short)) +
        geom_point(aes(size = gene_count, color = FoldEnrichment), alpha = 0.8) +
        scale_color_gradientn(
            colors = unified_colors,
            name = "Fold\nEnrichment",
            trans = "log10",
            breaks = c(2, 5, 10, 20, 50),
            labels = c("2", "5", "10", "20", "50"),
            limits = c(1, max_fold_enrich),
            guide = guide_colorbar(
                title.position = "top",
                title.hjust = 0.5,
                barwidth = 1,
                barheight = 6
            )
        ) +
        scale_size_continuous(
            name = "Gene\nCount",
            range = c(3, 10),
            limits = c(1, max_gene_count),
            breaks = pretty_breaks(n = 4),
            guide = guide_legend(
                title.position = "top",
                title.hjust = 0.5,
                override.aes = list(alpha = 0.8)
            )
        ) +
        scale_x_continuous(
            breaks = pretty_breaks(n = 5),
            limits = c(0, max_rich_factor * 1.05),
            expand = c(0.02, 0.02)
        ) +
        facet_wrap(~Database, scales = "free_y", ncol = 1) +
        labs(
            title = "Pathway Enrichment Analysis",
            x = "Rich Factor",
            y = NULL
        ) +
        theme_pub +
        theme(
            axis.text.y = element_text(size = 9),
            strip.text = element_text(size = 12, face = "bold"),
            legend.position = "right",
            legend.box = "vertical",
            panel.spacing = unit(1, "lines"),
            plot.title = element_text(size = 16)
        )
    
    return(p)
}

# Read all data
go_cc <- read_go_cc()
go_mf <- read_go_mf()
go_bp <- read_go_bp()

# Create unified GO plot
go_unified <- create_unified_go_plot(go_cc, go_mf, go_bp, 
                                     cc_top = 3, mf_top = 2, bp_top = 10)

# Display GO plot
print(go_unified)

# Save GO unified plot
ggsave("GO_Combined_Analysis_Unified.pdf", go_unified, 
       width = 14, height = 10, device = "pdf", bg = "white")

# Create pathway plots (only if CSV files exist)
if(file.exists("kegg.csv") && file.exists("reactome.csv")) {
    kegg_data <- read_kegg_data("kegg.csv")
    reactome_data <- read_reactome_data("reactome.csv")
    
    # Create unified pathway bubble plot
    pathway_plot <- create_unified_bubble_plot(kegg_data, reactome_data)
    
    # Save pathway plot
    ggsave("Pathway_Combined_Analysis_Unified.pdf", pathway_plot, 
           width = 14, height = 12, device = "pdf", bg = "white")
    
    print("Both GO and Pathway unified plots created successfully!")
    
    # Print summary
    cat("\nSummary:\n")
    cat("GO Combined plot: GO_Combined_Analysis_Unified.pdf\n")
    cat("Pathway Combined plot: Pathway_Combined_Analysis_Unified.pdf\n")
    cat("KEGG pathways:", nrow(kegg_data), "\n")
    cat("Reactome pathways:", nrow(reactome_data), "\n")
} else {
    print("GO unified plot created successfully!")
    print("Pathway CSV files not found - only GO plot generated.")
    cat("\nSummary:\n")
    cat("GO Combined plot: GO_Combined_Analysis_Unified.pdf\n")
}

cat("GO CC terms:", nrow(go_cc), "\n")
cat("GO MF terms:", nrow(go_mf), "\n") 
cat("GO BP terms:", nrow(go_bp), "\n")
