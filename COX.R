# Load required libraries
library(survival)    # For Cox proportional hazards models
library(ggplot2)     # For creating the forest plot
library(GEOquery)    # For parsing GEO data
library(dplyr)       # For data manipulation
library(broom)       # For tidying model outputs
gc()
rm(list = ls())
# Set working directory (modify as needed)
setwd("~/bioinformatics/gse13213")

# Step 1: Load and parse the GEO data
gse <- getGEO("GSE13213", destdir = getwd())
expr_data <- exprs(gse[[1]])
pheno_data <- pData(gse[[1]])

# Step 2: Prepare survival data and covariates
surv_data <- pheno_data[, c("characteristics_ch1.9", "characteristics_ch1.8", 
                            "characteristics_ch1.2", "characteristics_ch1.7", 
                            "characteristics_ch1.5")]

colnames(surv_data) <- c("time_raw", "status_raw", "age_raw", "stage_raw", "smoking_raw")

# Clean and convert data types
surv_data$time <- as.numeric(gsub("Survival \\(days\\): ", "", surv_data$time_raw))
surv_data$status <- ifelse(surv_data$status_raw == "Status: Dead", 1, 0)  # 1=event, 0=censored
surv_data$age <- as.numeric(gsub("Age: ", "", surv_data$age_raw))
surv_data$stage <- as.factor(gsub("Stage \\(Pathological \\): ", "", surv_data$stage_raw))
surv_data$smoking <- as.numeric(gsub("Smoking \\(BI\\): ", "", surv_data$smoking_raw))

# Select clean columns and remove NA values
surv_data_clean <- surv_data[, c("time", "status", "age", "stage", "smoking")]
surv_data_clean <- na.omit(surv_data_clean)

# Print status distribution to verify coding
print("Status distribution:")
print(table(surv_data_clean$status))

# Step 3: Define genes of interest and get probe mappings
genes_of_interest <- c("KIF20A", "PYCR1", "AURKA", "UHRF1")

# Use either platform table mapping or direct probe mapping
# First try platform table
platform_info <- getGEO("GPL6480", destdir = getwd())
platform_table <- NULL

tryCatch({
    platform_table <- Table(platform_info)
    if (!is.null(platform_table) && "GENE_SYMBOL" %in% colnames(platform_table) && "ID" %in% colnames(platform_table)) {
        # Find genes in platform table (case insensitive)
        found_rows <- grep(paste(genes_of_interest, collapse = "|"), platform_table$GENE_SYMBOL, ignore.case = TRUE)
        
        if (length(found_rows) > 0) {
            # Extract matching rows
            probes_for_genes <- platform_table[found_rows, c("ID", "GENE_SYMBOL")]
            
            # Extract valid probe IDs
            valid_probes <- probes_for_genes$ID
            
            # Find matching probes in expression data
            matching_expr_probes <- intersect(valid_probes, rownames(expr_data))
            
            if (length(matching_expr_probes) > 0) {
                # Extract expression data for these probes
                gene_expr <- expr_data[matching_expr_probes, ]
                
                # Get gene symbols for these probes
                probe_to_gene <- platform_table$GENE_SYMBOL[match(rownames(gene_expr), platform_table$ID)]
                
                # Rename rows with gene symbols
                rownames(gene_expr) <- probe_to_gene
            }
        }
    }
}, error = function(e) {
    print(paste("Error with platform table:", e$message))
})

# If platform mapping failed, use direct probe mapping
if (!exists("gene_expr") || nrow(gene_expr) == 0) {
    # Direct probe mapping for the genes
    direct_probes <- list(
        "KIF20A" = "A_23_P256956",
        "PYCR1" = "A_23_P130194",
        "AURKA" = "A_23_P131866", 
        "UHRF1" = "A_23_P208880"
    )
    
    # Extract expression data for direct probes
    gene_expr <- matrix(nrow = 0, ncol = ncol(expr_data))
    for (gene in names(direct_probes)) {
        probe <- direct_probes[[gene]]
        if (probe %in% rownames(expr_data)) {
            gene_expr <- rbind(gene_expr, expr_data[probe, , drop = FALSE])
            rownames(gene_expr)[nrow(gene_expr)] <- gene
        }
    }
}

# Transpose and merge gene expression with survival data
gene_df <- t(gene_expr)
gene_df <- as.data.frame(gene_df)
surv_data_final <- cbind(surv_data_clean, gene_df)

# Step 4: Get univariate and multivariate Cox regression results
# Create a list to store univariate results
univariate_results_df <- data.frame()

# Perform univariate Cox regression for each gene
for (gene in genes_of_interest) {
    if (gene %in% colnames(surv_data_final)) {
        formula_str <- as.formula(paste("Surv(time, status) ~", gene))
        model <- coxph(formula_str, data = surv_data_final)
        sum_model <- summary(model)
        
        # Extract key statistics
        temp_df <- data.frame(
            Variable = gene,
            Analysis = "Univariate",
            HR = sum_model$conf.int[1, "exp(coef)"],
            Lower_CI = sum_model$conf.int[1, "lower .95"],
            Upper_CI = sum_model$conf.int[1, "upper .95"],
            p_value = sum_model$coefficients[1, "Pr(>|z|)"],
            stringsAsFactors = FALSE
        )
        univariate_results_df <- rbind(univariate_results_df, temp_df)
    }
}

# Create univariate models for clinical variables
clinical_vars <- c("age", "smoking")
for (var in clinical_vars) {
    formula_str <- as.formula(paste("Surv(time, status) ~", var))
    model <- coxph(formula_str, data = surv_data_final)
    sum_model <- summary(model)
    
    # Extract key statistics
    temp_df <- data.frame(
        Variable = var,
        Analysis = "Univariate",
        HR = sum_model$conf.int[1, "exp(coef)"],
        Lower_CI = sum_model$conf.int[1, "lower .95"],
        Upper_CI = sum_model$conf.int[1, "upper .95"],
        p_value = sum_model$coefficients[1, "Pr(>|z|)"],
        stringsAsFactors = FALSE
    )
    univariate_results_df <- rbind(univariate_results_df, temp_df)
}

# Create univariate model for stage
formula_str <- as.formula(paste("Surv(time, status) ~ stage"))
model <- coxph(formula_str, data = surv_data_final)
sum_model <- summary(model)

# Extract key statistics for each stage level
for (i in 1:nrow(sum_model$coefficients)) {
    var_name <- rownames(sum_model$coefficients)[i]
    temp_df <- data.frame(
        Variable = var_name,
        Analysis = "Univariate",
        HR = sum_model$conf.int[i, "exp(coef)"],
        Lower_CI = sum_model$conf.int[i, "lower .95"],
        Upper_CI = sum_model$conf.int[i, "upper .95"],
        p_value = sum_model$coefficients[i, "Pr(>|z|)"],
        stringsAsFactors = FALSE
    )
    univariate_results_df <- rbind(univariate_results_df, temp_df)
}

# Perform multivariate Cox regression
covar_list <- c(genes_of_interest, "age", "stage", "smoking")
covar_list <- covar_list[covar_list %in% colnames(surv_data_final)]
formula_str_mv <- as.formula(paste("Surv(time, status) ~", paste(covar_list, collapse = " + ")))
multivariate_model <- coxph(formula_str_mv, data = surv_data_final)
multi_results <- summary(multivariate_model)

# Extract multivariate results
multi_results_df <- data.frame(
    Variable = rownames(multi_results$coefficients),
    Analysis = "Multivariate",
    HR = multi_results$conf.int[, "exp(coef)"],
    Lower_CI = multi_results$conf.int[, "lower .95"],
    Upper_CI = multi_results$conf.int[, "upper .95"],
    p_value = multi_results$coefficients[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
)

# Step 5: Combine univariate and multivariate results
combined_results <- rbind(univariate_results_df, multi_results_df)

# Improve variable naming
combined_results$Variable <- gsub("stage", "Stage: ", combined_results$Variable)

# Organize variables into categories
gene_vars <- genes_of_interest
stage_vars <- grep("Stage:", combined_results$Variable, value = TRUE)
other_vars <- setdiff(unique(combined_results$Variable), c(gene_vars, stage_vars))

# Create ordering variable for better organization
combined_results$Category <- "Other"
combined_results$Category[combined_results$Variable %in% gene_vars] <- "Gene"
combined_results$Category[combined_results$Variable %in% stage_vars] <- "Stage"

# Sort by category and hazard ratio within category
combined_results$Category <- factor(combined_results$Category, levels = c("Gene", "Stage", "Other"))
combined_results <- combined_results[order(combined_results$Category, 
                                           combined_results$Variable, 
                                           combined_results$Analysis), ]

# Set factor levels for consistent ordering
combined_results$Variable <- factor(combined_results$Variable, 
                                    levels = unique(combined_results$Variable[order(combined_results$Category)]))

# Add significance indicator
combined_results$Significant <- ifelse(combined_results$p_value < 0.05, "Yes", "No")

# Print the combined results table
print("Combined Cox regression results:")
# Modified code for the forest plot with improved text positioning

# First, let's examine HR values to assess where we need more space
print(combined_results %>% 
          dplyr::arrange(desc(Upper_CI)) %>% 
          dplyr::select(Variable, Analysis, HR, Upper_CI, p_display))

# Fix for the P-value overlap issue
forest_plot_vertical <- ggplot(combined_results, 
                               aes(x = HR, y = Variable, xmin = Lower_CI, xmax = Upper_CI, 
                                   color = Significant, shape = Analysis)) +
    # Reference line at HR=1
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", size = 0.5) +
    
    # Points and confidence intervals with increased size
    geom_pointrange(position = position_dodge(width = 0.8), size = 1) +  # INCREASED dodge width
    
    # Set up color and shape scales
    scale_color_manual(values = c("Yes" = "#E41A1C", "No" = "#377EB8"),
                       name = "P < 0.05") +
    scale_shape_manual(values = c("Univariate" = 16, "Multivariate" = 17),
                       name = "Analysis") +
    
    # Set up log scale with MORE SPACE on right side
    scale_x_log10(breaks = c(0.2, 0.5, 1, 2, 5, 10, 20, 40), 
                  labels = c(0.2, 0.5, 1, 2, 5, 10, 20, 40),
                  limits = c(0.15, 60)) +  # Extended upper limit for text
    
    # Add labels and titles
    labs(title = "Cox Regression Analysis of Survival in LUAD",
         subtitle = "GSE13213 Dataset",
         x = "Hazard Ratio (log10 scale)",
         y = "",
         caption = "Red: P < 0.05, Blue: P ≥ 0.05 | Circles: Univariate, Triangles: Multivariate") +
    
    # Apply theme with larger font size
    theme_minimal(base_size = 14) +
    theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "gray90", size = 0.3),
        
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        # Fix for vectorized input warning
        axis.text.y = element_text(face = "plain"),
        
        plot.title = element_text(size = 20, face = "bold", margin = margin(0, 0, 15, 0)),
        plot.subtitle = element_text(size = 16, margin = margin(0, 0, 20, 0)),
        plot.caption = element_text(size = 10, hjust = 0, margin = margin(15, 0, 0, 0)),
        
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.margin = margin(15, 0, 5, 0),
        legend.box.spacing = unit(1, "cm"),
        
        strip.text = element_text(size = 16, face = "bold"),
        strip.background = element_rect(fill = "gray95", color = NA),
        panel.spacing = unit(2, "lines"),
        
        # INCREASED right margin for p-values
        plot.margin = margin(20, 60, 20, 20)
    ) +
    
    # IMPROVED p-value positioning with better offset calculation
    geom_text(aes(
        label = p_display,
        # More intelligent positioning based on HR magnitude
        x = case_when(
            Upper_CI > 30 ~ 40,                        # Fixed position for extremely large values
            Upper_CI > 10 ~ Upper_CI * 1.3,            # More space for large values
            HR > 1 ~ Upper_CI * 1.25,                  # Standard right positioning
            TRUE ~ Lower_CI * 0.7                      # Standard left positioning
        )),
        position = position_dodge(width = 0.8),        # Match dodge width with points
        hjust = ifelse(combined_results$HR > 1, 0, 1), # Alignment
        size = 4.2,                                    # Slightly reduced size to prevent overlap
        vjust = 0.5
    ) +
    
    # Make gene names bold
    geom_text(
        data = subset(combined_results, Variable %in% genes_of_interest),
        aes(x = 0.18, y = Variable, label = ""),       # Empty label as a hack to apply formatting
        hjust = 1, fontface = "bold", size = 4.5,
        inherit.aes = FALSE
    ) +
    
    # Facet by category
    facet_grid(Category ~ ., scales = "free_y", space = "free",
               labeller = labeller(Category = c(Gene = "Gene Expression", 
                                                Stage = "Disease Stage", 
                                                Other = "Clinical Factors")))

# Save with adjusted dimensions to accommodate the wider right margin
ggsave("vertical_forest_plot_fixed.pdf", forest_plot_vertical, 
       width = 12, height = 12, units = "in", dpi = 600)  # INCREASED width
ggsave("vertical_forest_plot_fixed.png", forest_plot_vertical, 
       width = 12, height = 12, units = "in", dpi = 300)
# Step 7: Save plot in different formats with LARGER dimensions
ggsave("vertical_forest_plot.pdf", forest_plot_vertical, 
       width = 10, height = 12, units = "in", dpi = 600)
ggsave("vertical_forest_plot.png", forest_plot_vertical, 
       width = 10, height = 12, units = "in", dpi = 300)
ggsave("vertical_forest_plot.tiff", forest_plot_vertical,
       width = 10, height = 12, units = "in", dpi = 300, compression = "lzw")

# Display the plot
print(forest_plot_vertical)
