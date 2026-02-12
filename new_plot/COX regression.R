# Clear workspace and load libraries
rm(list = ls())
library(forestplot)
library(dplyr)
library(readr)
library(ggplot2)

# Read the Cox regression CSV files
kif20a_cox <- read_csv("/Users/yangpaul/Downloads/KIF20A Cox Regression Results.csv", show_col_types = FALSE)
uhrf1_cox <- read_csv("/Users/yangpaul/Downloads/UHRF1 Cox Regression Results (2).csv", show_col_types = FALSE)
aurka_cox <- read_csv("/Users/yangpaul/Downloads/ AURKA Cox Regression Results (1).csv", show_col_types = FALSE)

# Fixed function to extract gene-specific rows
extract_gene_values <- function(df, gene_name) {
    # Find the gene header row
    gene_header_idx <- which(grepl(gene_name, df$Characteristics, ignore.case = TRUE))
    
    if(length(gene_header_idx) > 0) {
        # Get the next few rows after the header
        start_idx <- gene_header_idx[1]
        if(start_idx < nrow(df)) {
            subsequent_rows <- df[(start_idx + 1):min(start_idx + 3, nrow(df)), ]
            
            # Look for "Low" or "High" in characteristics (non-reference rows)
            low_row <- subsequent_rows[grepl("Low", subsequent_rows$Characteristics, ignore.case = TRUE), ]
            high_row <- subsequent_rows[grepl("High", subsequent_rows$Characteristics, ignore.case = TRUE), ]
            
            # Return the non-reference row (check for "Reference" properly)
            if(nrow(low_row) > 0) {
                mv_val <- low_row$`HR(95% CI) Multivariate analysis`[1]
                if(!is.na(mv_val) && !grepl("Reference", mv_val)) {
                    return(low_row[1, ])
                }
            } 
            
            if(nrow(high_row) > 0) {
                mv_val <- high_row$`HR(95% CI) Multivariate analysis`[1]
                if(!is.na(mv_val) && !grepl("Reference", mv_val)) {
                    return(high_row[1, ])
                }
            }
        }
    }
    return(NULL)
}

# Extract gene data
kif20a_values <- extract_gene_values(kif20a_cox, "KIF20A")
uhrf1_values <- extract_gene_values(uhrf1_cox, "UHRF1")
aurka_values <- extract_gene_values(aurka_cox, "AURKA")

# Check what we extracted
print("=== Extracted Gene Data ===")
print("KIF20A:")
print(kif20a_values)
print("\nUHRF1:")
print(uhrf1_values)
print("\nAURKA:")
print(aurka_values)

# Function to parse HR and CI from multivariate analysis
parse_hr_ci <- function(hr_string) {
    if(is.na(hr_string) || hr_string == "" || hr_string == "NA") {
        return(list(hr = NA, lower = NA, upper = NA))
    }
    
    # Extract decimal numbers
    numbers <- as.numeric(unlist(regmatches(hr_string, gregexpr("\\d+\\.\\d+", hr_string))))
    
    if(length(numbers) >= 3) {
        return(list(hr = numbers[1], lower = numbers[2], upper = numbers[3]))
    } else {
        return(list(hr = NA, lower = NA, upper = NA))
    }
}

# Parse multivariate data (only if we have valid data)
if(!is.null(kif20a_values)) {
    kif20a_mv <- parse_hr_ci(kif20a_values$`HR(95% CI) Multivariate analysis`)
} else {
    kif20a_mv <- list(hr = NA, lower = NA, upper = NA)
}

if(!is.null(uhrf1_values)) {
    uhrf1_mv <- parse_hr_ci(uhrf1_values$`HR(95% CI) Multivariate analysis`)
} else {
    uhrf1_mv <- list(hr = NA, lower = NA, upper = NA)
}

if(!is.null(aurka_values)) {
    aurka_mv <- parse_hr_ci(aurka_values$`HR(95% CI) Multivariate analysis`)
} else {
    aurka_mv <- list(hr = NA, lower = NA, upper = NA)
}

print("\n=== Parsed Multivariate Results ===")
print(paste("KIF20A:", kif20a_mv$hr, "(", kif20a_mv$lower, "-", kif20a_mv$upper, ")"))
print(paste("UHRF1:", uhrf1_mv$hr, "(", uhrf1_mv$lower, "-", uhrf1_mv$upper, ")"))
print(paste("AURKA:", aurka_mv$hr, "(", aurka_mv$lower, "-", aurka_mv$upper, ")"))

# Create forest plot data
forest_data <- data.frame(
    Gene = c("KIF20A (Low vs High)", "UHRF1 (High vs Low)", "AURKA (High vs Low)"),
    HR = c(kif20a_mv$hr, uhrf1_mv$hr, aurka_mv$hr),
    Lower = c(kif20a_mv$lower, uhrf1_mv$lower, aurka_mv$lower),
    Upper = c(kif20a_mv$upper, uhrf1_mv$upper, aurka_mv$upper),
    P_value_char = c(
        ifelse(is.null(kif20a_values), NA, kif20a_values$`P value Multivariate analysis`),
        ifelse(is.null(uhrf1_values), NA, uhrf1_values$`P value Multivariate analysis`),
        ifelse(is.null(aurka_values), NA, aurka_values$`P value Multivariate analysis`)
    ),
    stringsAsFactors = FALSE
)

print("\n=== Forest Plot Data ===")
print(forest_data)

# Convert p-values to numeric
convert_p_value <- function(p_str) {
    if(is.na(p_str) || p_str == "" || p_str == "NA") return(NA)
    if(grepl("< 0.001", p_str)) return(0.0005)
    if(grepl("^< ", p_str)) return(as.numeric(gsub("< ", "", p_str)))
    return(as.numeric(p_str))
}

forest_data$P_numeric <- sapply(forest_data$P_value_char, convert_p_value)
forest_data$Significant <- forest_data$P_numeric < 0.05

# Filter out rows with missing data
forest_data_clean <- forest_data %>%
    filter(!is.na(HR) & !is.na(Lower) & !is.na(Upper))

print("\n=== Clean Forest Plot Data ===")
print(forest_data_clean)

# Only create plot if we have data
# Create the corrected forest plot (fixing ggplot2 deprecation issues)
if(nrow(forest_data_clean) > 0) {
    p_forest <- ggplot(forest_data_clean, aes(x = reorder(Gene, HR), y = HR)) +
        geom_point(aes(color = Significant), size = 5) +
        geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Significant), 
                      width = 0.2, linewidth = 1.2) +  # Changed size to linewidth
        geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.8, linewidth = 1) +
        scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "darkred"),
                           labels = c("FALSE" = "Not Significant (p ≥ 0.05)", 
                                      "TRUE" = "Significant (p < 0.05)")) +
        coord_flip() +
        scale_y_continuous(trans = "log10", 
                           breaks = c(0.3, 0.5, 0.7, 1.0, 1.5, 2.0),
                           labels = c("0.3", "0.5", "0.7", "1.0", "1.5", "2.0")) +
        theme_minimal() +
        theme(
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.title = element_text(size = 11, face = "bold"),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_line(color = "grey90"),  # Fixed alpha issue
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)  # Fixed size to linewidth
        ) +
        labs(
            title = "Gene Expression and Overall Survival",
            subtitle = "Multivariate Cox Proportional Hazards Regression",
            x = "Gene Expression Comparison",
            y = "Hazard Ratio (95% Confidence Interval)",
            color = "Statistical Significance",
            caption = "Red dashed line indicates HR = 1.0 (no effect)\nError bars represent 95% confidence intervals"
        )
    
    print(p_forest)
    
    # Create an annotated version with HR and p-values
    p_forest_annotated <- p_forest +
        geom_text(aes(label = paste0("HR = ", round(HR, 3), "\np = ", P_value_char)), 
                  nudge_x = 0.3, size = 3.5, hjust = 0, 
                  color = "black", fontface = "bold")
    
    print("\n=== Annotated Forest Plot ===")
    print(p_forest_annotated)
    
} else {
    print("No valid data found for forest plot")
}

# Create a summary table
library(kableExtra)

summary_table <- forest_data_clean %>%
    mutate(
        `HR (95% CI)` = paste0(round(HR, 3), " (", round(Lower, 3), "-", round(Upper, 3), ")"),
        `P-value` = P_value_char,
        `Interpretation` = case_when(
            Significant ~ "Significant prognostic factor",
            !Significant ~ "Not significant"
        )
    ) %>%
    select(Gene, `HR (95% CI)`, `P-value`, Interpretation)

print("\n=== Summary Table ===")
print(summary_table)

# Create a publication-ready table
kable(summary_table, 
      caption = "Multivariate Cox Regression: Gene Expression and Overall Survival",
      align = c("l", "c", "c", "l")) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = FALSE) %>%
    row_spec(which(forest_data_clean$Significant), 
             background = "#ffebee", bold = TRUE) %>%
    footnote(general = "Significant results (p < 0.05) are highlighted in red")

