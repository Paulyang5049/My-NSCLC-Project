# Clear workspace
rm(list = ls())

library(GEOquery)
gset <- getGEO('GSE118370', destdir=".", getGPL = T)
class(gset)
gset[[1]]

pdata <- pData(gset[[1]])
table(pdata$source_name_ch1)

library(stringr)
group_list <- ifelse(str_detect(pdata$source_name_ch1, "normal lung tissue"), "Control", 'Tumor')
group_list = factor (group_list,
                     levels = c("Control", "Tumor"))
table (group_list)


exp <-  exprs(gset[[1]])
boxplot(exp,outline=FALSE, notch=T, col=group_list, las=2)
library (limma)
exp=normalizeBetweenArrays(exp)
boxplot(exp,outline=FALSE, notch=T, col=group_list, las=2)
range(exp)
# Log-transform and visualize data
#exp <- log2(exp + 1) # Uncomment if log transformation is needed
boxplot(exp, col = "skyblue", main = "Boxplot of Raw Data")

# Save initial data
save_data(exp, clinical, "step1.RData")

# Reload data and process group information
rm(list = ls())
load_data("step1.RData")

# 6. Extract Group Information
extract_group_info <- function(clinical) {
  Group <- ifelse(str_detect(clinical$source_name_ch1, "Adjacent"), "Control", "Tumor")
  Group <- factor(Group, levels = c("Control", "Tumor"))
  return(Group)
}

# Extract group information
Group <- extract_group_info(clinical)

# 7. Probe ID Conversion
convert_probe_id <- function(exp, file_path) {
  id <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, comment.char = "#")
  id <- id[, c(1, 6, 10, 13)]
  exp <- exp[match(id$ID, rownames(exp)), ]
  rownames(exp) <- id$ILMN_Gene
  return(exp)
}

# Convert probe IDs
exp <- convert_probe_id(exp, "GPL6884-11607.txt")

# Save data after ID conversion
save(exp, clinical, id, Group, file = "step2.RData")

# Reload data for PCA and heatmap
rm(list = ls())
load("step2.RData")
library(scatterplot3d)
library(tidyverse)
library(openxlsx)
library(factoextra) # Extract and Visualize the Results of Multivariate Data Analyses
color.bin <- c("#00599F","#D01910")
res.pca.comp <- prcomp(exp, scale = F)
plot.data <- as.data.frame(res.pca.comp$rotation[, 1:10])
plot.data <- plot.data %>% 
  mutate(ID=rownames(plot.data),
         Type=Group,
         TypeColor=color.bin[as.numeric(as.factor(Group))])
pdf("./Plot/01. GEO/GSE32863/PCA3d.pdf", width = 7, height = 7)
scatterplot3d(x = plot.data$PC2, 
              y = plot.data$PC1, 
              z = plot.data$PC3,
              color = plot.data$TypeColor,
              pch = 16, cex.symbols = 1,
              scale.y = 0.7, angle = 45,
              xlab = "PC2", ylab = "PC1", zlab = "PC3",
              main="3D Scatter Plot of PCA",
              col.axis = "#444444", col.grid = "#CCCCCC")
legend("bottom", legend = levels(as.factor(Group)),
       col =  color.bin,  pch = 16,
       inset = -0.25, xpd = TRUE, horiz = TRUE)
dev.off()
write.xlsx(plot.data, "./Plot/01. GEO/GSE32863/PCA3d.xlsx", overwrite = T)
# variance percent of each PC
p <- fviz_eig(res.pca.comp)
var_explained <- get_eig(res.pca.comp)
# var_explained <- res.pca.comp$sdev^2 / sum(res.pca.comp$sdev^2)
ggsave("./Plot/01. GEO/GSE32863/PCA_percent.pdf",p,width = 5, height = 5)
write.xlsx(var_explained, "./Plot/01. GEO/GSE32863/PCA_percent.xlsx",rowNames=T, overwrite = T)


# # 8. Perform PCA
# perform_pca <- function(exp, Group) {
#   dat <- as.data.frame(t(exp))
#   dat.pca <- PCA(dat, graph = FALSE)
#   pca_plot <- fviz_pca_ind(dat.pca,
#                            geom.ind = "point",
#                            col.ind = Group,
#                            palette = c("#00AFBB", "#E7B800"),
#                            addEllipses = TRUE,
#                            legend.title = "Groups")
#   return(pca_plot)
# }
# 
# # Generate and save PCA plot
# pca_plot <- perform_pca(exp, Group)
# print(pca_plot)
# save(pca_plot, file = "GSE32863_pca_plot.Rdata")

# 9. Generate Heatmap
# Load necessary libraries
library(pheatmap)

# Generate Heatmap with Group Clustering
generate_heatmap <- function(exp, Group) {
  # Define custom color palette
  custom_palette <- colorRampPalette(c("#599CB4","#92B5CA","#AECFD4","#CCE4EF",
                                       "#F5DFDB","#EDB8B0","#E69191","#C25759"))(100)
  
  # Calculate standard deviation for each row and select top 1000 genes
  cg <- names(tail(sort(apply(exp, 1, sd)), 1000))
  n <- exp[cg, ]
  
  # Order columns by group to cluster Control and Tumor separately
  ordered_indices <- order(Group)
  n <- n[, ordered_indices]
  Group <- Group[ordered_indices]
  
  # Create annotation data frame
  annotation_col <- data.frame(group = Group)
  rownames(annotation_col) <- colnames(n)
  
  
  # Generate the heatmap
  heatmap <- pheatmap(n,
                      color = custom_palette,
                      show_colnames = FALSE,
                      show_rownames = FALSE,
                      annotation_col = annotation_col,
                      scale = "row",
                      breaks = seq(-3, 3, length.out = 100),
                      fontsize = 10, # Font size for readability
                      fontsize_row = 6,
                      fontsize_col = 6,
                      main = "GSE32863 Heatmap", # Title
                      clustering_distance_rows = "euclidean",
                      clustering_method = "complete", # Hierarchical clustering
                      cluster_cols = FALSE) # Do not cluster columns
  
  return(heatmap)
}

# Generate and save the heatmap
heatmap <- generate_heatmap(exp, Group)
print(heatmap)
dev.off()

# Reload data for differential expression analysis
rm(list = ls())
load_data("step2.RData")

# 10. Differential Expression Analysis
perform_DE_analysis <- function(exp, Group) {
  design <- model.matrix(~Group)
  fit <- lmFit(exp, design)
  fit <- eBayes(fit)
  deg <- topTable(fit, coef = 2, number = Inf)
  return(deg)
}

# Perform DE analysis
deg <- perform_DE_analysis(exp, Group)

# 11. Add Change Column
add_change_column <- function(deg, logFC_t = 1, P.Value_t = 0.05) {
  k1 <- (deg$P.Value < P.Value_t) & (deg$logFC < -logFC_t)
  k2 <- (deg$P.Value < P.Value_t) & (deg$logFC > logFC_t)
  deg <- mutate(deg, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
  return(deg)
}

# Update DE results with change column
deg <- add_change_column(deg)

# Make the 'ID' column unique and set as row names
deg$ID <- make.unique(as.character(deg$ID))
rownames(deg) <- deg$ID
deg <- deg[,-which(names(deg) == "ID")]

# Save DE results
save(deg, file = "GSE32863_deg.Rdata")

# 12. Draw Volcano Plot
draw_volcano <- function(data, logFC_cutoff = 1, p_value_cutoff = 0.05, title = "Volcano Plot") {
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
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    labs(x = "log2 Fold Change", y = "-log10 P Value") +
    ggtitle(title) +
    geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(p_value_cutoff), linetype = "dashed", color = "black") +
    theme(
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "gray90")
    )
  
  return(volcano)
}

# Load DE results and plot volcano
volcano_plot <- draw_volcano(deg, logFC_cutoff = 1, p_value_cutoff = 0.05, title = "GSE32863 Volcano Plot")
print(volcano_plot)
