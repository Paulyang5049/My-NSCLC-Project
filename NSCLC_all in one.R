# Single-cell RNA-seq - normalization
load('~/R/seurat_filtered.RData')
# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)
# Load cell cycle markers
load("~/R/cycle.rda")

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase,
                                 g2m.features = g2m_genes,
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells
View(seurat_phase@meta.data)

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase,
                                     selection.method = "vst",
                                     nfeatures = 2000,
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Identify the 15 most highly variable genes
ranked_variable_genes <- VariableFeatures(seurat_phase)
top_genes <- ranked_variable_genes[1:15]

# Plot the average expression and variance of these genes
# With labels to indicate which genes are in the top 15
p <- VariableFeaturePlot(seurat_phase)
LabelPoints(plot = p, points = top_genes, repel = TRUE)


# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")

# Check quartile values
summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio,
                                     breaks=c(-Inf, 0.0144, 0.0199, 0.0267, Inf),
                                     labels=c("Low","Medium","Medium high", "High"))

## DO NOT RUN CODE ##

# SCTranform
seurat_phase <- SCTransform(seurat_phase, vars.to.regress = c("mitoRatio"))

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "sample")


options(future.globals.maxSize = 4000 * 1024^2)


for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"), vst.flavor = "v2")
}
# Check which assays are stored in objects
split_seurat$IAC@assays

# Save the split seurat object
saveRDS(split_seurat, "~/R/data/split_seurat.rds")
saveRDS(seurat_phase, "~/R/data/seurat_phase.rds")










rm(list=ls())
# Load the split seurat object into the environment
split_seurat <- readRDS("~/R/data/split_seurat.rds")
seurat_phase <- readRDS("~/R/data/seurat_phase.rds")
seurat_filtered <- load("~/R/data/seurat_filtered.RData")
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(ggplot2)
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

seurat_phase <- RunUMAP(seurat_phase,
                        dims = 1:30,reduction = "pca")
DimPlot(seurat_phase)



# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat,
                                            nfeatures = 2000)


# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat,
                                   anchor.features = integ_features)


# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat,
                                        normalization.method = "SCT",
                                        anchor.features = integ_features)


# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors,
                                   normalization.method = "SCT")

###UMAP visualization
# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

# Plot PCA
PCAPlot(seurat_integrated)
PCAPlot(seurat_integrated,
        split.by = "sample")

# Set seed
set.seed(123456)

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated,
                             dims = 1:30,
                             reduction = "pca")

# Plot UMAP
DimPlot(seurat_integrated)

# Plot UMAP split by sample
DimPlot(seurat_integrated,
        split.by = "sample")


# Save integrated seurat object
saveRDS(seurat_integrated, "~/R/results/integrated_seurat.rds")







# Single-cell RNA-seq - clustering
rm(list=ls())
# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(ggplot2)
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

seurat_integrated <- readRDS("~/R/results/integrated_seurat.rds")
# Explore heatmap of PCs
DimHeatmap(seurat_integrated,
           dims = 1:9,
           cells = 500,
           balanced = TRUE)


# Printing out the most variable genes driving PCs
print(x = seurat_integrated[["pca"]],
      dims = 1:10,
      nfeatures = 5)

# Plot the elbow plot
ElbowPlot(object = seurat_integrated,
          ndims = 40)

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated,
                                   dims = 1:20,
                                   reduction = 'pca')
# Determine the clusters for various resolutions
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))


# Explore resolutions
seurat_integrated@meta.data %>%
  View()


# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

## Calculation of UMAP
## DO NOT RUN (calculated in the last lesson)

seurat_integrated <- RunUMAP(seurat_integrated,
                             reduction = "pca",
                             dims = 1:30)

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)


# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)



## QC of cluster analysis
# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_integrated,
                     vars = c("ident", "sample")) %>%
  dplyr::count(ident, sample)

# Barplot of number of cells per cluster by sample
ggplot(n_cells, aes(x=ident, y=n, fill=sample)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_text(aes(label=n), vjust = -.2, position=position_dodge(1))

# UMAP of cells in each cluster by sample
DimPlot(seurat_integrated,
        label = TRUE,
        split.by = "sample")  + NoLegend()

# Barplot of proportion of cells in each cluster by sample
ggplot(seurat_integrated@meta.data) +
  geom_bar(aes(x=integrated_snn_res.0.8, fill=sample), position=position_fill())

# Explore whether clusters segregate by cell cycle phase
DimPlot(seurat_integrated,
        label = TRUE,
        split.by = "Phase")  + NoLegend()

# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

FeaturePlot(seurat_integrated,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)


# Boxplot of nGene per cluster
ggplot(seurat_integrated@meta.data) +
  geom_boxplot(aes(x=integrated_snn_res.0.8, y=nGene, fill=integrated_snn_res.0.8)) +NoLegend()


# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:16),
             "ident",
             "umap_1", "umap_2")

# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated,
                     vars = columns)


# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated,
                        vars = c("ident", "umap_1", "umap_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(umap_1), y=mean(umap_2))

# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:16), function(pc){
  ggplot(pc_data,
         aes(umap_1, umap_2)) +
    geom_point(aes_string(color=pc),
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE,
                         low = "grey90",
                         high = "blue")  +
    geom_text(data=umap_label,
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>%
  plot_grid(plotlist = .)


# Examine PCA results
print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)




# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)
## CD14+ monocyte markers
FeaturePlot(seurat_integrated,
            reduction = "umap",
            features = c("CD14", "LYZ"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
## FCGR3A+ monocyte markers
FeaturePlot(seurat_integrated,
            reduction = "umap",
            features = c("FCGR3A", "MS4A7"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

## Macrophages
FeaturePlot(seurat_integrated,
            reduction = "umap",
            features = c("MARCO", "ITGAM", "ADGRE1"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

## Conventional dendritic cell markers
FeaturePlot(seurat_integrated,
            reduction = "umap",
            features = c("FCER1A", "CST3"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

## Plasmacytoid dendritic cell markers
FeaturePlot(seurat_integrated,
            reduction = "umap",
            features = c("IL3RA", "GZMB", "SERPINF1", "ITM2C"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

## Hypoxia related
FeaturePlot(seurat_integrated,
            reduction = "umap",
            features = c("HIF1A"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            split.by = "sample")


# List of known celltype markers
markers <- list()
markers[["CD14+ monocytes"]] <- c("CD14", "LYZ")
markers[["FCGR3A+ monocyte"]] <- c("FCGR3A", "MS4A7")
markers[["Macrophages"]] <- c("MARCO", "ITGAM", "ADGRE1")
markers[["Conventional dendritic"]] <- c("FCER1A", "CST3")
markers[["Plasmacytoid dendritic"]] <- c("IL3RA", "GZMB", "SERPINF1", "ITM2C")

# Create dotplot based on RNA expression
DotPlot(seurat_integrated, markers, assay="RNA")


# save the file
saveRDS(seurat_integrated, '~/R/data/seurat_integrated.rds')


rm(list=ls())

#加载R包
library(Seurat)
library(SeuratObject)
library(fastSave)
library(celldex)
library(assertthat)
library(monocle)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(tricycle)  
library(scattermore)
library(scater)
library(patchwork)
library(ggplot2)
library(CCA)
library(clustree)
library(cowplot)
library(SCpubr)
library(UCell)
library(irGSEA)
library(GSVA)
library(GSEABase)
library(harmony)
library(plyr)
library(randomcoloR)
library(CellChat)
library(future)
library(ggforce)
library(ggsci)
library(parallel)
library(doParallel)
library(data.table)
library(qs)
library(presto)
library(COSG)
library(starTracer)

seurat_integrated <- readRDS('~/R/GSE189357/data/seurat_integrated.rds')
annotations <- read.csv("~/R/GSE189357/data/annotation.csv")
seurat_integrated <- JoinLayers(seurat_integrated)


#多种算法计算 marker genes
# FindAllMarkers
# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(object = seurat_integrated,
                          only.pos = TRUE,
                          logfc.threshold = 0.25,
                          min.pct = 0.25)
write.csv(markers, file = "1.FindAllMarkers每个细胞Maeker基因.csv")

DefaultAssay(seurat_integrated) <- "RNA"
table(Idents(seurat_integrated))


# 0    1    2    3    4    5    6    7    8 
# 1038  739  340  322  227  195  173  141  131 
# 9   10   11   12 
# 61   58   45   39 

top5scRNA.markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

col <- c(ggsci::pal_npg()(9),
         ggsci::pal_jco()(9),
         ggsci::pal_jama()(7),
         ggsci::pal_nejm()(8))

pdf(file = "1.FindAllMarkers每个细胞前5基因热图.pdf",width =22,height = 16)

DoHeatmap(seurat_integrated, features = top5scRNA.markers$gene,
          group.colors = col) +
  ggsci::scale_colour_npg() +
  scale_fill_gradient2(low = '#0099CC',
                       mid = 'white',
                       high = '#CC0033',
                       name = 'Z-score')
dev.off()


cluster0_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                   ident.1 = 0,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 1)




# Combine markers with gene descriptions
cluster0_ann_markers <- cluster0_conserved_markers %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

View(cluster0_ann_markers)




# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}
# Iterate function across desired clusters
# conserved_markers <- map_dfr(c(4,0,6,2), get_conserved)

# Get all unique cluster identities
all_clusters <- unique(Idents(seurat_integrated))

# Iterate function across all clusters
conserved_markers <- map_dfr(all_clusters, get_conserved)


# Extract top 10 markers per cluster
top5 <- conserved_markers %>%
  mutate(avg_fc = (`AIS_avg_log2FC` + `MIA_avg_log2FC` +`IAC_avg_log2FC`) /3) %>%
  group_by(cluster_id) %>%
  top_n(n = 5,
        wt = avg_fc)
# Visualize top 5 markers per cluster
View(top5)


# Define a function to find markers for a specific cluster in the caries group
get_caries_markers <- function(cluster_id) {
  # Use FindMarkers on the caries-only subset for the specified cluster
  markers <- FindMarkers(caries_only, ident.1 = cluster_id, only.pos = TRUE, logfc.threshold = 1) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster_id, .)
  
  return(markers)
}



# write.csv(top10, file = "results/top10_markers_per_cluster.csv")



# # Plot interesting marker gene expression for cluster 0
# FeaturePlot(object = seurat_integrated,
#             features = c("PCK2", "GAPDHS", "HK1"),
#             order = TRUE,
#             min.cutoff = 'q10',
#             label = TRUE,
#             repel = TRUE,
#             split.by = )
# 
# # Vln plot - cluster 4
# VlnPlot(object = seurat_integrated,
#         features = c("PCK2", "GAPDHS", "HK1"))


## Identifying gene markers for each cluster

# # Determine differentiating markers for CD4+ T cell
# cd4_tcells <- FindMarkers(seurat_integrated,
#                           ident.1 = 2,
#                           ident.2 = c(0,4,6))
# 
# # Add gene symbols to the DE table
# cd4_tcells <- cd4_tcells %>%
#   rownames_to_column(var = "gene") %>%
#   left_join(y = unique(annotations[, c("gene_name", "description")]),
#             by = c("gene" = "gene_name"))
# 
# # Reorder columns and sort by padj
# cd4_tcells <- cd4_tcells[, c(1, 3:5,2,6:7)]
# 
# cd4_tcells <- cd4_tcells %>%
#   dplyr::arrange(p_val_adj)
# 
# # View data
# View(cd4_tcells)

# Rename all identities
seurat_integrated_labelled <- RenameIdents(object = seurat_integrated,
                                           "0" = "Tcm",
                                           "1" = "CD8+ T Cells",
                                           "2" = "B Cells",
                                           "3" = "NK Cells",
                                           "4" = "cDC2",
                                           "5" = "Club Cells",
                                           "6" = "Mast Cells",
                                           "7" = "NK Cells",
                                           "8" = "Immature Neutrophil",
                                           "9" = "Plasma Cells",
                                           "10" = "Macrophages",
                                           "11" = "Macrophages",
                                           "12" = "Endothelial Cells",
                                           "13" = "Malignant AT2",
                                           "14" = "eTreg",
                                           "15" = "CAF",
                                           "16" = "Mature Neutrophil",
                                           "17" = "AT2",
                                           "18" = "B Cells",
                                           "19" = "Ciliated Cells")
# library(celldex)
# library(SingleR)
# ls("package:celldex")
# f = "/Users/haoran/R working directory/scRNA_group/scRNA_group/day06/ref_BlueprintEncode.RData"
# #用了昨天文件夹里的数据，day6是昨天的文件夹名，按需修改
# if(!file.exists(f)){
#   ref <- celldex::BlueprintEncodeData()
#   save(ref,file = f)
# }
# ref <- get(load(f))
# library(BiocParallel)
# scRNA = seurat_integrated
# test = scRNA@assays$RNA$data
# pred.scRNA <- SingleR(test = test,
#                       ref = ref,
#                       labels = ref$label.main,
#                       clusters = scRNA@active.ident)
# pred.scRNA$pruned.labels
# new.cluster.ids <- pred.scRNA$pruned.labels
# names(new.cluster.ids) <- levels(scRNA)
# seurat_integrated <- RenameIdents(scRNA,new.cluster.ids)


# Plot the UMAP
DimPlot(object = seurat_integrated_labelled,
        reduction = "umap",
        label = TRUE,
        label.size = 3,
        repel = TRUE)

library(SCP)
seurat_integrated_labelled$celltype <- Idents(seurat_integrated_labelled)
CellDimPlot(seurat_integrated_labelled,group.by="celltype",reduction="umap")
#左下小箭头
CellDimPlot(seurat_integrated_labelled,group.by="celltype",
            reduction="umap",
            theme_use="theme_blank")

library(Nebulosa)
#绘制一个基因PTPRC
p1<-plot_density(seurat_integrated_labelled,c("KIF20A"),size=0.3)
p1
# #同时绘制多个基因
# p3<-plot_density(seurat_integrated_labelled,c("CD3D","CD4","CD8A","CD68"),size=0.3)
# p3
# #combine=FALSE只保留两个基因都表达的区域
# p4<-plot_density(seurat_integrated_labelled,c("CD8A","CCR7"),
#                  joint=TRUE,
#                  combine=FALSE)
# p4


# FeaturePlot(object = seurat_integrated,
#             features = c("PCK2", "GAPDHS", "HK1"),
#             order = TRUE,
#             min.cutoff = 'q10',
#             label = TRUE,
#             repel = TRUE)


# # Remove the unwanted cells
# seurat_integrated_labelled <- subset(seurat_integrated_labelled,
#                                      idents = c('B Cells','Plasma Cells'), invert = TRUE)
#
#
# # Re-visualize the clusters
# DimPlot(object = seurat_integrated_labelled,
#         reduction = "umap",
#         label = TRUE,
#         label.size = 3,
#         repel = TRUE)





# Add celltype annotation as a column in meta.data
qsave(seurat_integrated_labelled,
      "手动注释.qs",
      nthreads = detectCores())



rm(list = ls())
gc()
#加载R包
library(Seurat)
library(SeuratObject)
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(DOSE)
library(celldex)
library(assertthat)
library(monocle)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(tricycle)   
library(scattermore)
library(scater)
library(patchwork)
library(ggplot2)
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(SCpubr)
library(UCell)
library(irGSEA)
library(GSVA)
library(GSEABase)
library(harmony)
library(plyr)
library(randomcoloR)
library(CellChat)
library(future)
library(ggplot2)
library(ggforce)
library(ggsci)
library(parallel)
library(doParallel)
library(data.table)
library(qs)
#生成随机颜色
randomColor <- function() {
  paste0("#",paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE),collapse = ""))
}

# 生成100个随机颜色
randomColors <- replicate(100,randomColor())
scRNA=qread("手动注释.qs",nthreads = detectCores())#读取数据

#挑选每个细胞亚群中高表达的前200个基因
sample.markers <- FindAllMarkers(scRNA,
                                 only.pos = TRUE,
                                 min.pct = 0.25, logfc.threshold = 0.25)
top200 <- sample.markers %>% group_by(cluster) %>% top_n(n = 200,
                                                         wt = avg_log2FC) 
top200_gene<-unstack(top200,gene~cluster)
#将gene symbol转换成entriz gene id
#做KEGG富集分析需要
top200_entrez <- lapply(X = top200_gene, FUN = function(x) {
  x <- bitr(x, fromType="SYMBOL",
            toType="ENTREZID",
            OrgDb="org.Hs.eg.db")
  x <- x[,-1]
})


######################################
#对所有亚群做GO富集分析 Gene ontoloty
#####################################
######################################
# CC cellular component
######################################
all_cc <- compareCluster(top200_entrez, 
                         fun='enrichGO',
                         ont= 'CC',
                         OrgDb='org.Hs.eg.db')
write.csv(as.data.frame(all_cc),"1.所有亚群GO分析_cc.csv")

#绘制CC富集图
pdf('2.所有亚群GO分析_CC.pdf', width = 15, height = 6)
dotplot(all_cc, includeAll=FALSE,label_format=100)+ theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

##################################################
# MF molecular function
#################################################
all_mf <- compareCluster(top200_entrez, 
                         fun='enrichGO',
                         ont= 'MF',
                         OrgDb='org.Hs.eg.db')
write.csv(as.data.frame(all_mf),"3.所有亚群GO分析_mf.csv")

#绘制MF富集图
pdf('4.所有亚群GO分析_MF.pdf', width = 15, height = 6)
dotplot(all_mf, includeAll=FALSE,label_format=100)+ theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()
#############################################
# BP biological process
#############################################
all_bp <- compareCluster(top200_entrez, 
                         fun='enrichGO',
                         ont= 'BP',
                         OrgDb='org.Hs.eg.db')
write.csv(as.data.frame(all_bp),"5.所有亚群GO分析_bp.csv")

#绘制BP富集图
pdf('6.所有亚群GO分析_BP.pdf', width = 15, height = 6)
dotplot(all_bp, includeAll=FALSE,label_format=100)+ theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()


#批量对每个细胞进行GO富集分析
lapply(1:length(top200_entrez),function(i){
  all_go_cell  <- enrichGO(gene = top200_entrez[[i]],  
                           OrgDb = 'org.Hs.eg.db',
                           ont = 'ALL',
                           pvalueCutoff = 0.05,
                           qvalueCutoff =0.2)
  
  all_go_cell<-setReadable(all_go_cell, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  write.csv(as.data.frame(all_go_cell),paste0("7.",names(top200_entrez[i]),"GO分析.csv"))
  plotGO<-dotplot(all_go_cell, split="ONTOLOGY",label_format=100) + facet_grid(ONTOLOGY~., scale="free")
  pdf(paste0('8.',names(top200_entrez[i]),'GO分析.pdf'), width = 16, height = 8)
  print(plotGO)
  dev.off()
})
###################################
#对所有亚群做KEGG富集分析
###################################
all_kegg <- compareCluster(top200_entrez,
                           fun='enrichKEGG',
                           pvalueCutoff=0.05, 
                           organism="hsa"
)

write.csv(as.data.frame(all_kegg),"9.所有亚群_KEGG.csv")
pdf('10.所有亚群_KEGG.pdf', width = 15, height = 8)
dotplot(all_kegg, showCategory=5, includeAll=FALSE,label_format=100)+ theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

#批量对每个细胞进行KEGG富集分析
lapply(1:length(top200_entrez),function(i){
  kegg_cell <- enrichKEGG(gene = top200_entrez[[i]],
                          organism  = 'hsa', 
                          pvalueCutoff = 1,
                          qvalueCutoff =1)
  kegg_cell<-setReadable(kegg_cell, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  write.csv(kegg_cell,paste0("11.",names(top200_entrez[i]),"特定亚群KEGG.csv"))
  plotKEGG<-dotplot(kegg_cell,label_format=100)
  pdf(paste0('12.',names(top200_entrez[i]),'特定亚群KEGG.pdf'), width = 12, height = 6)
  print(plotKEGG)
  dev.off()
})



rm(list = ls())
gc()
#加载R包
library(Seurat)
library(SeuratObject)
library(celldex)
library(assertthat)
library(monocle3)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(tricycle)   
library(scattermore)
library(scater)
library(patchwork)
library(ggplot2)
library(CCA)
library(clustree)
library(cowplot)
library(monocle3)
library(tidyverse)
library(SCpubr)
library(UCell)
library(irGSEA)
library(GSVA)
library(GSEABase)
library(harmony)
library(plyr)
library(randomcoloR)
library(CellChat)
library(future)
library(ggplot2)
library(ggforce)
library(ggsci)
library(parallel)
library(doParallel)
library(data.table)
library(qs)

#生成随机颜色
randomColor <- function() {
  paste0("#",paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE),collapse = ""))
}

# 生成100个随机颜色
randomColors <- replicate(100,randomColor())
scRNA=qread("手动注释.qs",nthreads = detectCores())#读取数据
###拟时序分析
##随机抽取1%个细胞，仅作为演示，实际操作或高性能工作站直接运行scRNA_tpm<-scRNA
scRNA_tpm<-scRNA
data <-GetAssayData(scRNA_tpm, assay = 'RNA', slot = 'counts')
cell_metadata <- scRNA_tpm@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data),row.names = row.names(data))
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#标准化
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
#umap降维
#umap.fast_sgd=T进行随机降维,cores为多线程运算数量,这两个参数都能加快处理速度,但是windows系统下无效
cds <- reduce_dimension(cds,reduction_method = "UMAP",cores = 1,umap.fast_sgd=F)
plot_cells(cds)
colnames(colData(cds))
#以之前的Seurat坐标来添加颜色
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
p1
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(scRNA_tpm, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')
p2
#绘制与原来UMAP的对比
p = p1|p2
p
ggsave("1.Monocle3与seurat对比.pdf", plot = p, width = 10, height = 5)
#展示指定基因
genes <- c("KIF20A")
pdf("2.特定基因表达.pdf")
plot_cells(cds,
           genes=genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
dev.off()
## 识别轨迹
cds <- cluster_cells(cds,cluster_method = "louvain")
#如果出现Error in leidenbase::leiden_find_partition(graph_result[["g"]], partition_type = partition_type,  : REAL() can only be applied to a 'numeric', not a 'NULL'报错，将cluster_method设置为“louvain“，或者将 igraph 软件包降级到 1.4.3 版本
plot_cells(cds, color_cells_by = "celltype")
cds <- learn_graph(cds)
p = plot_cells(cds, color_cells_by = "celltype",
               label_groups_by_cluster=FALSE,
               label_leaves=FALSE,
               label_branch_points=FALSE)
ggsave("3.细胞轨迹.pdf", plot = p, width = 8, height = 6)
p=plot_cells(cds,
             color_cells_by = "celltype",
             label_cell_groups=FALSE,
             label_leaves=TRUE,
             label_branch_points=TRUE,
             graph_label_size=1.5)
ggsave("4.发育节点.pdf", plot = p, width = 8, height = 6)

#选择起始点
cds <- order_cells(cds)
p = plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
               label_leaves = FALSE,  label_branch_points = FALSE)
ggsave("5.手动节点发育时间.pdf", plot = p, width = 8, height = 6)
#筛选时序差异基因
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
write.csv(Track_genes, "6.时序差异基.csv", row.names = F)
#挑选前10个基因进行绘图
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
#绘制这10个基因的趋势图
p <- plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="celltype", 
                              min_expr=0.5, ncol = 2)
ggsave("7.基因趋势图.pdf", plot = p, width = 8, height = 6)
# Plot specific gene KIF20A
p <- plot_genes_in_pseudotime(cds["KIF20A",], 
                              color_cells_by="celltype", 
                              min_expr=0.5)
ggsave("8.KIF20A_趋势图.pdf", plot = p, width = 8, height = 6)





rm(list = ls())
gc()
scRNA=qread("手动注释.qs",nthreads = detectCores())
###细胞通讯分析
#获取表达矩阵
cellchat <- createCellChat(GetAssayData(scRNA, assay = "RNA", layer = "data"))
#获取注释细胞标签
meta <- data.frame(cellType = scRNA$celltype, row.names =  Cells(scRNA))
#把metadata信息加到CellChat对象中，添加细胞标签
cellchat <- addMeta(cellchat, meta = meta, meta.name = "celltype")
#把细胞标签设置成默认的ID
cellchat <- setIdent(cellchat, ident.use = "celltype") 
#统计每个细胞亚群中的细胞数目
groupSize <- as.numeric(table(cellchat@idents)) 

#CellChat提供了人和小鼠的配受体数据库，分别可以用CellChatDB.human,CellChatDB.mouse来导入
CellChatDB <- CellChatDB.human


#择特定的信息描述细胞间的相互作用，比用一个大的配体库更精细
#用Secreted Signaling来分析细胞通信
#查看有哪些可以选择的子数据库
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")


#将使用的数据库信息写入到cellchat对象中
cellchat@DB <- CellChatDB.use

#抽取信号通路相关基因的表达矩阵
cellchat <- subsetData(cellchat) 


#对表达数据进行预处理，用于细胞间的通信分析。
#首先在一个细胞组中识别过表达的配体或受体，然后将基因表达数据投射到蛋白-蛋白相互作用(PPI)网络上。
#如果配体或受体过表达，则识别过表达配体和受体之间的相互作用。
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat)

#相互作用推断
#为每个相互作用分配一个概率值并进行置换检验来推断生物意义上的细胞-细胞通信
cellchat <- computeCommunProb(cellchat)

#通过计算与每个信号通路相关的所有配体-受体相互作用的通信概率来推断信号通路水平上的通信概率
#通过计算链路的数量或汇总通信概率来计算细胞间的聚合通信网络
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

#可视化聚合的通讯网络
groupSize <- as.numeric(table(cellchat@idents))
pdf("1.细胞通讯网络_细胞交互次数.pdf")
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F,
                 title.name = "Number of interactions")
dev.off()

pdf("2.细胞通讯网络_细胞交互权重.pdf")
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$weight,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F,
                 title.name = "Interaction weights/strength")
dev.off()

#比较不同细胞通讯网络网络之间的权重
mat <- cellchat@net$weight
pdf("3.不同细胞通讯网络权重.pdf")


par(mfrow = c(2,3), xpd=TRUE)

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat),
                 ncol = ncol(mat),
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2,
                   vertex.weight = groupSize,
                   weight.scale = T,
                   edge.weight.max = max(mat),
                   title.name = rownames(mat)[i])
}
dev.off()


#显示重要通信的信号路径
cellchat@netP$pathways
levels(cellchat@idents) 
#在层次绘图的时候，第一列显示的细胞类型的数目，一般有几个细胞就选几个
vertex.receiver = seq(1,5)
# Check all available pathways in your CellChat object
available_pathways <- cellchat@netP$pathways
print("Available pathways:")
print(available_pathways)
#显示的信号通路
pathways.show <- "TGFb"

pdf("4.重要通路通讯_贝壳图TGFb.pdf")
par(mfrow=c(1,1),xpd=TRUE)
netVisual_aggregate(cellchat,
                    signaling = pathways.show,
                    layout = "circle")
dev.off()
pdf("5.重要通路通讯_弦图.pdf")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat,
                    signaling = pathways.show,
                    layout = "chord")
dev.off()
#绘制通讯热图
pdf("6.通讯热图.pdf")
netVisual_heatmap(cellchat,
                  signaling = pathways.show,
                  color.heatmap = "Reds")
dev.off()

#绘制配体受体对整个信号通路的贡献度
pdf("7.配体受体贡献度.pdf")
netAnalysis_contribution(cellchat,
                         signaling = pathways.show)
dev.off()

#可视化单个配体受体通讯网络
pairLR.CXCL <- extractEnrichedLR(cellchat,
                                 signaling = pathways.show,
                                 geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] 
vertex.receiver = seq(1,5) 
pdf("8.单个配体受体通讯网络_贝壳图.pdf")
netVisual_individual(cellchat,
                     signaling = pathways.show,
                     pairLR.use = LR.show,
                     layout = "circle")
dev.off()

pdf("9.单个配体受体通讯网络_弦图.pdf")
netVisual_individual(cellchat,
                     signaling = pathways.show,
                     pairLR.use = LR.show,
                     layout = "chord")
dev.off()

#可视化由多种配体-受体或信号通路介导的细胞间通讯
levels(cellchat@idents)
pdf("10.多配体受体介导气泡图.pdf")
netVisual_bubble(cellchat,
                 sources.use = 5,
                 remove.isolate = FALSE)
dev.off()
pdf("11.多配体受体介导弦图.pdf")
netVisual_chord_gene(cellchat,
                     sources.use = 5,
                     lab.cex = 0.5,
                     legend.pos.y = 30)
dev.off()

#识别细胞群体中的信号作用（例如，显性发送者、接收者）以及主要贡献的信号
cellchat <- netAnalysis_computeCentrality(cellchat,
                                          slot.name = "netP")
pdf("12.中心性得分.pdf")
netAnalysis_signalingRole_network(cellchat,
                                  signaling = pathways.show,
                                  width = 8,
                                  height = 2.5,
                                  font.size = 10)
dev.off()

pdf("13.发送者和接收者散点图.pdf")
netAnalysis_signalingRole_scatter(cellchat)
dev.off()

#识别对某些细胞群体出向或入向信号传递贡献最大的信号
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
pdf("14.信号通路在细胞通讯中的作用.pdf",width = 12,height = 6)
ht1 + ht2
dev.off()

#识别并可视化分泌细胞的输出通信模式
#推断模式数量
pdf("15.selectK推断_输出.pdf")
selectK(cellchat, pattern = "outgoing")
dev.off()
#输出模式数量为4时曲线骤降
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat,
                                          pattern = "outgoing",
                                          k = nPatterns)
pdf("16.通讯模式热图_输出.pdf")
identifyCommunicationPatterns(cellchat,
                              pattern = "outgoing",
                              k = nPatterns)
dev.off()
pdf("17.通讯模式桑基图_输出.pdf")
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()
pdf("18.通讯模式点图_输出.pdf")
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()

#识别并可视化分泌细胞的输入通信模式
pdf("19.selectK推断_输入.pdf")
selectK(cellchat, pattern = "incoming")
dev.off()
#输出模式数量为4时曲线骤降
nPatterns = 7
cellchat <- identifyCommunicationPatterns(cellchat,
                                          pattern = "incoming",
                                          k = nPatterns)
pdf("20.通讯模式热图_输入.pdf")
identifyCommunicationPatterns(cellchat,
                              pattern = "incoming",
                              k = nPatterns)
dev.off()
pdf("21.通讯模式桑基图_输入.pdf")
netAnalysis_river(cellchat,
                  pattern = "incoming")
dev.off()
pdf("22.通讯模式点图_输入.pdf")
netAnalysis_dot(cellchat,
                pattern = "incoming")
dev.off()

#信号网络的多维和分类学习分析
cellchat <- computeNetSimilarity(cellchat, type = "functional")

cellchat <- netEmbedding(cellchat,
                         type = "functional",
                         umap.method = "uwot")

cellchat <- netClustering(cellchat,
                          type = "functional")

pdf("23.二维分类图.pdf")
netVisual_embedding(cellchat, 
                    type = "functional",
                    label.size = 3.5)
dev.off()

#根据相似性识别信号组
cellchat <- computeNetSimilarity(cellchat,
                                 type = "structural")
cellchat <- netEmbedding(cellchat,
                         type = "structural",
                         umap.method = "uwot")

cellchat <- netClustering(cellchat,
                          type = "structural")
pdf("24.二维信号组图.pdf")
netVisual_embedding(cellchat,
                    type = "structural",
                    label.size = 3.5)
dev.off()
pdf("25.分组二维信号图.pdf")
netVisual_embeddingZoomIn(cellchat,
                          type = "structural", nCol = 2)
dev.off()
qsave(cellchat,
      file = "cellchat.qs",
      nthreads = detectCores())



library(Nebulosa)
library(Seurat)
library(ggplot2)
library(viridis)


sce.all.int <- readRDS("results/seurat_integrated_labelled.rds")

p1 <- plot_density(sce.all.int,
                   c("SCARA5","PYCR1","KIF20A","AURKA","UHRF1"),
                   size=0.3)
p1

p4 <- plot_density(sce.all.int, c("SCARA5","PYCR1","KIF20A","AURKA","UHRF1"), joint = TRUE,combine = FALSE)
p4





# Define genes
genes <- c("SCARA5", "PYCR1", "KIF20A", "AURKA", "UHRF1")

# Function to create publication-quality individual plot WITH axes

# Correct function - use the scale replacement method but suppress warnings
create_individual_plot <- function(seurat_obj, gene_name) {
  suppressWarnings({
    plot_density(seurat_obj, gene_name, size = 0.4) +
      scale_color_gradientn(
        colors = c("#000033", "#0D0887", "#6A00A8", "#B12A90", 
                   "#E16462", "#FCA636", "#F0F921"),  # Custom plasma-like with your dark start
        name = "Expression\nDensity",
        na.value = "#000033"
      ) +
      labs(
        title = gene_name,
        x = "UMAP 1",
        y = "UMAP 2"
      ) +
      theme_classic() +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 9),
        plot.margin = margin(15, 15, 15, 15),
        panel.grid = element_blank()
      )
  })
}


# Now create the plots
p_SCARA5 <- create_individual_plot(sce.all.int, "SCARA5")
p_PYCR1 <- create_individual_plot(sce.all.int, "PYCR1")
p_KIF20A <- create_individual_plot(sce.all.int, "KIF20A")
p_AURKA <- create_individual_plot(sce.all.int, "AURKA")
p_UHRF1 <- create_individual_plot(sce.all.int, "UHRF1")

# Display plots
print(p_SCARA5)
print(p_PYCR1)
print(p_KIF20A)
print(p_AURKA)
print(p_UHRF1)



# If you also want PNG versions alongside PDFs
for(gene in genes) {
  # Create plot
  p <- create_individual_plot(sce.all.int, gene)
  
  # Save PDF (60pt x 60pt)
  ggsave(
    filename = paste0("density_", gene, ".pdf"),
    plot = p,
    width = 200,
    height = 200,
    units = "mm",
    dpi = 600,
    device = "pdf"
  )
  
  # Save PNG (same size)
  ggsave(
    filename = paste0("density_", gene, ".png"),
    plot = p,
    width = 200,
    height = 200,
    units = "mm",
    dpi = 600,
    device = "png"
  )
  
  cat("Saved:", paste0("density_", gene, ".pdf"), "and", paste0("density_", gene, ".png"), "\n")
}





library(Nebulosa)
library(Seurat)
library(ggplot2)
library(viridis)

# Set different thresholds for each gene based on their expression patterns
thresholds <- list(
  "SCARA5" = 0.15,   # Higher threshold for SCARA5 (max ~0.4)
  "PYCR1" = 0.03,    # Medium threshold for PYCR1 (max ~0.12)
  "KIF20A" = 0.025,   # Medium threshold for KIF20A (max ~0.125)
  "AURKA" = 0.015,    # Lower threshold for AURKA (max ~0.04)
  "UHRF1" = 0.05     # Medium threshold for UHRF1 (max ~0.15)
)

for(gene in names(p4_enhanced)) {
  threshold <- thresholds[[gene]]
  
  p <- suppressWarnings({
    p4_enhanced[[gene]] + 
      scale_color_viridis_c(
        option = "plasma", 
        name = "Expression\nDensity",
        begin = 0.0,
        end = 0.98,
        limits = c(threshold, NA),
        na.value = "grey80",
        oob = scales::squish
      ) +
      labs(
        title = paste0(gene, " (threshold: ", threshold, ")"),
        x = "UMAP 1", 
        y = "UMAP 2"
      ) +
      theme_classic() +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        legend.position = "right",
        panel.grid = element_blank()
      )
  })
  
  print(p)
  ggsave(paste0("density_", gene, "_custom_threshold.pdf"), p, 
         width = 200, height = 200, units = "mm", dpi = 300)
}


test_threshold <- function(gene_name, threshold_value) {
  suppressWarnings({
    p4_enhanced[[gene_name]] + 
      scale_color_viridis_c(
        option = "plasma", 
        name = "Expression\nDensity",
        begin = 0.0,
        end = 0.98,
        limits = c(threshold_value, NA),
        na.value = "grey80",           # This will now work
        oob = scales::censor           # This creates NAs for out-of-bounds values
      ) +
      labs(
        title = paste0(gene_name),
        x = "UMAP 1", 
        y = "UMAP 2"
      ) +
      theme_classic() +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "right",
        panel.grid = element_blank()
      )
  })
}

# Test it
print(test_threshold("SCARA5", 0.1))
print(test_threshold("PYCR1", 0.025))
print(test_threshold("KIF20A", 0.025))
print(test_threshold("AURKA", 0.01))
print(test_threshold("UHRF1", 0.05))
# Save the plots with custom thresholds
for(gene in names(thresholds)) {
  threshold <- thresholds[[gene]]
  
  p <- test_threshold(gene, threshold)
  
  ggsave(
    filename = paste0("density_", gene, "_custom_threshold.pdf"),
    plot = p,
    width = 200,
    height = 200,
    units = "mm",
    dpi = 300,
    device = "pdf"
  )
  
  cat("Saved:", paste0("density_", gene, "_custom_threshold.pdf"), "\n")
}

