library(Nebulosa)
library(SCP)


sce.all.int <- readRDS("results/seurat_integrated_labelled.rds")

p1 <- plot_density(sce.all.int,c("SCARA5","PYCR1","KIF20A","AURKA","UHRF1"),size=0.3)
p1

p4 <- plot_density(sce.all.int, c("SCARA5","PYCR1","KIF20A","AURKA","UHRF1"), joint = TRUE,combine = FALSE)
p4


head(sce.all.int@meta.data)
CellDimPlot(sce.all.int, group.by = "celltype", reduction = "UMAP")

