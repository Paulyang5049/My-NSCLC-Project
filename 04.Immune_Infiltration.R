rm(list = ls())
load("common_genes.RData")
load("step2.RData")
load("geneset_cell.RData")
load("TISIDB肿瘤浸润淋巴细胞基因集.rdata")
# extract rows in exp that are in common_genes

exp_1 <- exp[common_genes$common_genes,]

# save exp_1 as a xlsx file
library(GSVA)
library(limma)
library(tibble)
library(magrittr)
library(reshape2)
library(ggplot2)
library(ggpubr)
ssGSEA_matrix <- ssgseaParam(exp_1,
                      tisidb_cell,
                      assay = NA_character_,
                      annotation = NA_character_,
                      minSize = 2,
                      maxSize = Inf,
                      alpha = 0.25,
                      normalize = TRUE)
gsva_es <- gsva(ssGSEA_matrix)


tmp1 <- gsva_es%>%t()%>%as.data.frame() %>%
  rownames_to_column("Sample")
tmp1$Group <- Group
tmp1 <- melt(tmp1)
colnames(tmp1) <- c("Sample","Group","Celltype","Score")

pdf("./Plot/04. Immune Infiltration/immune_infiltration_box.pdf",width = 10,height = 6)
ggplot(tmp1,aes(Celltype,Score)) + 
  geom_boxplot(aes(fill = Group),outlier.shape = 21)+
  theme_bw() +
  labs(x = NULL, y = "Score") +
  scale_fill_manual(values = c("#2ecc71", "#e67e22"))+
  stat_compare_means(aes(group = Group,label = after_stat(p.signif)),
                     method = "wilcox.test",
                     hide.ns = T)+
  theme(plot.margin=unit(c(1,1,1,1),'cm'),
        plot.title = element_text(size = 12,color="black",hjust = 0.5),
        axis.title = element_text(size = 12,color ="black"), 
        axis.text = element_text(size= 12,color = "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1 ),
        panel.grid=element_blank(),
        legend.position = "top",
        legend.text = element_text(size= 12),
        legend.title= element_text(size= 12))
dev.off()



# 加载必要的库
library(pheatmap)
library(RColorBrewer)
# 注释数据框
annotation_col <- data.frame(Group = factor(Group))
rownames(annotation_col) <- colnames(gsva_es)

# 注释颜色
annotation_colors <- list(
  Group = c(Control = "skyblue", Tumor = "salmon")
)
# 配色方案
my_colors <- colorRampPalette(c("blue", "white", "red"))(50)
# 绘制热图
pheatmap(gsva_es,
         color = my_colors,
         cluster_rows = TRUE,  # 行聚类
         cluster_cols = TRUE,  # 列聚类
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize = 10,
         fontsize_row = 10,
         fontsize_col = 10,
         scale = "row",  # 对行进行缩放
         main = "Immune Cell Infiltration"
)

# 保存热图到文件
pdf("./Plot/04. Immune Infiltration/immune_infiltration_heatmap.pdf", width = 8, height = 10)
