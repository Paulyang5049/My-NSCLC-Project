# Remove all objects from the current workspace
rm(list = ls())

# Load the data
load("GSE32863_deg.RData")
load("LUAD.RData")


# Show the number of changes in each dataset
table(deg$change)
table(LUAD_deg$change)

# Output
# geo_deg:
#   down stable     up 
#   936  47218    649

# LUAD_deg:
#   down stable     up 
#   2597  50231   5559

# Load libraries
library(VennDiagram)
library(gridExtra)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# Filter out genes that are "stable"
geo_deg <- deg %>%
  mutate(gene = rownames(deg)) %>%
  filter(change != "stable")

tcga_deg <- LUAD_deg %>%
  mutate(gene = rownames(LUAD_deg)) %>%
  filter(change != "stable")

# Extract genes that are "down" or "up"
geo_genes <- geo_deg$gene
tcga_genes <- tcga_deg$gene

# Venn Diagram
venn.plot <- venn.diagram(
  x = list(
    `TCGA-LUAD` = tcga_genes,
    `GSE32863` = geo_genes
  ),
  category.names = c("TCGA-LUAD", "GSE32863"),
  fill = c("#88CEEA", "#FFC0CA"),  # Colors: Custom colors
  alpha = 0.4,                     # Transparency
  cex = 2,                         # Text size for numbers
  cat.cex = 2,                     # Text size for category labels
  cat.col = c("#88CEEA", "#FFC0CA"), # Colors for category labels
  cat.pos = c(-20, 20),            # Position of category names
  cat.dist = c(0.055, 0.055),      # Distance of category names from circles
  lwd = 2,                         # Line width for circles
  lty = "solid",                   # Line type
  filename = NULL                  # Prevents automatic saving
)

# Plot the Venn Diagram
grid.newpage()
grid.draw(venn.plot)

# Extract common genes
common_genes <- intersect(geo_genes, tcga_genes)


#下面我们进行GO和KEGG，我们只需要一串gene symbol即可（需要转换为Gene ID）
gene <- rownames(common_genes)

gene = bitr(common_genes,, 
            fromType="SYMBOL", 
            toType="ENTREZID", 
            OrgDb="org.Hs.eg.db")
table(duplicated(gene))


#GO/KEGG Cluster Analysis
GO_database <- 'org.Hs.eg.db'
KEGG_database <- 'hsa'
class(gene)
GO<-enrichGO(gene$ENTREZID,#GO富集分析
             OrgDb = GO_database,
             keyType = "ENTREZID",#设定读取的gene ID类型
             ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
             pvalueCutoff = 0.05,#设定p值阈值
             qvalueCutoff = 0.05,#设定q值阈值
             readable = T)
GO
x <- data.frame(GO@result)

KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
y <- data.frame(KEGG@result)
common_genes <- as.data.frame(common_genes)
save(common_genes, geo_deg, tcga_deg, GO, KEGG, file = "common_genes.RData")
#一些可视化

barplot(GO,split="ONTOLOGY")+
  facet_grid(ONTOLOGY~., scale= 'free')+ 
  scale_fill_gradient(low = 'turquoise', high = 'coral')+
  theme(strip.text = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 12, hjust = 1))

barplot(KEGG,showCategory = 15,title = 'KEGG Pathway')+
  theme(axis.text.y = element_text(size = 12))+
  scale_fill_gradient(low = 'turquoise', high = 'coral')


dotplot(GO, split="ONTOLOGY",showCategory = 5)+
  facet_grid(ONTOLOGY~., scale="free")+
  scale_fill_gradient(low = 'turquoise', high = 'coral')





#GSEA函数
rm(list = ls())
load("deg.RData")
data <- rownames(allDiff)
data = bitr(data, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
data <- dplyr::distinct(data,SYMBOL,.keep_all=TRUE)
data01 <- data.frame(logFC=allDiff$logFC,SYMBOL = rownames(allDiff))
data01 <- merge(data01,data,by="SYMBOL")
#目的就是获得一个包含symbol GeneID以及logFC的数据框

geneList <- data01$logFC
names(geneList) = data01$ENTREZID
geneList = sort(geneList, decreasing = TRUE)
geneList[1:10]


#https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
hallmarks <- read.gmt("h.all.v2023.1.Hs.entrez.gmt")
res <- GSEA(geneList,TERM2GENE =hallmarks)

res1<- setReadable(res,OrgDb=org.Hs.eg.db, keyType = "ENTREZID")

dotplot(res1,
        showCategory=12,split=".sign")+
  facet_grid(~.sign)

data <- data.frame(res@result)
gseaplot2(res,"HALLMARK_MYC_TARGETS_V2",color = "red",pvalue_table = T)
gseaplot2(res,6,color = "red",pvalue_table = T)
ridgeplot(res,showCategory = 10)
gseaplot2(res, geneSetID = 1:5)






####demo04:GSVA####
#options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
#if(!require("pheatmap")) install.packages("pheatmap")
#if(!require("limma")) BiocManager::install("limma",update = F,ask = F)
rm(list = ls())
library(GSVA)
library(limma)
#https://zhuanlan.zhihu.com/p/518145829

load("data_expr.Rdata")

kegggmt <- read.gmt("c2.cp.v2023.1.Hs.symbols.gmt")
colnames(kegggmt)
kegg_list = split(kegggmt$gene, kegggmt$term)
kegg1 <- gsva(expr=as.matrix(data), kegg_list, kcdf="Gaussian",method = "gsva",parallel.sz=1)

exprSet <- kegg1

group <- c(rep("con",18),rep("treat",19)) 
group <- factor(group,levels = c("con","treat"),ordered = F)

design <- model.matrix(~group)

colnames(design) <- levels(group)


fit <- lmFit(exprSet,design)
fit2 <- eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf) 