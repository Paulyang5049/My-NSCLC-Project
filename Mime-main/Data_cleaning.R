library(readxl)
Driver <- read_excel("/Volumes/Haoran's SSD/Scientific research/My research/LUAD_我的第一篇文章/Mime-main/Driver.xlsx")
Suppressor <- read_excel("/Volumes/Haoran's SSD/Scientific research/My research/LUAD_我的第一篇文章/Mime-main/Suppressor.xlsx")

Driver <- Driver[[1]]
Suppressor <- Suppressor[[1]]
# Combine them into one vector
Ferro_genes <- c(Driver, Suppressor)

save(Driver, Suppressor, Ferro_genes, file = "genelist.Rdata")

rm(list=ls())
load("genelist.Rdata")
load('LUAD_prog_data.Rdata')

# # Assuming samples are columns and genes are rows in tpm_data
# tumor_samples <- mrna_expr_tpm[, Group_tcga == "Tumor"]
# tumor_samples <- t(tumor_samples)

library(tidyverse)
luad_surv_gdc <- clinical_indexed %>%
  dplyr::select(submitter_id, bcr_patient_barcode, vital_status,days_to_last_follow_up,
         days_to_death) %>%
  dplyr::mutate(OS_time = if_else(vital_status == "Dead", days_to_death,
                           days_to_last_follow_up),
         vital_status = if_else(vital_status == "Dead",1,0)) %>%
  dplyr::select(sample=submitter_id, OS=vital_status,
         `LUAD_PATIENT`=bcr_patient_barcode,OS.time=OS_time
  )

rownames(luad_surv_gdc) <- NULL
luad_surv_gdc %>% distinct(`LUAD_PATIENT`,.keep_all = T) %>% dim()

log_TPM <- log2(mrna_expr_tpm + 1)
t_tpm <- t(log_TPM)
t_tpm <- as.data.frame(t_tpm)

# 将LUAD_count的行名设置为第一列
t_tpm <- cbind(ID = rownames(t_tpm), t_tpm)
# 假设 luad_surv_gdc 和 t_tpm 已经加载为数据框
# 截取 t_tpm 中 ID 列的前 12 个字符
t_tpm$ID <- substr(t_tpm$ID, 1, 12)

# 使用前 12 个字符与 luad_surv_gdc 的 sample 列匹配
matched_rows <- match(luad_surv_gdc$sample, t_tpm$ID)

# 提取 t_tpm 中匹配的行
t_tpm_matched <- t_tpm[matched_rows, ]

# 确保 luad_surv_gdc 和 t_tpm_matched 是正确匹配的顺序
# 将 luad_surv_gdc 的 OS.time 和 OS 列加入到 t_tpm_matched 的第二和第三列
t_tpm_matched <- cbind(t_tpm_matched[, 1], luad_surv_gdc$OS.time, luad_surv_gdc$OS, t_tpm_matched[, -1])

# 重命名列
colnames(t_tpm_matched)[2:3] <- c("OS.time", "OS")

# 查看结果
View(t_tpm_matched)

# 将第一列重命名为 ID
colnames(t_tpm_matched)[1] <- "ID"

# Assuming samples are columns and genes are rows in tpm_data
tumor_samples <- t_tpm_matched[Group_tcga == "Tumor",]
# tumor_samples <- t(tumor_samples)
# 删除含有缺失值的行
tumor_samples_clean <- na.omit(tumor_samples)
tumor_samples_clean <- tumor_samples_clean[tumor_samples_clean$OS.time != 0, ]
tumor_samples <- tumor_samples_clean

# randomly select 100 rows in tumor_samples
set.seed(520)
tumor_samples <- tumor_samples[sample(nrow(tumor_samples_clean), 100),]

Dataset1 <- tumor_samples
rownames(Dataset1) <- NULL
Dataset1$OS.time <- as.double(Dataset1$OS.time)
save(Dataset1, file = "Dataset1.Rdata")







# GEO
library(GEOquery)#直接调用
rm(list=ls())

eSet <- getGEO("GSE13213",
               destdir = '.',
               getGPL = T)

exprSet = exprs(eSet[[1]])#表达量矩阵
fdata = fData(eSet[[1]])#需要getGPL = T,平台信息,注释文件包含探针与基因的对应关系
pdata = pData(eSet[[1]]) #表型信息，里面包含了需要的分组信息
samplenames <- sampleNames(eSet[[1]])

dim(exprSet)

dim(fdata)

dim(pdata)

length(samplenames)


#找出ID对应的基因
genes <- fdata[,c("ID","GENE_SYMBOL")]
exprSet <- data.frame(exprSet)
exprSet$ID <- rownames(exprSet)
exprSet <- merge(exprSet,genes,by.x = "ID",by.y = 'ID')
exprSet <- exprSet[,-1]
dim(exprSet)    #最后一列为基因

# Remove duplicate rows based on GENE_SYMBOL
exprSet <- exprSet[!duplicated(exprSet[, "GENE_SYMBOL"]), ]

# Assign row names
rownames(exprSet) <- exprSet[, "GENE_SYMBOL"]

# Optionally, remove the 118th column if it was just for row names
exprSet <- exprSet[, -118]


t_exp <- t(exprSet)
t_exp <- as.data.frame(t_exp)


GEO_surv <- pdata %>%
  dplyr::select(ID = geo_accession, OS = characteristics_ch1.8,
         OS.time = characteristics_ch1.9) %>%
  dplyr::mutate(OS = if_else(OS == "Status: Dead", 1, 0)) %>%
  dplyr::select(ID, OS.time, OS)

# Assuming GEO_surv is your dataframe
GEO_surv$OS.time <- gsub("Survival \\(days\\): ", "", GEO_surv$OS.time)

# View the updated dataframe
View(GEO_surv)
as.numeric(GEO_surv$OS.time)


# Assuming 'exprSet_na' is your matrix or dataframe and 'pdata$geo_accession' contains the accession IDs
AD_col <- exprSet[, colnames(exprSet) %in% pdata$geo_accession]

# View the extracted columns
View(AD_col)

# Assuming 'GEO_surv' has row names and 'pdata$geo_accession' contains the row names to match
AD_surv <- GEO_surv[rownames(GEO_surv) %in% pdata$geo_accession, ]

# View the extracted rows
View(AD_surv)

Dataset2 <- t(AD_col)
Dataset2 <- as.data.frame(Dataset2)
# Assuming AD_surv and Dataset2 are your dataframes

# Combine the columns from AD_surv and the first three columns of Dataset2
combined_dataset <- cbind(AD_surv, Dataset2)

# View the combined dataset
View(combined_dataset)

Dataset2 <- combined_dataset
Dataset2$OS.time <- as.double(Dataset2$OS.time)
# rename the first column as ID
colnames(Dataset2)[1] <- "ID"
rownames(Dataset2) <- NULL
Dataset2 <- Dataset2[, c(1, 3, 2, 4:ncol(Dataset2))]
missing_ferro_genes <- setdiff(Ferro_genes, colnames(Dataset2))

save(Dataset2, file = "Dataset2.Rdata")
