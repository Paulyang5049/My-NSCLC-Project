rm(list=ls())
library(easyTCGA)
getclinical("TCGA-LUAD")
getmrnaexpr("TCGA-LUAD")


load('genelist.Rdata')
save(Dataset1, file = "Dataset1.Rdata")













