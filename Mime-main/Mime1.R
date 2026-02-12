# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

depens<-c('GSEABase', 'GSVA', 'cancerclass', 'mixOmics', 'sparrow', 'sva' , 'ComplexHeatmap' )
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))  BiocManager::install(depen,update = FALSE)
}

if (!requireNamespace("CoxBoost", quietly = TRUE))
  devtools::install_github("binderh/CoxBoost")

if (!requireNamespace("fastAdaboost", quietly = TRUE))
  devtools::install_github("souravc83/fastAdaboost")

if (!requireNamespace("Mime", quietly = TRUE))
  devtools::install_github("l-magnificence/Mime")

library(Mime1)

load("External data/Example.cohort.Rdata")
list_train_vali_Data[["Dataset1"]][1:5,1:5]
#>               ID    OS.time OS   MT-CO1   MT-CO3
#>  TCGA.DH.A66B.01 1281.65322  0 13.77340 13.67931
#>  TCGA.HT.7607.01   96.19915  1 14.96535 14.31857
#>  TCGA.DB.A64Q.01  182.37755  0 13.90659 13.65321
#>  TCGA.DU.8167.01  471.97707  0 14.90695 14.59776
#>  TCGA.HT.7610.01 1709.53901  0 15.22784 14.62756



load("External data/Example.ici.Rdata")
list_train_vali_Data[["training"]][1:5,1:5]
#>                                   ID Var      FTH1   EEF1A1      ACTB
#>                      SAMf2ce197162ce   N 10.114846 4.817746 11.230180
#>                           ERR2208915   Y  2.044180 5.038854  3.977902
#> G138701_RCCBMS-00141-T_v1_RNA_OnPrem   Y  5.406008 5.341635  5.366668
#>                      SAMe41b1e773582   N  9.215794 4.707360 11.412721
#>                      SAM5ffd7e4cd794   N  9.003710 3.908884 10.440559

load("External data/genelist.Rdata")
class(genelist)
#> [1] "MYC"    "CTNNB1" "JAG2"   "NOTCH1" "DLL1"   "AXIN2"  "PSEN2"  "FZD1"   "NOTCH4" "LEF1"   "AXIN1"  "NKD1"   "WNT5B"
#>[14] "CUL1"   "JAG1"   "MAML1"  "KAT2A"  "GNAI1"  "WNT6"   "PTCH1"  "NCOR2"  "DKK4"   "HDAC2"  "DKK1"   "TCF7"   "WNT1"
#>[27] "NUMB"   "ADAM17" "DVL2"   "PPARD"  "NCSTN"  "HDAC5"  "CCND2"  "FRAT1"  "CSNK1E" "RBPJ"   "FZD8"   "TP53"   "SKP2"
#>[40] "HEY2"   "HEY1"   "HDAC11"

rm(list=ls())
library(Mime1)
load("Dataset1.Rdata")
load("Dataset2.Rdata")
load("genelist.Rdata")
list_Data <- list(Train_tcga = Dataset1, Test_geo = Dataset2)
save(list_Data, file = "My_cohort.Rdata")
# load("External data/Example.cohort.Rdata")
# load("External data/example.genelist.Rdata")
res <- ML.Dev.Prog.Sig(
  train_data = list_Data$Train_tcga,
  list_train_vali_Data = list_Data,
  unicox.filter.for.candi = TRUE,
  unicox_p_cutoff = 0.15,  # Increase from 0.05 to retain ≥2 genes
  candidate_genes = Ferro_genes,
  mode = 'all',
  nodesize = 5,
  seed = 2676
)


cindex_dis_all(res,
               validate_set = names(list_Data)[-1],
               order =names(list_Data),
               width = 0.35)


cindex_dis_select(res,
                  model="StepCox[forward] + RSF",
                  order= names(list_Data))

survplot <- vector("list",2)
for (i in c(1,2)) {
  print(survplot[[i]]<-rs_sur(res, model_name = "StepCox[forward] + RSF",
                              dataset = names(list_Data)[i],
                              #color=c("blue","green"),
                              median.line = "hv",
                              cutoff = 0.75,
                              conf.int = T,
                              xlab="Day",pval.coord=c(2000,0.9)))
}

aplot::plot_list(gglist=survplot,ncol=2)



all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_Data[["Train_tcga"]],
                             inputmatrix.list = list_Data,mode = 'all',AUC_time = 1,
                             auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_Data[["Train_tcga"]],
                             inputmatrix.list = list_Data,mode = 'all',AUC_time = 3,
                             auc_cal_method="KM")
all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_Data[['Train_tcga']],
                             inputmatrix.list = list_Data,mode = 'all',AUC_time = 5,
                             auc_cal_method="KM")


# 1 Year AUC predicted by all models
auc_dis_all(all.auc.1y,
            dataset = names(list_Data),
            validate_set=names(list_Data)[-1],
            order= names(list_Data),
            width = 0.35,
            year=1)
# 3 Year AUC predicted by all models
auc_dis_all(all.auc.3y,
            dataset = names(list_Data),
            validate_set=names(list_Data)[-1],
            order= names(list_Data),
            width = 0.35,
            year=3)
# 5 Year AUC predicted by all models
auc_dis_all(all.auc.5y,
            dataset = names(list_Data),
            validate_set=names(list_Data)[-1],
            order= names(list_Data),
            width = 0.35,
            year=5)
# Plot ROC of specific model among different datasets:
roc_vis(all.auc.1y,
        model_name = "StepCox[forward] + RSF",
        dataset = names(list_Data),
        order= names(list_Data),
        anno_position=c(0.65,0.55),
        year=1)

roc_vis(all.auc.3y,
        model_name = "StepCox[forward] + RSF",
        dataset = names(list_Data),
        order= names(list_Data),
        anno_position=c(0.65,0.55),
        year=3)

roc_vis(all.auc.5y,
        model_name = "StepCox[forward] + RSF",
        dataset = names(list_Data),
        order= names(list_Data),
        anno_position=c(0.65,0.55),
        year=5)



# Plot 1, 3, and 5-year AUC of specific model among different datasets
auc_dis_select(list(all.auc.1y,all.auc.3y,all.auc.5y),
               model_name="StepCox[forward] + RSF",
               dataset = names(list_Data),
               order= names(list_Data),
               year=c(1,3,5))


# Meta-analysis of univariate cox regression for specific model
unicox.rs.res <- cal_unicox_ml_res(res.by.ML.Dev.Prog.Sig = res,
                                   optimal.model = "StepCox[forward] + RSF",
                                   type ='categorical')
metamodel <- cal_unicox_meta_ml_res(input = unicox.rs.res)
meta_unicox_vis(metamodel,
                dataset = names(list_Data))


# Comparison with previously pblished models

rs.glioma.lgg.gbm <- cal_RS_pre.prog.sig(use_your_own_collected_sig = F,type.sig = c('LGG','GBM','Glioma'),
                                         list_input_data = list_Data)

HR_com(rs.glioma.lgg.gbm,
       res,
       model_name="StepCox[forward] + RSF",
       dataset=names(list_train_vali_Data),
       type = "categorical")


cc.glioma.lgg.gbm <- cal_cindex_pre.prog.sig(use_your_own_collected_sig = F,type.sig = c('Glioma','LGG','GBM'),
                                             list_input_data = list_train_vali_Data)


cindex_comp(cc.glioma.lgg.gbm,
            res,
            model_name="StepCox[forward] + RSF",
            dataset=names(list_train_vali_Data))


auc.glioma.lgg.gbm.1 <- cal_auc_pre.prog.sig(use_your_own_collected_sig = F,
                                             type.sig = c('Glioma','LGG','GBM'),
                                             list_input_data = list_train_vali_Data,AUC_time = 1,
                                             auc_cal_method = 'KM')


auc_comp(auc.glioma.lgg.gbm.1,
         all.auc.1y,
         model_name="StepCox[forward] + plsRcox",
         dataset=names(list_train_vali_Data))



#Immune infiltration analysis
devo <- TME_deconvolution_all(list_Data)
immuno_heatmap(res,
               devo,
               model_name="StepCox[forward] + RSF",
               dataset="Train_tcga")





#2. Construct predicting models for response
load("./Example.ici.Rdata")
load("./genelist.Rdata")
res.ici <- ML.Dev.Pred.Category.Sig(train_data = list_Data$Train_tcga,
                                    list_train_vali_Data = list_Data,
                                    candidate_genes = Ferro_genes,
                                    methods = c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass'),
                                    seed = 5201314,
                                    cores_for_parallel = 60
)



auc_vis_category_all(res.ici,dataset = c("training","validation"),
                     order= c("training","validation"))




plot_list<-list()
methods <- c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass')
for (i in methods) {
  plot_list[[i]]<-roc_vis_category(res.ici,model_name = i,dataset = c("training","validation"),
                                   order= c("training","validation"),
                                   anno_position=c(0.4,0.25))
}
aplot::plot_list(gglist=plot_list,ncol=3)


auc.other.pre <- cal_auc_previous_sig(list_train_vali_Data = list_train_vali_Data,seed = 5201314,
                                      train_data = list_train_vali_Data$training,
                                      cores_for_parallel = 32)



auc_category_comp(res.ici,
                  auc.other.pre,
                  model_name="svmRadialWeights",
                  dataset=names(list_train_vali_Data))



#3. Core feature selection
rm(list=ls())
load("My_cohort.Rdata")
load("genelist.Rdata")
load("External data/Example.cohort.Rdata")
load("External data/example.genelist.Rdata")

res.feature.all <- ML.Corefeature.Prog.Screen(InputMatrix = list_Data$Test_geo,
                                              candidate_genes = Ferro_genes,
                                              mode = 'all_without_SVM',
                                              nodesize =5,
                                              seed = 5201314)

core_feature_select(res.feature.all)

core_feature_rank(res.feature.all, top=5)


dataset_col<-c("#3182BDFF","#E6550DFF")
corplot <- list()
for (i in c(1:2)) {
  print(corplot[[i]]<-cor_plot(list_Data[[i]],
                               dataset=names(list_Data)[i],
                               color = dataset_col[i],
                               feature1="SCARA5",
                               feature2="PYCR1",
                               method="pearson"))
}
aplot::plot_list(gglist=corplot,ncol=2)


survplot <- vector("list",2)
for (i in c(1:2)) {
  print(survplot[[i]]<-core_feature_sur("SCARA5",
                                        InputMatrix=list_Data[[i]],
                                        dataset = names(list_Data)[i],
                                        #color=c("blue","green"),
                                        median.line = "hv",
                                        cutoff = 0.5,
                                        conf.int = T,
                                        xlab="Day",pval.coord=c(2000,0.9)))
}
aplot::plot_list(gglist=survplot,ncol=2)









