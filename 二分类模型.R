library(openxlsx)
library(seqinr)
library(plyr)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(caret)
library(mboost)
library(e1071)
library(BART)
library(MASS)
library(snowfall)
library(xgboost)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)
library(readr)
library(VIM)
library(caret)
library(rpart)
library(rpart.plot)
library(Metrics)
library(stringr)
library(rpart)
library(tibble)
library(bitops)
library(rattle)
library(rpart.plot)
library(RColorBrewer)
library(tidyverse)
library(pheatmap)
library(visNetwork)
library(ggpol)
library(ggplot2)
library(sparkline)
library(randomForest)
library(venn)
library(sparkline)
library(dplyr)
library(tidyverse)
library(caret)
library(DALEX)
library(gbm)
library(caret)
library(glmnet) 
library(xgboost)
library(DALEX)
library(gbm)
library(VennDiagram)
library(limma)  
library(neuralnet)
library(NeuralNetTools)
library(kernlab)
library(doParallel)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(data.table) 
library(ggplot2)
library(ggsci)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(ggbreak)
library(tidyr)
library(ggbreak)
library(edgeR)
library(limma)
library(survival)
library(survminer)
library(stringi)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(beepr)
library(pheatmap)
library(data.table)
library(ggsignif) 
library(RColorBrewer)
library(future.apply)
library(gplots)
library(DESeq2)
library(ggrepel)
library(Rcpp)
library(survivalsvm)
library(dplyr)
library(rms)
library(pec)
library(ggDCA)
library(glmnet)
library(foreign)
library(regplot)
library(randomForestSRC)
library(timeROC)
library(tidyr)
library(tibble)
library(caret)
library(regplot)
library(gbm)
library(tidyverse)
library(gbm)
library(obliqueRSF)
library(remotes)
library(aorsf)
library(xgboost)
library(party)
library(partykit)
library(UpSetR)
#####################################################
######读取数据
#####################################################
data0=read.table("data.txt",header = T,sep = "\t",check.names = F,row.names = 1)
#随机拆分训练集与测试集
data0=t(data0)
group0=sapply(strsplit(rownames(data0),"\\_"), "[", 2)
inTrain<-createDataPartition(y=data0[,2],p=0.7,list=F)
train<-data0[inTrain,]
train=as.data.frame(train)
rownames(train)=make.names(rownames(train))
test<-data0[-inTrain,]
test=as.data.frame(test)
rownames(test)=make.names(rownames(test))
#数据分组
grouptrain=sapply(strsplit(rownames(train),"\\_"), "[", 2)
grouptest=sapply(strsplit(rownames(test),"\\_"), "[", 2)
train=data.frame(group=grouptrain,train)
test=data.frame(group=grouptest,test)
datalist=list(Train=train,Test=test)

trainexp=as.matrix(train[,2:ncol(train)])
testexp=as.matrix(test[,2:ncol(test)])
result=data.frame()
modelgenes=list()
dir.create("ROC")
#####################################################
######1.Enet
#####################################################
for (alpha in seq(0.1,0.9,0.1)) {
cv.fit = cv.glmnet(x = trainexp,
                   y = grouptrain,
                   family = "binomial", alpha = alpha, nfolds = 100)
fit = glmnet(x = trainexp,
             y = grouptrain,
             family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
Enetgenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Enetgenes=Enetgenes[2:length(Enetgenes)]
modelname=paste0('Enet','[alpha=',alpha,']')
modelgenes[[modelname]]=Enetgenes
rs <- lapply(datalist,function(x){cbind(x[,1:2],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
}
#####################################################
######1.2Enet+glm
#####################################################
for (alpha in seq(0.1,0.9,0.1)) {
  cv.fit = cv.glmnet(x = trainexp,
                     y = grouptrain,
                     family = "binomial", alpha = alpha, nfolds = 100)
  fit = glmnet(x = trainexp,
               y = grouptrain,
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(data0)
  Enetgenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
  Enetgenes=Enetgenes[2:length(Enetgenes)]
  trainexp2=as.data.frame(trainexp[,Enetgenes])
  testexp2=as.data.frame(testexp[,Enetgenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test2=data.frame(group=grouptest,testexp2)
  datalist2=list(Train=train2,Test=test2)
  fit <- step(glm(formula = ifelse(grouptrain=="con",0,1) ~ .,
                  family = "binomial", 
                  data = as.data.frame(trainexp2)),trace = 0)
  fit$subFeature = colnames(trainexp2)
  glmgenes=names(coef(fit))[2:length(names(coef(fit)))]
  modelname=paste0('Enet','[alpha=',alpha,']+glm')
  modelgenes[[modelname]]=glmgenes
  rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(fit, type = 'response', as.data.frame(x)))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
    rownames_to_column('ID')
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
    plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
    text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
    dev.off()
  }
  AUCs$Model <- modelname
  result <- rbind(result,AUCs)
}
#####################################################
######1.3Enet+svm
#####################################################
for (alpha in seq(0.1,0.9,0.1)) {
  cv.fit = cv.glmnet(x = trainexp,
                     y = grouptrain,
                     family = "binomial", alpha = alpha, nfolds = 100)
  fit = glmnet(x = trainexp,
               y = grouptrain,
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(data0)
  Enetgenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
  Enetgenes=Enetgenes[2:length(Enetgenes)]
  trainexp2=as.data.frame(trainexp[,Enetgenes])
  testexp2=as.data.frame(testexp[,Enetgenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test2=data.frame(group=grouptest,testexp2)
  datalist2=list(Train=train2,Test=test2)
data <- as.data.frame(trainexp2)
cl <- makeCluster(16)
registerDoParallel(cl)
fit=rfe(x=data,
        y=as.numeric(as.factor(grouptrain)),
        rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
        methods="svmRadial")
stopCluster(cl)
modelname=paste0('Enet','[alpha=',alpha,']+svm')
modelgenes[[modelname]]=fit[["optVariables"]]
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(fit, as.data.frame(x)))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
}
#####################################################
######1.4Enet+plsRglmmodel
#####################################################
for (alpha in seq(0.1,0.9,0.1)) {
  cv.fit = cv.glmnet(x = trainexp,
                     y = grouptrain,
                     family = "binomial", alpha = alpha, nfolds = 100)
  fit = glmnet(x = trainexp,
               y = grouptrain,
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(data0)
  Enetgenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
  Enetgenes=Enetgenes[2:length(Enetgenes)]
  trainexp2=as.data.frame(trainexp[,Enetgenes])
  testexp2=as.data.frame(testexp[,Enetgenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test2=data.frame(group=grouptest,testexp2)
  datalist2=list(Train=train2,Test=test2)
cv.plsRglm.res = cv.plsRglm(formula = grouptrain ~ ., 
                            data = as.data.frame(trainexp2),
                            nt=10, verbose = FALSE)
fit <- plsRglm(ifelse(grouptrain=="con",0,1), 
               as.data.frame(trainexp2), 
               modele = "pls-glm-logistic",
               verbose = F, sparse = T)
fit$subFeature = colnames(trainexp2)
plsRglmgenes=rownames(fit$Coeffs)[fit$Coeffs!=0]
modelname=paste0('Enet','[alpha=',alpha,']+plsRglmmodel')
modelgenes[[modelname]]=plsRglmgenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(fit, type = "response", as.data.frame(x[,2:ncol(x)])))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
}
#####################################################
######1.5Enet+RF
#####################################################
for (alpha in seq(0.1,0.9,0.1)) {
  cv.fit = cv.glmnet(x = trainexp,
                     y = grouptrain,
                     family = "binomial", alpha = alpha, nfolds = 100)
  fit = glmnet(x = trainexp,
               y = grouptrain,
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(data0)
  Enetgenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
  Enetgenes=Enetgenes[2:length(Enetgenes)]
  trainexp2=as.data.frame(trainexp[,Enetgenes])
  testexp2=as.data.frame(testexp[,Enetgenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test2=data.frame(group=grouptest,testexp2)
  datalist2=list(Train=train2,Test=test2)
rf=randomForest(as.factor(grouptrain)~., data=trainexp2, ntree=500)
optionTrees=which.min(rf$err.rate[,1])
rf2=randomForest(as.factor(grouptrain)~., data=trainexp2, ntree=optionTrees)
importance=importance(x=rf2)
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes=names(rfGenes[rfGenes>0]) 
modelname=paste0('Enet','[alpha=',alpha,']+RF')
modelgenes[[modelname]]=rfGenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(rf2, type = "response", as.data.frame(x[,2:ncol(x)])))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
}
#####################################################
######1.6Enet+rpart
#####################################################
for (alpha in seq(0.1,0.9,0.1)) {
  cv.fit = cv.glmnet(x = trainexp,
                     y = grouptrain,
                     family = "binomial", alpha = alpha, nfolds = 100)
  fit = glmnet(x = trainexp,
               y = grouptrain,
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(data0)
  Enetgenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
  Enetgenes=Enetgenes[2:length(Enetgenes)]
  trainexp2=as.data.frame(trainexp[,Enetgenes])
  testexp2=as.data.frame(testexp[,Enetgenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test2=data.frame(group=grouptest,testexp2)
  datalist2=list(Train=train2,Test=test2)
mod1<-rpart(as.factor(grouptrain)~.,data = trainexp2,method = "class")
importances <- varImp(mod1)
jcsimportances=as.matrix(importances %>%
                           arrange(desc(Overall)))
JCSGenes=jcsimportances[order(jcsimportances[,"Overall"], decreasing = TRUE),]
JCSGenes=names(JCSGenes[JCSGenes>0])
modelname=paste0('Enet','[alpha=',alpha,']+rpart')
modelgenes[[modelname]]=JCSGenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(mod1, as.data.frame(x[,2:ncol(x)])))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
}
#####################################################
######1.7Enet+GBM
#####################################################
for (alpha in seq(0.1,0.9,0.1)) {
  cv.fit = cv.glmnet(x = trainexp,
                     y = grouptrain,
                     family = "binomial", alpha = alpha, nfolds = 100)
  fit = glmnet(x = trainexp,
               y = grouptrain,
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(data0)
  Enetgenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
  Enetgenes=Enetgenes[2:length(Enetgenes)]
  trainexp2=as.data.frame(trainexp[,Enetgenes])
  testexp2=as.data.frame(testexp[,Enetgenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test2=data.frame(group=grouptest,testexp2)
  datalist2=list(Train=train2,Test=test2)
fit <- gbm(formula = ifelse(grouptrain=="con",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 6)
best <- which.min(fit$cv.error)
fit <- gbm(formula = ifelse(grouptrain=="con",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001, n.cores = 8)
fit$subFeature = colnames(trainexp2)
GBMgenes=rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf>0]
modelname=paste0('Enet','[alpha=',alpha,']+GBM')
modelgenes[[modelname]]=GBMgenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(fit, type = "response", as.data.frame(x)))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
}
#####################################################
######1.8Enet+Ridge
#####################################################
for (alpha in seq(0.1,0.9,0.1)) {
  cv.fit = cv.glmnet(x = trainexp,
                     y = grouptrain,
                     family = "binomial", alpha = alpha, nfolds = 100)
  fit = glmnet(x = trainexp,
               y = grouptrain,
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(data0)
  Enetgenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
  Enetgenes=Enetgenes[2:length(Enetgenes)]
  trainexp2=as.data.frame(trainexp[,Enetgenes])
  testexp2=as.data.frame(testexp[,Enetgenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test2=data.frame(group=grouptest,testexp2)
  datalist2=list(Train=train2,Test=test2)
cv.fit = cv.glmnet(x = trainexp2,
                   y = grouptrain,
                   family = "binomial", alpha = 0, nfolds = 1000)
fit = glmnet(x = trainexp2,
             y = grouptrain,
             family = "binomial", alpha = 0, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
Ridgegenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Ridgegenes=Ridgegenes[2:length(Ridgegenes)]
modelname=paste0('Enet','[alpha=',alpha,']+Ridge')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
}
#####################################################
######1.9Enet+Lasso
#####################################################
for (alpha in seq(0.1,0.9,0.1)) {
  cv.fit = cv.glmnet(x = trainexp,
                     y = grouptrain,
                     family = "binomial", alpha = alpha, nfolds = 100)
  fit = glmnet(x = trainexp,
               y = grouptrain,
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(data0)
  Enetgenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
  Enetgenes=Enetgenes[2:length(Enetgenes)]
  trainexp2=as.data.frame(trainexp[,Enetgenes])
  testexp2=as.data.frame(testexp[,Enetgenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test2=data.frame(group=grouptest,testexp2)
  datalist2=list(Train=train2,Test=test2)
cv.fit = cv.glmnet(x = trainexp2,
                   y = grouptrain,
                   family = "binomial", alpha = 1, nfolds = 1000)
fit = glmnet(x = trainexp2,
             y = grouptrain,
             family = "binomial", alpha = 1, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
Lassogenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Lassogenes=Lassogenes[2:length(Lassogenes)]
modelname=paste0('Enet','[alpha=',alpha,']+Lasso')
modelgenes[[modelname]]=Lassogenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
}

#####################################################
######2.glm
#####################################################
fit <- step(glm(formula = ifelse(grouptrain=="con",0,1) ~ .,
                family = "binomial", 
                data = as.data.frame(trainexp)),trace = 0)
fit$subFeature = colnames(trainexp)
glmgenes=names(coef(fit))[2:length(names(coef(fit)))]
modelname=paste0('glm')
modelgenes[[modelname]]=glmgenes
rs <- lapply(datalist,function(x){cbind(x[,1:2],RS=as.data.frame(predict(fit, type = 'response', as.data.frame(x)))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
trainexp2=as.data.frame(trainexp[,glmgenes])
testexp2=as.data.frame(testexp[,glmgenes])
train2=data.frame(group=grouptrain,trainexp2)
test2=data.frame(group=grouptest,testexp2)
datalist2=list(Train=train2,Test=test2)

#####################################################
######2.2.glm+SVM
#####################################################
data <- as.data.frame(trainexp2)
cl <- makeCluster(16)
registerDoParallel(cl)
fit=rfe(x=data,
        y=as.numeric(as.factor(grouptrain)),
        rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
        methods="svmRadial")
stopCluster(cl)
modelname=paste0("glm+svm")
modelgenes[[modelname]]=fit[["optVariables"]]
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(fit, as.data.frame(x)))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)

#####################################################
######2.3.glm+plsRglmmodel
#####################################################
cv.plsRglm.res = cv.plsRglm(formula = grouptrain ~ ., 
                            data = as.data.frame(trainexp2),
                            nt=10, verbose = FALSE)
fit <- plsRglm(ifelse(grouptrain=="con",0,1), 
               as.data.frame(trainexp2), 
               modele = "pls-glm-logistic",
               verbose = F, sparse = T)
fit$subFeature = colnames(trainexp2)
plsRglmgenes=rownames(fit$Coeffs)[fit$Coeffs!=0]
modelname=paste0('glm+plsRglmmodel')
modelgenes[[modelname]]=plsRglmgenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(fit, type = "response", as.data.frame(x[,2:ncol(x)])))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)

#####################################################
######2.4.glm+RF
#####################################################
rf=randomForest(as.factor(grouptrain)~., data=trainexp2, ntree=500)
optionTrees=which.min(rf$err.rate[,1])
rf2=randomForest(as.factor(grouptrain)~., data=trainexp2, ntree=optionTrees)
importance=importance(x=rf2)
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes=names(rfGenes[rfGenes>0]) 
modelname=paste0('glm+RF')
modelgenes[[modelname]]=rfGenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(rf2, type = "response", as.data.frame(x[,2:ncol(x)])))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)

#####################################################
######2.5.glm+rpart
#####################################################
mod1<-rpart(as.factor(grouptrain)~.,data = trainexp2,method = "class")
importances <- varImp(mod1)
jcsimportances=as.matrix(importances %>%
                           arrange(desc(Overall)))
JCSGenes=jcsimportances[order(jcsimportances[,"Overall"], decreasing = TRUE),]
JCSGenes=names(JCSGenes[JCSGenes>0])
modelname=paste0('glm+rpart')
modelgenes[[modelname]]=JCSGenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(mod1, as.data.frame(x[,2:ncol(x)])))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)

#####################################################
######2.6.glm+GBM
#####################################################
fit <- gbm(formula = ifelse(grouptrain=="con",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 6)
best <- which.min(fit$cv.error)
fit <- gbm(formula = ifelse(grouptrain=="con",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001, n.cores = 8)
fit$subFeature = colnames(trainexp2)
GBMgenes=rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf>0]
modelname=paste0('glm+GBM')
modelgenes[[modelname]]=GBMgenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(fit, type = "response", as.data.frame(x)))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
#####################################################
######2.7.glm+Ridge
#####################################################
cv.fit = cv.glmnet(x = trainexp2,
                   y = grouptrain,
                   family = "binomial", alpha = 0, nfolds = 1000)
fit = glmnet(x = trainexp2,
             y = grouptrain,
             family = "binomial", alpha = 0, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
Ridgegenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Ridgegenes=Ridgegenes[2:length(Ridgegenes)]
modelname=paste0('glm+Ridge')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)

#####################################################
######2.8.glm+Lasso
#####################################################
cv.fit = cv.glmnet(x = trainexp2,
                   y = grouptrain,
                   family = "binomial", alpha = 1, nfolds = 1000)
fit = glmnet(x = trainexp2,
             y = grouptrain,
             family = "binomial", alpha = 1, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
lassogenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
lassogenes=lassogenes[2:length(lassogenes)]
modelname=paste0('glm+Lasso')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)

#####################################################
######3.SVM
#####################################################
data <- as.data.frame(trainexp)
cl <- makeCluster(16)
registerDoParallel(cl)
fit=rfe(x=data,
        y=as.numeric(as.factor(grouptrain)),
        rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
        methods="svmRadial")
stopCluster(cl)
modelname=paste0("glm+svm")
svmgenes=fit[["optVariables"]]
modelgenes[[modelname]]=svmgenes
rs <- lapply(datalist,function(x){cbind(x[,1:2],RS=as.data.frame(predict(fit, as.data.frame(x)))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
trainexp2=as.data.frame(trainexp[,svmgenes])
testexp2=as.data.frame(testexp[,svmgenes])
train2=data.frame(group=grouptrain,trainexp2)
test2=data.frame(group=grouptest,testexp2)
datalist2=list(Train=train2,Test=test2)

#####################################################
######3.2.SVM+plsRglmmodel 欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
cv.plsRglm.res = cv.plsRglm(formula = grouptrain ~ ., 
                            data = as.data.frame(trainexp2),
                            nt=10, verbose = FALSE)
fit <- plsRglm(ifelse(grouptrain=="con",0,1), 
               as.data.frame(trainexp2), 
               modele = "pls-glm-logistic",
               verbose = F, sparse = T)
fit$subFeature = colnames(trainexp2)
plsRglmgenes=rownames(fit$Coeffs)[fit$Coeffs!=0]
modelname=paste0('SVM+plsRglmmodel')
modelgenes[[modelname]]=plsRglmgenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(fit, type = "response", as.data.frame(x[,2:ncol(x)])))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
#####################################################
######3.3.SVM+RF 欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
rf=randomForest(as.factor(grouptrain)~., data=trainexp2, ntree=500)
optionTrees=which.min(rf$err.rate[,1])
rf2=randomForest(as.factor(grouptrain)~., data=trainexp2, ntree=optionTrees)
importance=importance(x=rf2)
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes=names(rfGenes[rfGenes>0]) 
modelname=paste0('SVM+RF')
modelgenes[[modelname]]=rfGenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(rf2, type = "response", as.data.frame(x[,2:ncol(x)])))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
#####################################################
######3.4.SVM+rpart 欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
mod1<-rpart(as.factor(grouptrain)~.,data = trainexp2,method = "class")
importances <- varImp(mod1)
jcsimportances=as.matrix(importances %>%
                           arrange(desc(Overall)))
JCSGenes=jcsimportances[order(jcsimportances[,"Overall"], decreasing = TRUE),]
JCSGenes=names(JCSGenes[JCSGenes>0])
modelname=paste0('SVM+rpart')
modelgenes[[modelname]]=JCSGenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(mod1, as.data.frame(x[,2:ncol(x)])))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
#####################################################
######3.5.SVM+GBM 欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
fit <- gbm(formula = ifelse(grouptrain=="con",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 6)
best <- which.min(fit$cv.error)
fit <- gbm(formula = ifelse(grouptrain=="con",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001, n.cores = 8)
fit$subFeature = colnames(trainexp2)
GBMgenes=rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf>0]
modelname=paste0('SVM+GBM')
modelgenes[[modelname]]=GBMgenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(fit, type = "response", as.data.frame(x)))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
#####################################################
######3.6.SVM+Ridge 欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
cv.fit = cv.glmnet(x = trainexp2,
                   y = grouptrain,
                   family = "binomial", alpha = 0, nfolds = 1000)
fit = glmnet(x = trainexp2,
             y = grouptrain,
             family = "binomial", alpha = 0, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
Ridgegenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Ridgegenes=Ridgegenes[2:length(Ridgegenes)]
modelname=paste0('SVM+Ridge')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
#####################################################
######3.6.SVM+Lasso 欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
cv.fit = cv.glmnet(x = trainexp2,
                   y = grouptrain,
                   family = "binomial", alpha = 1, nfolds = 1000)
fit = glmnet(x = trainexp2,
             y = grouptrain,
             family = "binomial", alpha = 1, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
Lassogenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Lassogenes=Enetgenes[2:length(Lassogenes)]
modelname=paste0('SVM+Lasso')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)

#####################################################
######4.plsRglmmodel 欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
cv.plsRglm.res = cv.plsRglm(formula = grouptrain ~ ., 
                            data = as.data.frame(trainexp),
                            nt=10, verbose = FALSE)
fit <- plsRglm(ifelse(grouptrain=="con",0,1), 
               as.data.frame(trainexp), 
               modele = "pls-glm-logistic",
               verbose = F, sparse = T)
fit$subFeature = colnames(trainexp)
plsRglmgenes=rownames(fit$Coeffs)[fit$Coeffs!=0][-1]
modelname=paste0('plsRglmmodel')
modelgenes[[modelname]]=plsRglmgenes
rs <- lapply(datalist,function(x){cbind(x[,1:2],RS=as.data.frame(predict(fit, type = "response", as.data.frame(x[,2:ncol(x)])))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
trainexp2=as.data.frame(trainexp[,plsRglmgenes])
testexp2=as.data.frame(testexp[,plsRglmgenes])
train2=data.frame(group=grouptrain,trainexp2)
test2=data.frame(group=grouptest,testexp2)
datalist2=list(Train=train2,Test=test2)
#####################################################
######4.2.plsRglmmodel+RF 欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
rf=randomForest(as.factor(grouptrain)~., data=trainexp2, ntree=500)
optionTrees=which.min(rf$err.rate[,1])
rf2=randomForest(as.factor(grouptrain)~., data=trainexp2, ntree=optionTrees)
importance=importance(x=rf2)
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes=names(rfGenes[rfGenes>0]) 
modelname=paste0('plsRglmmodel+RF')
modelgenes[[modelname]]=rfGenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(rf2, type = "response", as.data.frame(x[,2:ncol(x)])))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)

#####################################################
######4.3.plsRglmmodel+rpart 欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
mod1<-rpart(as.factor(grouptrain)~.,data = trainexp2,method = "class")
importances <- varImp(mod1)
jcsimportances=as.matrix(importances %>%
                           arrange(desc(Overall)))
JCSGenes=jcsimportances[order(jcsimportances[,"Overall"], decreasing = TRUE),]
JCSGenes=names(JCSGenes[JCSGenes>0])
modelname=paste0('plsRglmmodel+rpart')
modelgenes[[modelname]]=JCSGenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(mod1, as.data.frame(x[,2:ncol(x)])))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
#####################################################
######4.4.plsRglmmodel+GBM 欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
fit <- gbm(formula = ifelse(grouptrain=="con",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 6)
best <- which.min(fit$cv.error)
fit <- gbm(formula = ifelse(grouptrain=="con",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001, n.cores = 8)
fit$subFeature = colnames(trainexp2)
GBMgenes=rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf>0]
modelname=paste0('plsRglmmodel+GBM')
modelgenes[[modelname]]=GBMgenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(fit, type = "response", as.data.frame(x)))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)

#####################################################
######4.5.plsRglmmodel+Ridge 欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
cv.fit = cv.glmnet(x = trainexp2,
                   y = grouptrain,
                   family = "binomial", alpha = 0, nfolds = 1000)
fit = glmnet(x = trainexp2,
             y = grouptrain,
             family = "binomial", alpha = 0, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
Ridgegenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Ridgegenes=Ridgegenes[2:length(Ridgegenes)]
modelname=paste0('plsRglmmodel+Ridge')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
#####################################################
######4.6.plsRglmmodel+Lasso 欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
cv.fit = cv.glmnet(x = trainexp2,
                   y = grouptrain,
                   family = "binomial", alpha = 1, nfolds = 1000)
fit = glmnet(x = trainexp2,
             y = grouptrain,
             family = "binomial", alpha = 1, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
Lassogenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Lassogenes=Lassogenes[2:length(Lassogenes)]
modelname=paste0('plsRglmmode+Lasso')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)

#####################################################
######5.RF 欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
rf=randomForest(as.factor(grouptrain)~., data=trainexp, ntree=500)
optionTrees=which.min(rf$err.rate[,1])
rf2=randomForest(as.factor(grouptrain)~., data=trainexp, ntree=optionTrees)
importance=importance(x=rf2)
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes=names(rfGenes[rfGenes>0]) 
modelname=paste0('RF')
modelgenes[[modelname]]=rfGenes
rs <- lapply(datalist,function(x){cbind(x[,1:2],RS=as.data.frame(predict(rf2, type = "response", as.data.frame(x[,2:ncol(x)])))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
trainexp2=as.data.frame(trainexp[,rfGenes])
testexp2=as.data.frame(testexp[,rfGenes])
train2=data.frame(group=grouptrain,trainexp2)
test2=data.frame(group=grouptest,testexp2)
datalist2=list(Train=train2,Test=test2)

#####################################################
######5.2.RF+rpart 欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
mod1<-rpart(as.factor(grouptrain)~.,data = trainexp2,method = "class")
importances <- varImp(mod1)
jcsimportances=as.matrix(importances %>%
                           arrange(desc(Overall)))
JCSGenes=jcsimportances[order(jcsimportances[,"Overall"], decreasing = TRUE),]
JCSGenes=names(JCSGenes[JCSGenes>0])
modelname=paste0('RF+rpart')
modelgenes[[modelname]]=JCSGenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(mod1, as.data.frame(x[,2:ncol(x)])))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
#####################################################
######5.3.RF+GBM 欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
fit <- gbm(formula = ifelse(grouptrain=="con",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 6)
best <- which.min(fit$cv.error)
fit <- gbm(formula = ifelse(grouptrain=="con",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001, n.cores = 8)
fit$subFeature = colnames(trainexp2)
GBMgenes=rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf>0]
modelname=paste0('RF+GBM')
modelgenes[[modelname]]=GBMgenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(fit, type = "response", as.data.frame(x)))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
#####################################################
######5.4.RF+Ridge 欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
cv.fit = cv.glmnet(x = trainexp2,
                   y = grouptrain,
                   family = "binomial", alpha = 0, nfolds = 1000)
fit = glmnet(x = trainexp2,
             y = grouptrain,
             family = "binomial", alpha = 0, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
Ridgegenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Ridgegenes=Ridgegenes[2:length(Ridgegenes)]
modelname=paste0('RF+Ridge')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
#####################################################
######5.5.RF+Lasso 欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
cv.fit = cv.glmnet(x = trainexp2,
                   y = grouptrain,
                   family = "binomial", alpha = 1, nfolds = 1000)
fit = glmnet(x = trainexp2,
             y = grouptrain,
             family = "binomial", alpha = 1, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
Ridgegenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Ridgegenes=Ridgegenes[2:length(Ridgegenes)]
modelname=paste0('RF+Lasso')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
#####################################################
######6.rpart 欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
mod1<-rpart(as.factor(grouptrain)~.,data = as.data.frame(trainexp),method = "class")
importances <- varImp(mod1)
jcsimportances=as.matrix(importances %>%
                           arrange(desc(Overall)))
JCSGenes=jcsimportances[order(jcsimportances[,"Overall"], decreasing = TRUE),]
JCSGenes=names(JCSGenes[JCSGenes>0])
modelname=paste0('rpart')
modelgenes[[modelname]]=JCSGenes
rs <- lapply(datalist,function(x){cbind(x[,1:2],RS=as.data.frame(predict(mod1, as.data.frame(x[,2:ncol(x)])))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
trainexp2=as.data.frame(trainexp[,JCSGenes])
testexp2=as.data.frame(testexp[,JCSGenes])
train2=data.frame(group=grouptrain,trainexp2)
test2=data.frame(group=grouptest,testexp2)
datalist2=list(Train=train2,Test=test2)
#####################################################
######6.2.rpart+GBM欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
fit <- gbm(formula = ifelse(grouptrain=="con",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 6)
best <- which.min(fit$cv.error)
fit <- gbm(formula = ifelse(grouptrain=="con",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001, n.cores = 8)
fit$subFeature = colnames(trainexp2)
GBMgenes=rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf>0]
modelname=paste0('rpart+GBM')
modelgenes[[modelname]]=GBMgenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=as.data.frame(predict(fit, type = "response", as.data.frame(x)))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
#####################################################
######6.3.rpart+Ridge欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
cv.fit = cv.glmnet(x = trainexp2,
                   y = grouptrain,
                   family = "binomial", alpha = 0, nfolds = 1000)
fit = glmnet(x = trainexp2,
             y = grouptrain,
             family = "binomial", alpha = 0, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
Ridgegenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Ridgegenes=Enetgenes[2:length(Ridgegenes)]
modelname=paste0('rpart+Ridge')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
#####################################################
######6.4.rpart+Lasso欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
cv.fit = cv.glmnet(x = trainexp2,
                   y = grouptrain,
                   family = "binomial", alpha = 1, nfolds = 1000)
fit = glmnet(x = trainexp2,
             y = grouptrain,
             family = "binomial", alpha = 1, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
Lassogenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Lassogenes=Lassogenes[2:length(Lassogenes)]
modelname=paste0('rpart+Lasso')
modelgenes[[modelname]]=Lassogenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
#####################################################
######7.GBM欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
fit <- gbm(formula = ifelse(grouptrain=="con",0,1) ~ .,
           data = as.data.frame(trainexp),
           distribution = 'bernoulli',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 6)
best <- which.min(fit$cv.error)
fit <- gbm(formula = ifelse(grouptrain=="con",0,1) ~ .,
           data = as.data.frame(trainexp),
           distribution = 'bernoulli',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001, n.cores = 8)
fit$subFeature = colnames(trainexp2)
GBMgenes=rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf>0]
modelname=paste0('GBM')
modelgenes[[modelname]]=GBMgenes
rs <- lapply(datalist,function(x){cbind(x[,1:2],RS=as.data.frame(predict(fit, type = "response", as.data.frame(x)))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
trainexp2=as.data.frame(trainexp[,GBMgenes])
testexp2=as.data.frame(testexp[,GBMgenes])
train2=data.frame(group=grouptrain,trainexp2)
test2=data.frame(group=grouptest,testexp2)
datalist2=list(Train=train2,Test=test2)
#####################################################
######7.2.GBM+Ridge欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
cv.fit = cv.glmnet(x = trainexp2,
                   y = grouptrain,
                   family = "binomial", alpha = 0, nfolds = 1000)
fit = glmnet(x = trainexp2,
             y = grouptrain,
             family = "binomial", alpha = 0, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
Ridgegenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Ridgegenes=Ridgegenes[2:length(Ridgegenes)]
modelname=paste0('GBM+Ridge')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
#####################################################
######7.3.GBM+Lasso欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
cv.fit = cv.glmnet(x = trainexp2,
                   y = grouptrain,
                   family = "binomial", alpha = 1, nfolds = 1000)
fit = glmnet(x = trainexp2,
             y = grouptrain,
             family = "binomial", alpha = 1, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
Ridgegenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Ridgegenes=Enetgenes[2:length(Ridgegenes)]
modelname=paste0('GBM+Lasso')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
#####################################################
######8.Ridge欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
cv.fit = cv.glmnet(x = trainexp,
                   y = grouptrain,
                   family = "binomial", alpha = 0, nfolds = 1000)
fit = glmnet(x = trainexp,
             y = grouptrain,
             family = "binomial", alpha = 0, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
Ridgegenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Ridgegenes=Ridgegenes[2:length(Ridgegenes)]
modelname=paste0('Ridge')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist,function(x){cbind(x[,1:2],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
trainexp2=as.data.frame(trainexp[,Ridgegenes])
testexp2=as.data.frame(testexp[,Ridgegenes])
train2=data.frame(group=grouptrain,trainexp2)
test2=data.frame(group=grouptest,testexp2)
datalist2=list(Train=train2,Test=test2)
#####################################################
######8.2.Ridge+Lasso欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
cv.fit = cv.glmnet(x = trainexp2,
                   y = grouptrain,
                   family = "binomial", alpha = 1, nfolds = 1000)
fit = glmnet(x = trainexp2,
             y = grouptrain,
             family = "binomial", alpha = 1, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
Ridgegenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Ridgegenes=Ridgegenes[2:length(Ridgegenes)]
modelname=paste0('Ridge+Lasso')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist2,function(x){cbind(x[,1:2],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)
#####################################################
######9.Lasso欢迎关注《叉叉滴同学的生信笔记》######
#####################################################
cv.fit = cv.glmnet(x = trainexp,
                   y = grouptrain,
                   family = "binomial", alpha = 1, nfolds = 1000)
fit = glmnet(x = trainexp,
             y = grouptrain,
             family = "binomial", alpha = 1, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
Lassogenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Lassogenes=Lassogenes[2:length(Lassogenes)]
modelname=paste0('Lasso')
modelgenes[[modelname]]=Lassogenes
rs <- lapply(datalist,function(x){cbind(x[,1:2],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')
for(i in 1:length(rs)){
  rocdata=rs[[i]]
  roc1=roc(rocdata$group, as.numeric(rocdata$RS))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("./ROC/",modelname,"_",names(rs)[i],"_ROC.pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(modelname,"_",names(rs)[i]))
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
AUCs$Model <- modelname
result <- rbind(result,AUCs)


######绘制AUC热图######
result$auc=round(result$auc,4)
result2 <- result 
result2=setDT(result2)  # 将数据框转换为data.table  
result2_train=result2[result2$ID=="Train",]
result2_train=as.data.frame(result2_train)
result2_test=result2[result2$ID=="Test",]
result2_test=as.data.frame(result2_test)
aucnums=data.frame(Modle=result2_train$Model,Train=result2_train$auc,Test=result2_test$auc)
aucnums[,-1] <- apply(aucnums[,-1], 2, as.numeric)
aucnums$All <- apply(aucnums[,2:3], 1, mean)
# 根据auc排序
aucnums <- aucnums[order(aucnums$All, decreasing = T),]
#输出C指数结果
write.table(aucnums,"out_AUC.txt", col.names = T, row.names = F, sep = "\t", quote = F)
nums <- aucnums[, 2:3]%>%as.matrix()
rownames(nums)=aucnums$Modle
##热图绘制
auc_mat=nums
# 计算每种算法在所有队列中平均C-index
avg_auc <- apply(auc_mat, 1, mean)     
# 对各算法C-index由高到低排序
avg_auc <- sort(avg_auc, decreasing = T)    
# 对C-index矩阵排序
auc_mat <- auc_mat[names(avg_auc), ]      
# 保留三位小数
avg_auc <- as.numeric(format(avg_auc, digits = 3, nsmall = 3)) 
row_ha = rowAnnotation(bar = anno_barplot(avg_auc, bar_width = 0.8, border = FALSE,
                                          gp = gpar(fill = "steelblue", col = NA),
                                          add_numbers = T, numbers_offset = unit(-10, "mm"),
                                          axis_param = list("labels_rot" = 0),
                                          numbers_gp = gpar(fontsize = 9, col = "white"),
                                          width = unit(3, "cm")),
                       show_annotation_name = F)
CohortCol <- c('#BD3C29','#0172B6')
names(CohortCol) <- colnames(auc_mat)
col_ha = columnAnnotation("Cohort" = colnames(auc_mat),
                          col = list("Cohort" = CohortCol),
                          show_annotation_name = F)

cellwidth = 1
cellheight = 0.5
hm <- Heatmap(as.matrix(auc_mat), name = "AUC",
              right_annotation = row_ha, 
              top_annotation = col_ha,
              col = c("#4195C1", "#FFFFFF", "#CB5746"), 
              rect_gp = gpar(col = "black", lwd = 1), # 边框设置为黑色
              cluster_columns = FALSE, cluster_rows = FALSE, # 不进行聚类，无意义
              show_column_names = FALSE, 
              show_row_names = TRUE,
              row_names_side = "left",
              width = unit(cellwidth * ncol(auc_mat) + 2, "cm"),
              height = unit(cellheight * nrow(auc_mat), "cm"),
              column_split = factor(colnames(auc_mat), levels = colnames(auc_mat)), 
              column_title = NULL,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(label = format(auc_mat[i, j], digits = 3, nsmall = 3),
                          x, y, gp = gpar(fontsize = 10))
              }
)

pdf(file.path( "auc.pdf"), width = cellwidth * ncol(auc_mat) + 5, height = cellheight * nrow(auc_mat) * 0.45)
draw(hm)
invisible(dev.off())

######多算法筛选特征基因######
#1.Enet
cv.fit = cv.glmnet(x = data0,
                   y = group0,
                   family = "binomial", alpha = 0.9, nfolds = 100)
fit = glmnet(x = data0,
             y = group0,
             family = "binomial", alpha = 0.9, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
Enetgenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Enetgenes=Enetgenes[2:length(Enetgenes)]

#2.glm
fit <- step(glm(formula = ifelse(group0=="con",0,1) ~ .,
                family = "binomial", 
                data = as.data.frame(data0)),trace = 0)
fit$subFeature = colnames(data0)
glmgenes=names(coef(fit))[2:length(names(coef(fit)))]

#3.SVM
data <- as.data.frame(data0)
cl <- makeCluster(16)
registerDoParallel(cl)
fit=rfe(x=data,
        y=as.numeric(as.factor(group0)),
        rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
        methods="svmRadial")
stopCluster(cl)
SVMgene=fit[["optVariables"]]

#4.plsRglmmodel
cv.plsRglm.res = cv.plsRglm(formula = group0 ~ ., 
                            data = as.data.frame(data0),
                            nt=10, verbose = FALSE)
fit <- plsRglm(ifelse(group0=="con",0,1), 
               as.data.frame(data0), 
               modele = "pls-glm-logistic",
               verbose = F, sparse = T)
fit$subFeature = colnames(data0)
plsRglmgenes=rownames(fit$Coeffs)[fit$Coeffs!=0]

#5.RF
rf=randomForest(as.factor(group0)~., data=data0, ntree=500)
optionTrees=which.min(rf$err.rate[,1])
rf2=randomForest(as.factor(group0)~., data=data0, ntree=optionTrees)
importance=importance(x=rf2)
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes=names(rfGenes[rfGenes>0]) 

#6.rpart
mod1<-rpart(as.factor(group0)~.,data = as.data.frame(data0),method = "class")
importances <- varImp(mod1)
jcsimportances=as.matrix(importances %>%
                           arrange(desc(Overall)))
JCSGenes=jcsimportances[order(jcsimportances[,"Overall"], decreasing = TRUE),]
JCSGenes=names(JCSGenes[JCSGenes>0])

#7,GBM
fit <- gbm(formula = ifelse(group0=="con",0,1) ~ .,
           data = as.data.frame(data0),
           distribution = 'bernoulli',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 6)
best <- which.min(fit$cv.error)
fit <- gbm(formula = ifelse(group0=="con",0,1) ~ .,
           data = as.data.frame(data0),
           distribution = 'bernoulli',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001, n.cores = 8)
fit$subFeature = colnames(data0)
GBMgenes=rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf>0]

#8.Ridge
cv.fit = cv.glmnet(x = data0,
                   y = group0,
                   family = "binomial", alpha = 0, nfolds = 1000)
fit = glmnet(x = data0,
             y = group0,
             family = "binomial", alpha = 0, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
Ridgegenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Ridgegenes=Ridgegenes[2:length(Ridgegenes)]

#9.Lasso
cv.fit = cv.glmnet(x = data0,
                   y = group0,
                   family = "binomial", alpha = 1, nfolds = 1000)
fit = glmnet(x = data0,
             y = group0,
             family = "binomial", alpha = 1, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(data0)
Lassogenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Lassogenes=Lassogenes[2:length(Lassogenes)]

######特征交集######
mutationsML=data.frame(ID=colnames(data0),
                       Enet=ifelse(colnames(data0)  %in% Enetgenes,1,0),
                       GLM=ifelse(colnames(data0)  %in% glmgenes,1,0),
                       SVM=ifelse(colnames(data0)  %in% SVMgene,1,0),
                       plsRglmmodel=ifelse(colnames(data0)  %in% plsRglmgenes,1,0),
                       RF=ifelse(colnames(data0)  %in% rfGenes,1,0),
                       rpart=ifelse(colnames(data0)  %in% JCSGenes,1,0),
                       GBM=ifelse(colnames(data0)  %in% GBMgenes,1,0),
                       Ridge=ifelse(colnames(data0)  %in% Ridgegenes,1,0),
                       Lasso=ifelse(colnames(data0)  %in% Lassogenes,1,0))
upset(mutationsML,sets.bar.color = "#FB8072" ,matrix.color = "#80B1D3",color.pal = "#80B1D3",main.bar.color = "#80B1D3",
      nsets =9, 
      order.by = "freq",#交集排列方式，freq为升序，degree为降序
      decreasing = T#变量排列方式，T降序排列
)
queries<-list(list(query=intersects,
                   params=list("Enet","GLM","plsRglmmodel","RF","rpart","GBM","Ridge","Lasso"),
                   active=T,
                   color="red")) 
pdf("upset.pdf")
upset(mutationsML, sets.bar.color = "#FB8072" ,matrix.color = "#80B1D3",color.pal = "#80B1D3",main.bar.color = "#80B1D3",
      nsets =9,
      order.by = "freq",#交集排列方式，freq为升序，degree为降序
      decreasing = T,#变量排列方式，T降序排列
      queries=queries
)
dev.off()

#####        欢迎关注《叉叉滴同学的生信笔记》    #####
#####                            #####
#####      唯一微信：sciwink     #####
#####      唯一微信：sciwink     #####
#####                            #####
#####        欢迎关注《叉叉滴同学的生信笔记》   #####