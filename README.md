# My NSCLC Project

## Overview
This repository hosts R scripts and supporting artifacts for a non-small cell lung cancer (NSCLC, LUAD-focused) research workflow. The scripts cover GEO/TCGA preprocessing, differential analysis, enrichment, immune infiltration, survival modeling, and downstream visualization.

## Key R Scripts
- GEO microarray processing: [01.GEO芯片处理.R](01.GEO芯片处理.R)
- TCGA-LUAD processing: [02.TCGA_LUAD.R](02.TCGA_LUAD.R)
- Common genes and enrichment: [03.Common Genes & Enrichment analysis.R](03.Common%20Genes%20&%20Enrichment%20analysis.R)
- Immune infiltration analysis: [04.Immune_Infiltration.R](04.Immune_Infiltration.R)
- Cox and survival modeling: [COX.R](COX.R) and [new_plot/COX regression.R](new_plot/COX%20regression.R)
- GEO clinical metadata: [GEO_clinical_data.R](GEO_clinical_data.R)
- GO/KEGG summarization: [GOKEGG_comprehensive.R](GOKEGG_comprehensive.R)
- Ferritin-related analysis: [Ferr_differ_in_TCGA.R](Ferr_differ_in_TCGA.R)
- Binary classification model: [二分类模型.R](二分类模型.R)
- Integrated GSE189357 workflow: [GSE189357/NSCLC_all in one.R](GSE189357/NSCLC_all%20in%20one.R)
- GSE189357 visualization: [GSE189357/UMAP_plot.R](GSE189357/UMAP_plot.R)

## Repository Layout
- [Mime-main/](Mime-main/): Packaged methods and utilities used across analyses.
- [new_plot/](new_plot/): Plotting scripts and outputs for figures and tables.
- [Plot/](Plot/): Additional plotting outputs and legacy figures.
- [Date_13.08.25/](Date_13.08.25/): Time-stamped snapshots and related plots.
- [GSE189357/](GSE189357/): Dataset-specific analysis scripts and results.

## Notes
- Scripts assume local data files present in this workspace; adjust file paths as needed.
- Add package versioning (e.g., `renv`) if strict reproducibility is required.

## Dependencies (Detected in Scripts)
Core packages commonly used across workflows:
- Seurat, SeuratObject, SingleCellExperiment, monocle, monocle3, slingshot
- CellChat, GSVA, GSEABase, clusterProfiler, org.Hs.eg.db, DOSE, ReactomePA
- limma, edgeR, DESeq2, survival, survminer, timeROC
- dplyr, tidyr, tidyverse, data.table, readr, readxl, tibble
- ggplot2, patchwork, cowplot, ggpubr, ggforce, ggrepel, ggsci

Extended packages referenced in the scripts:
- aorsf, aplot, assertthat, BART, beepr, BiocManager, BiocParallel, bitops, broom
- cancerclass, caret, CCA, celldex, circlize, clustree, compareC, ComplexHeatmap, corrplot
- COSG, CoxBoost, DALEX, depen, doParallel, easyTCGA, enrichplot, factoextra, fastAdaboost
- fastSave, forestplot, forestploter, future, future.apply, gbm, GEOquery, ggbreak
- ggcorrplot, ggDCA, ggpol, ggsignif, glmnet, GOSemSim, gplots, grid, gridExtra
- harmony, Hmisc, igraph, immunedeconv, IOBR, irGSEA, kableExtra, kernlab, lme4
- lmerTest, magrittr, MASS, Matrix, matrixStats, mboost, meta, Metrics, MIME, Mime, Mime1
- miscTools, mixOmics, mixtools, Nebulosa, neuralnet, NeuralNetTools, obliqueRSF
- ontologyIndex, openxlsx, parallel, party, partykit, pec, pheatmap, plotrix, plsRcox
- plsRglm, plyr, ppcor, preprocessCore, presto, pROC, qs, randomcoloR, randomForest
- randomForestSRC, rattle, RColorBrewer, Rcpp, RCurl, recipes, regplot, remotes, reshape2
- rms, ROCit, ROCR, rpart, rpart.plot, scales, scater, scattermore, scatterplot3d
- scde, SCP, SCpubr, seqinr, SingleR, snowfall, sparkline, sparrow, starTracer, stats
- stringi, stringr, superpc, survivalROC, survivalsvm, sva, testthat, tradeSeq, tricycle
- UCell, UpSetR, venn, VennDiagram, VIM, viridis, visNetwork, xgboost

## License
Specify the license for this repository here.

## Contact
For questions, please contact the repository owner.
