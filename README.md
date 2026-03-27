# WGCNA-ML-Drug-Screening: Integrated Metastasis Drug Discovery Pipeline 🧬💊

![R](https://img.shields.io/badge/Language-R-blue.svg) ![License](https://img.shields.io/badge/License-MIT-green.svg) ![Bioinformatics](https://img.shields.io/badge/Domain-Bioinformatics-orange.svg)

## 📖 Overview (项目概述)

This repository contains a comprehensive bioinformatics analysis pipeline designed to uncover key driver genes involved in **cancer metastasis** and identify potential therapeutic drugs.

The workflow integrates **WGCNA** (Weighted Gene Co-expression Network Analysis) for module detection, **LASSO** regression for feature selection, and multiple **Machine Learning** algorithms (Random Forest, SVM, XGBoost) to construct a robust diagnostic model. Furthermore, it constructs a **multi-level regulatory network** (TF-miRNA-mRNA) and performs **drug screening** with safety (ADMET) assessment.

> **Key Dataset**: GSE81558 (Training), GSE49355 (Validation)

## ⚙️ Workflow (分析流程)

The pipeline consists of 10 major steps:

1.  **Data Preprocessing**: Quality control and normalization using `limma`.
2.  **WGCNA Analysis**: Construction of co-expression networks to identify metastasis-related modules.
3.  **Functional Enrichment**: GO and KEGG pathway analysis for biological insights.
4.  **Data Validation**: Independent validation using an external dataset (GSE49355).
5.  **Feature Selection**: LASSO regression to screen for robust biomarkers.
6.  **Model Construction**: Building diagnostic models using **Random Forest**, **SVM**, and **XGBoost**.
7.  **Signature Visualization**: Heatmaps and expression plots of key genes.
8.  **Regulatory Network**: Construction of upstream TF and miRNA regulatory networks using TRRUST and multiMiR.
9.  **Drug Screening**: Multi-level drug identification (Direct targets, TF-targeting, miRNA therapeutics).
10. **Safety Assessment**: Drug safety evaluation using SIDER and ADMET prediction (Lipinski's Rule).

## 🛠️ Prerequisites (环境依赖)

Ensure you have R (>= 4.0.0) installed. The following major packages are required:

```r
install.packages(c("WGCNA", "limma", "tidyverse", "glmnet", "randomForest", "e1071", "xgboost", "pROC", "pheatmap", "igraph", "enrichR", "webchem"))
# For Bioconductor packages:
BiocManager::install(c("GEOquery", "clusterProfiler", "org.Hs.eg.db", "multiMiR"))
