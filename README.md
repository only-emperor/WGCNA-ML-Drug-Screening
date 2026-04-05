# Integrated Multi-Omics & AI-Driven Drug Discovery Pipeline for Metastasis 🧬💊

[![R](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)
[![Bioinformatics](https://img.shields.io/badge/Field-Bioinformatics-green.svg)](https://en.wikipedia.org/wiki/Bioinformatics)
[![AIDD](https://img.shields.io/badge/Focus-AI_Drug_Discovery-orange.svg)](https://en.wikipedia.org/wiki/Drug_discovery)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 📖 Overview

This repository contains a comprehensive **computational drug repurposing framework**. The pipeline integrates **WGCNA (Weighted Gene Co-expression Network Analysis)** with **Machine Learning algorithms** to identify key driver genes involved in Colorectal Cancer (CRC) metastasis and screen for potential therapeutic agents.

The workflow moves beyond simple differential expression by constructing a multi-level regulatory network (TF-miRNA-mRNA) and performing **in silico drug screening** with rigorous safety (ADMET) assessments.

**Key Datasets:**
*   **Training:** GSE81558
*   **Validation:** GSE49355

---

## ⚙️ Workflow & Methodology

The analysis pipeline consists of 10 major modules, designed to bridge the gap between omics data and clinical application:

### Phase 1: Biomarker Identification (Omics Integration)
1.  **Data Preprocessing**: Quality control, background correction, and normalization using `limma`.
2.  **WGCNA Analysis**: Construction of scale-free co-expression networks to identify metastasis-specific gene modules.
3.  **Functional Enrichment**: GO and KEGG pathway analysis using `clusterProfiler` for biological interpretation.
4.  **Robust Feature Selection**: Applying **LASSO regression** to screen for the most predictive biomarkers.

### Phase 2: Predictive Modeling (Machine Learning)
5.  **Model Construction**: Benchmarking multiple ML algorithms to build a diagnostic classifier:
    *   Random Forest (RF)
    *   Support Vector Machines (SVM)
    *   XGBoost
6.  **Validation**: Evaluated model performance (AUC) on independent external datasets (GSE49355).

### Phase 3: AIDD & Regulatory Network
7.  **Regulatory Network Construction**: Building upstream TF and miRNA regulatory networks using **TRRUST** and **multiMiR**.
8.  **Drug Screening (Repurposing)**: Multi-level identification of drug candidates:
    *   Direct target interactions
    *   TF-targeting drugs
    *   miRNA-based therapeutics
9.  **Safety Assessment (ADMET)**: Evaluating drug safety profiles (Absorption, Distribution, Metabolism, Excretion, and Toxicity) using **SIDER** and **Lipinski's Rule of Five**.

---

## 📊 Key Features

- **Automated Pipeline**: Scripts are modularized for reproducibility.
- **Deep Learning Ready**: The processed features from this pipeline can serve as high-quality inputs for Graph Neural Networks (GNNs).
- **Visualization**: Includes scripts for generating:
    - WGCNA Module Heatmaps
    - ROC Curves for Model Evaluation
    - Drug-Target Interaction Networks

---

## 🛠️ Prerequisites

Ensure you have **R (>= 4.0.0)** installed.

### Required R Packages
You can install the dependencies using the following commands:

```r
# CRAN Packages
install.packages(c("WGCNA", "tidyverse", "glmnet", "randomForest", 
                   "e1071", "xgboost", "pROC", "pheatmap", "igraph", "enrichR"))

# Bioconductor Packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("limma", "GEOquery", "clusterProfiler", 
                       "org.Hs.eg.db", "multiMiR"))
## 🚧 Project Status / TODO

Currently, this repository is under active development. **Please note that some underlying code and modular scripts still need to be uploaded and retained in the main branch.** Upcoming updates will include:
- [ ] Full code for Phase 3 ADMET screening.
- [ ] Automated bash scripts for the entire pipeline.
