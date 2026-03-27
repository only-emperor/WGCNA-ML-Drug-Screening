# ======================================================
# =========          第一部分: 差异分析         =========
# ======================================================
cat("\n")
cat("========================================\n")
cat("  第一部分: 差异表达分析开始\n")
cat("========================================\n\n")

# 1. 加载包
library(GEOquery)
library(tidyverse)
library(limma)
library(stringr)
library(dplyr)
library(WGCNA)

# 📌 设置主工作目录（训练集目录）
main_dir <- "D:/dow/GSE81558"
setwd(main_dir)

gset <- getGEO("GSE81558", destdir=".", AnnotGPL = F, getGPL = F)
exp <- exprs(gset[[1]])
exp <- normalizeBetweenArrays(exp)

# 加载注释文件
gp1 <- read.table("GPL15207-17536.txt",
                  header = TRUE, fill = TRUE, sep = "\t",
                  quote = "", comment.char = "#", stringsAsFactors = FALSE)

ids <- gp1[, c("ID", "Gene.Symbol")]
colnames(ids) <- c("probe_id", "symbol")

# 注释清理
ids <- ids[ids$symbol != "" & !is.na(ids$symbol), ]
ids$symbol <- str_split(ids$symbol, " /// ", simplify = TRUE)[, 1]

# 匹配探针与表达矩阵
exp_filtered <- exp[rownames(exp) %in% ids$probe_id, ]
ids_filtered <- ids[ids$probe_id %in% rownames(exp_filtered), ]
ids_filtered <- ids_filtered[match(rownames(exp_filtered), ids_filtered$probe_id), ]

# 多探针-单基因
exp_anno <- as.data.frame(exp_filtered)
exp_anno$gene <- ids_filtered$symbol

exp_gene <- exp_anno %>%
  group_by(gene) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

exp_gene_mat <- as.matrix(exp_gene %>% select(-gene))
rownames(exp_gene_mat) <- exp_gene$gene

# 样本信息处理
sample_info <- pData(gset[[1]])
sample_info$patient_id <- str_extract(sample_info$title, "(F?\\d{3,4})")
sample_info$patient_id <- gsub("^F", "", sample_info$patient_id)

sample_info$group <- case_when(
  grepl("Primary", sample_info$title, ignore.case = TRUE) ~ "Primary",
  grepl("Metastasis", sample_info$title, ignore.case = TRUE) ~ "Metastasis",
  grepl("Non-tumoral", sample_info$title, ignore.case = TRUE) ~ "Normal",
  TRUE ~ "Other"
)

paired_patients <- sample_info %>%
  group_by(patient_id) %>%
  summarise(groups = paste(unique(group), collapse=",")) %>%
  filter(grepl("Primary", groups) & grepl("Metastasis", groups)) %>%
  pull(patient_id)

paired_samples <- sample_info %>%
  filter(patient_id %in% paired_patients) %>%
  arrange(patient_id, group) %>%
  dplyr::select(geo_accession, title, patient_id, group)

paired_samples2 <- paired_samples %>%
  filter(group %in% c("Primary","Metastasis"))

exp_paired <- exp_gene_mat[, paired_samples2$geo_accession]

# 差异分析
group <- factor(paired_samples2$group, levels = c("Primary","Metastasis"))
patient <- factor(paired_samples2$patient_id)
design <- model.matrix(~ patient + group)

fit <- lmFit(exp_paired, design)
fit2 <- eBayes(fit)
results <- topTable(fit2, coef = "groupMetastasis", n = Inf)

deg <- results %>%
  filter(abs(logFC) > 0.585, adj.P.Val < 0.05)

write.csv(results, "DEG_all_results.csv", row.names = TRUE)
write.csv(deg, "DEG_filtered_results.csv", row.names = TRUE)

cat("✅ 差异分析完成！\n")
cat(paste0("   总差异基因数: ", nrow(deg), "\n\n"))

# ======================================================
# =========      第二部分: WGCNA     =========
# ======================================================

cat("\n")
cat("========================================\n")
cat("  第二部分: WGCNA 分析开始\n")
cat("========================================\n\n")

# 📌 创建输出目录
wgcna_dir <- file.path(main_dir, "WGCNA_output")
if (!dir.exists(wgcna_dir)) {
  dir.create(wgcna_dir)
}

data <- exp_paired[apply(exp_paired, 1, sd) > 0.5, ]
datExpr0 <- t(data)

gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

traitData <- paired_samples2 %>%
  mutate(
    Primary = ifelse(group == "Primary", 1, 0),
    Metastasis = ifelse(group == "Metastasis", 1, 0)
  ) %>%
  dplyr::select(Primary, Metastasis)

fpkmSamples <- rownames(datExpr0)
traitSamples <- rownames(traitData)
sameSample <- intersect(fpkmSamples, traitSamples)
datExpr0 <- datExpr0[sameSample, ]
datTraits <- traitData[sameSample, ]

enableWGCNAThreads()
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
softPower <- sft$powerEstimate

if (is.na(softPower)) {
  softPower <- 6
  cat("手动设置 softPower =", softPower, "\n")
}

net <- blockwiseModules(datExpr0, power = softPower,
                        TOMType = "unsigned", minModuleSize = 60,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "TOM",
                        verbose = 3)

moduleColors <- labels2colors(net$colors)

nGenes <- ncol(datExpr0)
nSamples <- nrow(datExpr0)
MEs0 <- moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# 📊 【论文图1】
pdf(file = file.path(wgcna_dir, "Fig1_Module_trait_heatmap.pdf"), width = 6, height = 8)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships"))
dev.off()

modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep="")
names(MMPvalue) <- paste("p.MM", modNames, sep="")

traitNames <- names(datTraits)
geneTraitSignificance <- as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", traitNames, sep="")
names(GSPvalue) <- paste("p.GS.", traitNames, sep="")

probes <- colnames(datExpr0)
geneInfo <- data.frame(
  GeneSymbol = probes,
  ModuleColor = moduleColors
)
geneInfo <- cbind(geneInfo, geneTraitSignificance, GSPvalue, geneModuleMembership, MMPvalue)
write.csv(geneInfo, file = file.path(wgcna_dir, "WGCNA_geneInfo.csv"), row.names = FALSE)

for (module in unique(moduleColors)) {
  modGenes <- probes[moduleColors == module]
  write.table(modGenes,
              file = file.path(wgcna_dir, paste0("module_", module, "_genes.txt")),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

cat("✅ WGCNA 分析完成！\n\n")

# ======================================================
# =========   第三部分: 联合分析与功能富集   =========
# ======================================================

cat("\n")
cat("========================================\n")
cat("  第三部分: 联合分析与功能富集\n")
cat("========================================\n\n")

library(VennDiagram)
library(grid)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

target_trait <- "Metastasis"

if (!target_trait %in% colnames(moduleTraitPvalue)) {
  stop("错误: 'Metastasis' 性状未找到")
}
target_module_index <- which.min(moduleTraitPvalue[, target_trait])
target_module_name_ME <- rownames(moduleTraitPvalue)[target_module_index]
target_module_color <- sub("ME", "", target_module_name_ME)

cat(paste0("🔍 最相关模块: '", target_module_color, "'\n"))

module_genes <- geneInfo %>%
  filter(ModuleColor == target_module_color) %>%
  pull(GeneSymbol)

degs_genes <- rownames(deg)
interGenes <- intersect(module_genes, degs_genes)

cat(paste0("📊 DEGs 基因数: ", length(degs_genes), "\n"))
cat(paste0("📊 模块基因数: ", length(module_genes), "\n"))
cat(paste0("🧬 交集基因数: ", length(interGenes), "\n\n"))

# 📊 【论文图2】
venn_obj <- venn.diagram(
  x = list(
    DEGs = degs_genes,
    ModuleGenes = module_genes
  ),
  filename = NULL,
  main = paste("Overlap between DEGs and", target_module_color, "Module"),
  main.cex = 1.5,
  category.names = c("DEGs", paste(target_module_color, "Module")),
  cat.cex = 1.2,
  cat.fontface = "bold",
  lwd = 2,
  lty = 'blank',
  fill = c("#E69F00", "#56B4E9"),
  cex = 1.5,
  fontface = "bold",
  margin = 0.1
)

venn_output_file <- file.path(wgcna_dir, "Fig2_Venn_DEG_vs_Module.pdf")
pdf(venn_output_file, width = 8, height = 8)
grid.draw(venn_obj)
dev.off()
cat("✅ Venn 图已保存:", venn_output_file, "\n\n")

write.table(interGenes, file=file.path(wgcna_dir, "interGenes_list.txt"),
            quote=F, row.names=F, col.names=F)

# 富集分析
interGenes_entrez <- bitr(interGenes, fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)

if (nrow(interGenes_entrez) > 0) {
  
  go_results <- enrichGO(gene = interGenes_entrez$ENTREZID,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID",
                         ont = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.01,
                         qvalueCutoff = 0.1)
  
  if (!is.null(go_results) && nrow(go_results) > 0) {
    go_plot <- dotplot(go_results, showCategory=15, split="ONTOLOGY") +
      facet_grid(ONTOLOGY~., scale="free") +
      ggtitle("GO Enrichment Analysis")
    ggsave(file.path(wgcna_dir, "Fig3_GO_enrichment_dotplot.pdf"), 
           plot = go_plot, width = 12, height = 15)
    
    write.csv(as.data.frame(go_results), 
              file.path(wgcna_dir, "GO_enrichment_results.csv"))
  }
  
  kegg_results <- enrichKEGG(gene = interGenes_entrez$ENTREZID,
                             organism = 'hsa',
                             pvalueCutoff = 0.01,
                             qvalueCutoff = 0.1)
  
  if (!is.null(kegg_results) && nrow(kegg_results) > 0) {
    kegg_plot <- barplot(kegg_results, showCategory=20, 
                         title="KEGG Enrichment Analysis")
    ggsave(file.path(wgcna_dir, "Fig4_KEGG_enrichment_barplot.pdf"), 
           plot = kegg_plot, width = 10, height = 8)
    
    write.csv(as.data.frame(kegg_results), 
              file.path(wgcna_dir, "KEGG_enrichment_results.csv"))
  }
}

cat("✅ 联合分析完成！\n\n")

# ======================================================
# =========    第四部分: 测试集预处理    =========
# ======================================================

cat("\n")
cat("========================================\n")
cat("  第四部分: 测试集预处理\n")
cat("========================================\n\n")

test_dir <- 'D:/dow/GSE49355'
setwd(test_dir)

gset_test <- getGEO('GSE49355', destdir=".", AnnotGPL = F, getGPL = F)

exp <- exprs(gset_test[[1]])
exp <- normalizeBetweenArrays(exp)
exp <- log2(exp+1)

gp1 <- read.table('GPL96-57554.txt', 
                  header = TRUE, fill = TRUE, sep = "\t",
                  quote = "", comment.char = "#", stringsAsFactors = FALSE)

ids <- gp1[,c('ID','Gene.Symbol')]
colnames(ids) <- c('probe_id','symbol')
ids <- ids[ids$symbol != "" & !is.na(ids$symbol), ]
ids$symbol <- str_split(ids$symbol, " /// ", simplify = TRUE)[, 1]

exp_filtered <- exp[rownames(exp) %in% ids$probe_id, ]
ids_filtered <- ids[ids$probe_id %in% rownames(exp_filtered), ]
ids_filtered <- ids_filtered[match(rownames(exp_filtered), ids_filtered$probe_id), ]

pdata <- pData(gset_test[[1]])
tumor_met <- pdata[pdata$`characteristics_ch1` %in% 
                     c("organism part: Primary Tumor", "organism part: Liver Metastasis"), ]

patient_counts <- table(tumor_met$`patient number:ch1`,
                        tumor_met$`characteristics_ch1`)
paired_patients <- rownames(patient_counts)[
  patient_counts[,"organism part: Liver Metastasis"] > 0 &
    patient_counts[,"organism part: Primary Tumor"] > 0
]

paired_samples_test <- tumor_met[tumor_met$`patient number:ch1` %in% paired_patients, ] %>%
  arrange(`patient number:ch1`, `characteristics_ch1`)
rownames(paired_samples_test) <- paired_samples_test$geo_accession

exp_anno_test <- as.data.frame(exp_filtered)
exp_anno_test$gene <- ids_filtered$symbol
exp_gene_test <- exp_anno_test %>%
  group_by(gene) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")
exp_gene_mat_test <- as.matrix(exp_gene_test %>% dplyr::select(-gene))
rownames(exp_gene_mat_test) <- exp_gene_test$gene

exp_test_paired <- exp_gene_mat_test[, rownames(paired_samples_test)]

group_test <- ifelse(grepl("Primary Tumor", paired_samples_test$`characteristics_ch1`), 
                     "Primary", "Metastasis")
group_test <- factor(group_test, levels = c("Primary", "Metastasis"))

cat("✅ 测试集预处理完成！\n\n")

setwd(main_dir)

# ======================================================
# =========     第五部分: LASSO 特征筛选     =========
# ======================================================

cat("\n")
cat("========================================\n")
cat("  第五部分: LASSO 特征筛选\n")
cat("========================================\n\n")

library(glmnet)
library(caret)
library(randomForest)
library(e1071)
library(xgboost)
library(pROC)

common_genes <- intersect(interGenes, rownames(exp_test_paired))
cat(paste0("🎯 共同基因数: ", length(common_genes), "\n"))

if (length(common_genes) < 10) {
  common_genes <- intersect(degs_genes, rownames(exp_test_paired))
  cat(paste0("使用差异基因: ", length(common_genes), " 个\n"))
}

X_train <- t(exp_paired[common_genes, ])
y_train <- ifelse(group == "Metastasis", 1, 0)
X_test <- t(exp_test_paired[common_genes, ])
y_test <- ifelse(group_test == "Metastasis", 1, 0)

X_train <- as.matrix(X_train)
X_test <- as.matrix(X_test)

set.seed(553)
cv_fit <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 1, nfolds = 5)

pdf(file.path(wgcna_dir, "Fig5_LASSO_cross_validation.pdf"), width = 8, height = 6)
plot(cv_fit)
dev.off()

optimal_lambda <- cv_fit$lambda.1se
lasso_coefs <- coef(cv_fit, s = optimal_lambda)
lasso_genes <- lasso_coefs@Dimnames[[1]][which(lasso_coefs != 0)]
lasso_genes <- lasso_genes[lasso_genes != "(Intercept)"]

cat(paste0("🎯 LASSO 筛选基因数: ", length(lasso_genes), "\n"))
write.table(lasso_genes, file=file.path(wgcna_dir, "LASSO_selected_genes.txt"),
            quote=F, row.names=F, col.names=F)

# 保存LASSO系数（关键：用于后续药物筛选）
lasso_coefs_full <- coef(cv_fit, s = optimal_lambda)
lasso_coefs_df <- data.frame(
  Gene = lasso_coefs_full@Dimnames[[1]][lasso_coefs_full@i + 1],
  Coefficient = lasso_coefs_full@x
)
lasso_coefs_df <- lasso_coefs_df[lasso_coefs_df$Gene != "(Intercept)", ]

# 分离上调和下调基因（关键：用于药物筛选）
upregulated_genes <- lasso_coefs_df$Gene[lasso_coefs_df$Coefficient > 0]
downregulated_genes <- lasso_coefs_df$Gene[lasso_coefs_df$Coefficient < 0]

cat(paste0("📈 上调基因数（促进转移）: ", length(upregulated_genes), "\n"))
cat(paste0("📉 下调基因数（抑制转移）: ", length(downregulated_genes), "\n"))

write.csv(lasso_coefs_df, 
          file = file.path(wgcna_dir, "LASSO_genes_with_coefficients.csv"),
          row.names = FALSE)

cat("\n✅ LASSO 特征筛选完成！\n\n")

# ======================================================
# =========     第六部分: 机器学习模型     =========
# ======================================================

cat("\n")
cat("========================================\n")
cat("  第六部分: 机器学习模型构建与验证\n")
cat("========================================\n\n")

X_train_lasso <- X_train[, lasso_genes, drop = FALSE]
X_test_lasso <- X_test[, lasso_genes, drop = FALSE]

train_data <- data.frame(X_train_lasso, 
                         Class = factor(y_train, levels = c(0, 1), 
                                        labels = c("Primary", "Metastasis")))
test_data <- data.frame(X_test_lasso,
                        Class = factor(y_test, levels = c(0, 1), 
                                       labels = c("Primary", "Metastasis")))

set.seed(123)
model_rf <- randomForest(Class ~ ., data = train_data, ntree = 500)

set.seed(123)
model_svm <- svm(Class ~ ., data = train_data, probability = TRUE, kernel = "radial")

dtrain <- xgb.DMatrix(data = X_train_lasso, label = y_train)
dtest <- xgb.DMatrix(data = X_test_lasso, label = y_test)
params <- list(objective = "binary:logistic", eval_metric = "auc", eta = 0.1)
set.seed(123)
model_xgb <- xgb.train(params, dtrain, nrounds = 100, 
                       watchlist = list(train = dtrain), verbose = 0)

pred_prob_rf <- predict(model_rf, newdata = test_data, type = "prob")[, "Metastasis"]
pred_prob_svm_attrs <- attr(predict(model_svm, newdata = test_data, 
                                    probability = TRUE), "probabilities")
pred_prob_svm <- pred_prob_svm_attrs[, "Metastasis"]
pred_prob_xgb <- predict(model_xgb, newdata = X_test_lasso)

roc_rf <- roc(test_data$Class, pred_prob_rf, levels = c("Primary", "Metastasis"))
roc_svm <- roc(test_data$Class, pred_prob_svm, levels = c("Primary", "Metastasis"))
roc_xgb <- roc(test_data$Class, pred_prob_xgb, levels = c("Primary", "Metastasis"))

cat(paste0("📊 RF AUC: ", round(auc(roc_rf), 3), "\n"))
cat(paste0("📊 SVM AUC: ", round(auc(roc_svm), 3), "\n"))
cat(paste0("📊 XGBoost AUC: ", round(auc(roc_xgb), 3), "\n\n"))

roc_plot <- ggroc(list(RF = roc_rf, SVM = roc_svm, XGBoost = roc_xgb), 
                  legacy.axes = TRUE) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color="grey", linetype="dashed") +
  labs(
    title = "ROC Curves on Independent Test Set",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  annotate("text", x = 0.6, y = 0.25, 
           label = paste("RF AUC =", round(auc(roc_rf), 3)), hjust = 0) +
  annotate("text", x = 0.6, y = 0.15, 
           label = paste("SVM AUC =", round(auc(roc_svm), 3)), hjust = 0) +
  annotate("text", x = 0.6, y = 0.05,
           label = paste("XGBoost AUC =", round(auc(roc_xgb), 3)), hjust = 0) +
  theme_bw() +
  coord_equal() +
  scale_color_discrete(name = "Model")

ggsave(file.path(wgcna_dir, "Fig6_Model_ROC_curves.pdf"), 
       plot = roc_plot, width = 7, height = 7)

pred_class_rf <- factor(ifelse(pred_prob_rf > 0.5, "Metastasis", "Primary"), 
                        levels = c("Primary", "Metastasis"))
cm_rf <- confusionMatrix(pred_class_rf, test_data$Class, positive = "Metastasis")
print(cm_rf)

cat("\n✅ 机器学习模型完成！\n\n")

# ======================================================
# =========     第七部分: 关键基因热图     =========
# ======================================================

cat("\n")
cat("========================================\n")
cat("  第七部分: 关键基因热图绘制\n")
cat("========================================\n\n")

library(pheatmap)

combined_data <- rbind(X_train_lasso, X_test_lasso)

annotation_col <- data.frame(
  Group = c(as.character(train_data$Class), as.character(test_data$Class)),
  Dataset = rep(c("Train", "Test"), 
                times = c(nrow(X_train_lasso), nrow(X_test_lasso)))
)
rownames(annotation_col) <- rownames(combined_data)

pheatmap(t(combined_data),
         scale = "row",
         annotation_col = annotation_col,
         show_colnames = FALSE,
         main = paste("Expression of", length(lasso_genes), "Signature Genes"),
         filename = file.path(wgcna_dir, "Fig7_Signature_Genes_Heatmap.pdf"),
         width = 10, height = 8
)

cat("✅ 关键基因热图完成！\n\n")

# ======================================================
# =========  第八部分: 上游调控网络分析（优化位置）  =========
# ======================================================

cat("\n")
cat("========================================\n")
cat("  第八部分: 上游调控网络分析\n")
cat("========================================\n\n")

library(igraph)
library(httr)
library(jsonlite)

# 创建调控网络输出目录
regulation_dir <- file.path(wgcna_dir, "Regulatory_Network")
if (!dir.exists(regulation_dir)) {
  dir.create(regulation_dir)
}

# ============================================
# 8.1 转录因子预测（基于TRRUST v2数据库）
# ============================================
cat("📊 Step 1: 转录因子预测分析...\n")

trrust_url <- "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"
trrust_file <- file.path(regulation_dir, "trrust_rawdata.human.tsv")

if (!file.exists(trrust_file)) {
  tryCatch({
    download.file(trrust_url, destfile = trrust_file, method = "auto", quiet = FALSE)
    cat("✅ TRRUST数据库下载成功\n")
  }, error = function(e) {
    cat("⚠️ 自动下载失败，请手动下载 TRRUST 数据库\n")
    cat("   URL:", trrust_url, "\n")
    cat("   保存到:", trrust_file, "\n")
  })
}

# 读取TRRUST数据
if (file.exists(trrust_file)) {
  trrust_data <- read.table(trrust_file, sep = "\t", header = FALSE, 
                            stringsAsFactors = FALSE, quote = "")
  colnames(trrust_data) <- c("TF", "Target", "Regulation", "PMID")
  
  # 筛选与关键基因相关的TF
  tf_target_pairs <- trrust_data %>%
    filter(Target %in% lasso_genes) %>%
    dplyr::select(TF, Target, Regulation)
  
  cat(paste0("🎯 发现 ", nrow(tf_target_pairs), " 个TF-基因调控关系\n"))
  cat(paste0("🎯 涉及 ", length(unique(tf_target_pairs$TF)), " 个转录因子\n"))
  
  # 统计每个TF调控的基因数量
  tf_stats <- tf_target_pairs %>%
    group_by(TF) %>%
    summarise(
      n_targets = n(),
      targets = paste(Target, collapse = ", "),
      regulation_types = paste(unique(Regulation), collapse = ", ")
    ) %>%
    arrange(desc(n_targets))
  
  write.csv(tf_target_pairs, 
            file.path(regulation_dir, "TF_target_pairs.csv"), 
            row.names = FALSE)
  write.csv(tf_stats, 
            file.path(regulation_dir, "TF_statistics.csv"), 
            row.names = FALSE)
  
  # 📊 【论文图8A】TF调控网络图
  if (nrow(tf_target_pairs) > 0) {
    tf_edges <- tf_target_pairs %>%
      dplyr::select(from = TF, to = Target, type = Regulation)
    
    tf_graph <- graph_from_data_frame(tf_edges, directed = TRUE)
    
    # 节点属性
    V(tf_graph)$node_type <- ifelse(V(tf_graph)$name %in% lasso_genes, 
                                    "Target Gene", "Transcription Factor")
    V(tf_graph)$color <- ifelse(V(tf_graph)$node_type == "Target Gene", 
                                "#E69F00", "#56B4E9")
    V(tf_graph)$size <- ifelse(V(tf_graph)$node_type == "Target Gene", 12, 8)
    
    # 边属性
    E(tf_graph)$color <- ifelse(E(tf_graph)$type == "Activation", 
                                "darkgreen", "darkred")
    E(tf_graph)$arrow.size <- 0.6
    
    # 使用PNG格式避免字体问题
    png(file.path(regulation_dir, "Fig8A_TF_Regulatory_Network.png"), 
        width = 1400, height = 1200, res = 120)
    
    par(mar = c(1, 1, 3, 1))
    set.seed(123)
    
    plot(tf_graph,
         vertex.label.cex = 0.75,
         vertex.label.color = "black",
         vertex.label.family = "sans",
         vertex.frame.color = "white",
         edge.curved = 0.2,
         layout = layout_with_fr(tf_graph),
         main = "Transcription Factor Regulatory Network")
    
    legend("bottomleft",
           legend = c("Target Gene", "Transcription Factor", 
                      "Activation", "Repression"),
           col = c("#E69F00", "#56B4E9", "darkgreen", "darkred"),
           pch = c(19, 19, NA, NA),
           lty = c(NA, NA, 1, 1),
           pt.cex = 2,
           cex = 0.9,
           bty = "n")
    
    dev.off()
    cat("✅ Fig8A: TF调控网络图已生成 (PNG格式)\n")
  }
}

# ============================================
# 8.2 miRNA预测（使用multiMiR包）
# ============================================
cat("\n📊 Step 2: microRNA预测分析...\n")

if (!require("multiMiR", quietly = TRUE)) {
  tryCatch({
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("multiMiR", ask = FALSE)
  }, error = function(e) {
    cat("⚠️ multiMiR安装失败，请手动安装\n")
  })
}

if (require("multiMiR", quietly = TRUE)) {
  cat("🔍 查询miRNA靶基因数据库...\n")
  
  tryCatch({
    mirna_results <- get_multimir(
      org = "hsa",
      target = lasso_genes,
      table = "validated",
      summary = FALSE
    )
    
    if (!is.null(mirna_results) && nrow(mirna_results@data) > 0) {
      mirna_data <- mirna_results@data
      
      cat(paste0("🎯 发现 ", nrow(mirna_data), " 个miRNA-基因相互作用\n"))
      cat(paste0("🎯 涉及 ", length(unique(mirna_data$mature_mirna_id)), " 个miRNA\n"))
      
      # 数据清理
      mirna_target_pairs <- mirna_data %>%
        dplyr::select(miRNA = mature_mirna_id, 
                      Target = target_symbol,
                      Database = database) %>%
        distinct()
      
      # 统计
      mirna_stats <- mirna_target_pairs %>%
        group_by(miRNA) %>%
        summarise(
          n_targets = n(),
          targets = paste(unique(Target), collapse = ", "),
          databases = paste(unique(Database), collapse = ", ")
        ) %>%
        arrange(desc(n_targets))
      
      write.csv(mirna_target_pairs, 
                file.path(regulation_dir, "miRNA_target_pairs.csv"), 
                row.names = FALSE)
      write.csv(mirna_stats, 
                file.path(regulation_dir, "miRNA_statistics.csv"), 
                row.names = FALSE)
      
      # 📊 【论文图8B】miRNA调控网络图
      # 只保留调控≥2个基因的miRNA
      key_mirnas <- mirna_stats %>%
        filter(n_targets >= 2) %>%
        pull(miRNA)
      
      if (length(key_mirnas) == 0) {
        key_mirnas <- head(mirna_stats$miRNA, 15)
      }
      
      mirna_edges_filtered <- mirna_target_pairs %>%
        filter(miRNA %in% key_mirnas) %>%
        dplyr::select(from = miRNA, to = Target)
      
      mirna_graph <- graph_from_data_frame(mirna_edges_filtered, directed = TRUE)
      
      # 节点属性
      V(mirna_graph)$node_type <- ifelse(V(mirna_graph)$name %in% lasso_genes, 
                                         "Target Gene", "miRNA")
      V(mirna_graph)$color <- ifelse(V(mirna_graph)$node_type == "Target Gene", 
                                     "#E69F00", "#009E73")
      V(mirna_graph)$size <- ifelse(V(mirna_graph)$node_type == "Target Gene", 12, 8)
      
      # 使用PNG避免字体问题
      png(file.path(regulation_dir, "Fig8B_miRNA_Regulatory_Network.png"), 
          width = 1400, height = 1200, res = 120)
      
      par(mar = c(1, 1, 3, 1))
      set.seed(456)
      
      plot(mirna_graph,
           vertex.label.cex = 0.7,
           vertex.label.color = "black",
           vertex.label.family = "sans",
           vertex.frame.color = "white",
           edge.color = "grey60",
           edge.arrow.size = 0.6,
           edge.curved = 0.2,
           layout = layout_with_fr(mirna_graph),
           main = paste("miRNA Regulatory Network\n(",
                        length(key_mirnas), "key miRNAs)"))
      
      legend("bottomleft",
             legend = c("Target Gene", "miRNA"),
             col = c("#E69F00", "#009E73"),
             pch = 19,
             pt.cex = 2,
             cex = 0.9,
             bty = "n")
      
      dev.off()
      cat("✅ Fig8B: miRNA调控网络图已生成 (PNG格式)\n")
      
    } else {
      cat("⚠️ multiMiR查询未返回数据\n")
    }
    
  }, error = function(e) {
    cat("⚠️ miRNA分析出错:", e$message, "\n")
  })
}

# ============================================
# 8.3 整合TF和miRNA调控网络
# ============================================
cat("\n📊 Step 3: 整合调控网络...\n")

if (exists("tf_target_pairs") && exists("mirna_target_pairs")) {
  
  # 合并数据
  tf_edges <- tf_target_pairs %>%
    mutate(Regulator_Type = "TF") %>%
    dplyr::select(Regulator = TF, Target, Type = Regulation, Regulator_Type)
  
  mirna_edges <- mirna_target_pairs %>%
    mutate(Regulation = "Repression", Regulator_Type = "miRNA") %>%
    dplyr::select(Regulator = miRNA, Target, Type = Regulation, Regulator_Type) %>%
    distinct()
  
  # 筛选关键调控因子（调控≥2个基因）
  regulator_counts <- rbind(
    tf_edges %>% count(Regulator, Regulator_Type),
    mirna_edges %>% count(Regulator, Regulator_Type)
  ) %>%
    arrange(desc(n))
  
  key_regulators <- regulator_counts %>%
    filter(n >= 2) %>%
    pull(Regulator)
  
  if (length(key_regulators) < 10) {
    key_regulators <- head(regulator_counts$Regulator, 20)
  }
  
  integrated_edges <- rbind(tf_edges, mirna_edges) %>%
    filter(Regulator %in% key_regulators)
  
  if (nrow(integrated_edges) > 0) {
    # 创建网络
    integrated_graph <- graph_from_data_frame(
      integrated_edges %>% dplyr::select(from = Regulator, to = Target), 
      directed = TRUE
    )
    
    # 节点属性
    node_attrs <- data.frame(
      name = V(integrated_graph)$name,
      stringsAsFactors = FALSE
    ) %>%
      left_join(
        integrated_edges %>% 
          dplyr::select(name = Regulator, Regulator_Type) %>% 
          distinct(),
        by = "name"
      ) %>%
      mutate(
        node_type = case_when(
          name %in% lasso_genes ~ "Target Gene",
          Regulator_Type == "TF" ~ "Transcription Factor",
          Regulator_Type == "miRNA" ~ "miRNA",
          TRUE ~ "Unknown"
        )
      )
    
    V(integrated_graph)$node_type <- node_attrs$node_type
    V(integrated_graph)$color <- case_when(
      V(integrated_graph)$node_type == "Target Gene" ~ "#E69F00",
      V(integrated_graph)$node_type == "Transcription Factor" ~ "#56B4E9",
      V(integrated_graph)$node_type == "miRNA" ~ "#009E73",
      TRUE ~ "grey"
    )
    V(integrated_graph)$size <- ifelse(V(integrated_graph)$node_type == "Target Gene", 12, 8)
    
    # 📊 【论文图8C】整合调控网络（PNG格式）
    png(file.path(regulation_dir, "Fig8C_Integrated_Regulatory_Network.png"), 
        width = 1800, height = 1600, res = 120)
    
    par(mar = c(1, 1, 3, 1))
    set.seed(999)
    
    plot(integrated_graph,
         vertex.label.cex = 0.65,
         vertex.label.color = "black",
         vertex.label.family = "sans",
         vertex.frame.color = "white",
         edge.color = "grey60",
         edge.arrow.size = 0.5,
         edge.curved = 0.2,
         layout = layout_with_fr(integrated_graph),
         main = "Integrated TF-miRNA Regulatory Network")
    
    legend("bottomleft",
           legend = c("Target Gene", "Transcription Factor", "miRNA"),
           col = c("#E69F00", "#56B4E9", "#009E73"),
           pch = 19,
           pt.cex = 2,
           cex = 0.9,
           bty = "n")
    
    dev.off()
    cat("✅ Fig8C: 整合调控网络图已生成 (PNG格式)\n")
    
    # 保存数据
    write.csv(integrated_edges, 
              file.path(regulation_dir, "Integrated_regulatory_pairs.csv"), 
              row.names = FALSE)
    
    # 导出Cytoscape格式
    write.table(integrated_edges %>% 
                  dplyr::select(source = Regulator, target = Target, 
                                interaction = Type, regulator_type = Regulator_Type),
                file.path(regulation_dir, "Cytoscape_edges.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    node_table <- data.frame(
      id = V(integrated_graph)$name,
      node_type = V(integrated_graph)$node_type
    )
    write.table(node_table,
                file.path(regulation_dir, "Cytoscape_nodes.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    cat("✅ Cytoscape网络文件已导出\n")
  }
  
} else if (exists("tf_target_pairs")) {
  cat("⚠️ 仅TF数据可用，已生成TF网络\n")
} else if (exists("mirna_target_pairs")) {
  cat("⚠️ 仅miRNA数据可用，已生成miRNA网络\n")
}

# ============================================
# 8.4 识别关键调控因子（供药物筛选使用）
# ============================================
cat("\n📊 Step 4: 识别关键调控因子...\n")

# 提取关键TF
if (exists("tf_stats")) {
  key_TFs <- tf_stats %>% 
    filter(n_targets >= 2) %>% 
    pull(TF)
  
  if (length(key_TFs) == 0) {
    key_TFs <- head(tf_stats$TF, 10)
  }
  
  cat(paste0("🎯 识别到 ", length(key_TFs), " 个关键转录因子\n"))
  cat("   Top 10 关键TFs:\n")
  print(head(tf_stats, 10))
}

# 提取关键miRNA
if (exists("mirna_stats")) {
  key_miRNAs <- mirna_stats %>% 
    filter(n_targets >= 2) %>% 
    pull(miRNA)
  
  if (length(key_miRNAs) == 0) {
    key_miRNAs <- head(mirna_stats$miRNA, 15)
  }
  
  cat(paste0("\n🎯 识别到 ", length(key_miRNAs), " 个关键miRNA\n"))
  cat("   Top 15 关键miRNAs:\n")
  print(head(mirna_stats, 15))
}

cat("\n✅ 上游调控网络分析完成！\n\n")

# ======================================================
# =========  第九部分: 多层次药物筛选（优化位置）  =========
# ======================================================

cat("\n")
cat("========================================\n")
cat("  第九部分: 多层次药物筛选策略\n")
cat("========================================\n\n")

library(enrichR)
library(ggrepel)

# 工具函数：提取药物名称
extract_drug_name <- function(term) {
  drug <- str_split(term, "-")[[1]][1]
  drug <- tolower(trimws(drug))
  return(drug)
}

# ============================================
# 9.1 Layer 1: 直接靶向核心基因的药物
# ============================================
cat("========================================\n")
cat("  Layer 1: 直接靶向核心基因的药物\n")
cat("========================================\n\n")

drug_db <- "DSigDB"

cat(paste0("📈 上调基因数（促进转移）: ", length(upregulated_genes), "\n"))
cat(paste0("📉 下调基因数（抑制转移）: ", length(downregulated_genes), "\n"))
cat(paste0("   上调基因: ", paste(upregulated_genes, collapse = ", "), "\n"))
cat(paste0("   下调基因: ", paste(downregulated_genes, collapse = ", "), "\n\n"))

# 保存基因列表
write.table(upregulated_genes, 
            file = file.path(wgcna_dir, "L1_Upregulated_genes.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(downregulated_genes, 
            file = file.path(wgcna_dir, "L1_Downregulated_genes.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# 分析上调基因（寻找抑制性药物）
drugs_for_up <- NULL
if (length(upregulated_genes) > 0) {
  cat("🔍 分析上调基因（寻找抑制性药物）...\n")
  enriched_up <- enrichr(upregulated_genes, drug_db)
  
  if (!is.null(enriched_up[[drug_db]]) && nrow(enriched_up[[drug_db]]) > 0) {
    write.csv(enriched_up[[drug_db]], 
              file = file.path(wgcna_dir, "L1_EnrichR_Upregulated_DSigDB.csv"),
              row.names = FALSE)
    
    drugs_for_up <- enriched_up[[drug_db]] %>%
      filter(Adjusted.P.value < 0.05) %>%
      mutate(
        Drug = sapply(Term, extract_drug_name),
        Target_Type = "Inhibit_Upregulated",
        Target_Genes = Genes,
        Target_Level = "Direct_Gene"
      ) %>%
      dplyr::select(Drug, Term, Target_Type, Target_Level, Target_Genes, 
                    Adjusted.P.value, Combined.Score)
    
    cat(paste0("   ✓ 找到 ", nrow(drugs_for_up), " 个显著条目\n"))
    cat(paste0("   ✓ 涉及 ", length(unique(drugs_for_up$Drug)), " 个不同药物\n"))
  }
}

# 分析下调基因（寻找激活性药物）
drugs_for_down <- NULL
if (length(downregulated_genes) > 0) {
  cat("\n🔍 分析下调基因（寻找激活性药物）...\n")
  enriched_down <- enrichr(downregulated_genes, drug_db)
  
  if (!is.null(enriched_down[[drug_db]]) && nrow(enriched_down[[drug_db]]) > 0) {
    write.csv(enriched_down[[drug_db]], 
              file = file.path(wgcna_dir, "L1_EnrichR_Downregulated_DSigDB.csv"),
              row.names = FALSE)
    
    drugs_for_down <- enriched_down[[drug_db]] %>%
      filter(Adjusted.P.value < 0.05) %>%
      mutate(
        Drug = sapply(Term, extract_drug_name),
        Target_Type = "Activate_Downregulated",
        Target_Genes = Genes,
        Target_Level = "Direct_Gene"
      ) %>%
      dplyr::select(Drug, Term, Target_Type, Target_Level, Target_Genes, 
                    Adjusted.P.value, Combined.Score)
    
    cat(paste0("   ✓ 找到 ", nrow(drugs_for_down), " 个显著条目\n"))
    cat(paste0("   ✓ 涉及 ", length(unique(drugs_for_down$Drug)), " 个不同药物\n"))
  }
}

# 寻找双向调控药物
dual_action_drugs <- NULL
if (!is.null(drugs_for_up) && !is.null(drugs_for_down)) {
  common_drugs <- intersect(unique(drugs_for_up$Drug), 
                            unique(drugs_for_down$Drug))
  
  if (length(common_drugs) > 0) {
    cat(paste0("\n🎯 找到 ", length(common_drugs), " 个双向调控药物\n"))
    
    dual_action_drugs <- data.frame()
    for (drug in common_drugs) {
      info_up <- drugs_for_up %>% 
        filter(Drug == drug) %>% 
        arrange(Adjusted.P.value) %>%
        head(1)
      
      info_down <- drugs_for_down %>% 
        filter(Drug == drug) %>% 
        arrange(Adjusted.P.value) %>%
        head(1)
      
      combined_info <- data.frame(
        Drug = drug,
        Target_Level = "Direct_Gene",
        Inhibit_Term = info_up$Term,
        Inhibit_Targets = info_up$Target_Genes,
        Inhibit_Pvalue = info_up$Adjusted.P.value,
        Inhibit_Score = info_up$Combined.Score,
        Activate_Term = info_down$Term,
        Activate_Targets = info_down$Target_Genes,
        Activate_Pvalue = info_down$Adjusted.P.value,
        Activate_Score = info_down$Combined.Score,
        Mean_Score = mean(c(info_up$Combined.Score, info_down$Combined.Score)),
        Mean_Pvalue = mean(c(info_up$Adjusted.P.value, info_down$Adjusted.P.value))
      )
      
      dual_action_drugs <- rbind(dual_action_drugs, combined_info)
    }
    
    dual_action_drugs <- dual_action_drugs %>%
      arrange(desc(Mean_Score), Mean_Pvalue)
    
    write.csv(dual_action_drugs,
              file = file.path(wgcna_dir, "L1_Dual_Action_Drugs_DSigDB.csv"),
              row.names = FALSE)
    
    cat("Top 10 双向调控候选药物:\n")
    print(dual_action_drugs[1:min(10, nrow(dual_action_drugs)), 
                            c("Drug", "Mean_Score", "Mean_Pvalue")])
  }
}

cat("\n✅ Layer 1 完成！\n\n")

# ============================================
# 9.2 Layer 2: 靶向关键转录因子的药物
# ============================================
cat("========================================\n")
cat("  Layer 2: 靶向关键转录因子的药物\n")
cat("========================================\n\n")

tf_drugs <- NULL

if (exists("key_TFs") && length(key_TFs) > 0) {
  
  cat(paste0("🔍 分析 ", length(key_TFs), " 个关键TF...\n"))
  cat(paste0("   关键TFs: ", paste(head(key_TFs, 15), collapse = ", "), "\n\n"))
  
  # 使用多个数据库
  tf_drug_dbs <- c("DSigDB", "ChEA_2016")
  
  tf_drug_results <- enrichr(key_TFs, tf_drug_dbs)
  
  # 处理DSigDB结果
  if (!is.null(tf_drug_results[["DSigDB"]]) && 
      nrow(tf_drug_results[["DSigDB"]]) > 0) {
    
    tf_drugs <- tf_drug_results[["DSigDB"]] %>%
      filter(Adjusted.P.value < 0.05) %>%
      mutate(
        Drug = sapply(Term, extract_drug_name),
        Target_Level = "Transcription_Factor",
        TF_Targets = Genes,
        Mechanism = "Targeting TF"
      ) %>%
      arrange(Adjusted.P.value) %>%
      head(50)
    
    write.csv(tf_drugs, 
              file.path(wgcna_dir, "L2_Drugs_targeting_TFs.csv"),
              row.names = FALSE)
    
    cat(paste0("   ✓ 找到 ", nrow(tf_drugs), " 个靶向TF的药物候选\n"))
    cat(paste0("   ✓ 涉及 ", length(unique(tf_drugs$Drug)), " 个不同药物\n"))
  }
  
  # 处理ChEA结果
  if (!is.null(tf_drug_results[["ChEA_2016"]]) && 
      nrow(tf_drug_results[["ChEA_2016"]]) > 0) {
    
    chea_results <- tf_drug_results[["ChEA_2016"]] %>%
      filter(Adjusted.P.value < 0.05)
    
    write.csv(chea_results,
              file.path(wgcna_dir, "L2_ChEA_TF_enrichment.csv"),
              row.names = FALSE)
    
    cat(paste0("   ✓ ChEA富集: ", nrow(chea_results), " 个条目\n"))
  }
  
} else {
  cat("⚠️ 未找到关键TF，跳过Layer 2\n")
}

cat("\n✅ Layer 2 完成！\n\n")

# ============================================
# 9.3 Layer 3: miRNA治疗策略
# ============================================
cat("========================================\n")
cat("  Layer 3: miRNA治疗策略\n")
cat("========================================\n\n")

mirna_therapeutic_strategy <- NULL

if (exists("key_miRNAs") && length(key_miRNAs) > 0 && 
    exists("mirna_target_pairs")) {
  
  cat(paste0("🔍 分析 ", length(key_miRNAs), " 个关键miRNA...\n\n"))
  
  mirna_therapeutic_strategy <- data.frame()
  
  for (mirna in key_miRNAs) {
    # 获取该miRNA调控的核心基因
    targets <- mirna_target_pairs %>%
      filter(miRNA == mirna) %>%
      pull(Target)
    
    targets <- targets[targets %in% lasso_genes]
    
    if (length(targets) == 0) next
    
    # 判断这些靶基因在LASSO中的方向
    target_directions <- lasso_coefs_df %>%
      filter(Gene %in% targets) %>%
      summarise(
        n_upregulated = sum(Coefficient > 0),
        n_downregulated = sum(Coefficient < 0),
        mean_coef = mean(Coefficient)
      )
    
    # 制定策略
    # miRNA通常抑制靶基因表达
    # 如果靶基因上调(促转移) -> 需要用miRNA mimic抑制它们
    # 如果靶基因下调(抑转移) -> 需要用antagomiR抑制miRNA，让基因上调
    strategy <- if (target_directions$mean_coef > 0) {
      "Use_miRNA_Mimic"  # 靶基因上调 -> 用mimic抑制基因
    } else {
      "Use_AntagomiR"    # 靶基因下调 -> 用antagomiR阻断miRNA
    }
    
    mirna_therapeutic_strategy <- rbind(
      mirna_therapeutic_strategy,
      data.frame(
        miRNA = mirna,
        n_targets = length(targets),
        Strategy = strategy,
        Targets = paste(targets, collapse = ", "),
        Mean_Target_Coef = target_directions$mean_coef,
        n_upregulated = target_directions$n_upregulated,
        n_downregulated = target_directions$n_downregulated
      )
    )
  }
  
  if (nrow(mirna_therapeutic_strategy) > 0) {
    write.csv(mirna_therapeutic_strategy,
              file.path(wgcna_dir, "L3_miRNA_therapeutic_strategy.csv"),
              row.names = FALSE)
    
    cat(paste0("   ✓ 制定 ", nrow(mirna_therapeutic_strategy), 
               " 个miRNA治疗策略\n"))
    
    # 统计
    n_mimic <- sum(mirna_therapeutic_strategy$Strategy == "Use_miRNA_Mimic")
    n_antagomir <- sum(mirna_therapeutic_strategy$Strategy == "Use_AntagomiR")
    
    cat(paste0("      - 需使用miRNA mimic: ", n_mimic, " 个\n"))
    cat(paste0("      - 需使用antagomiR: ", n_antagomir, " 个\n\n"))
    
    cat("Top 10 miRNA治疗策略:\n")
    print(head(mirna_therapeutic_strategy %>% 
                 arrange(desc(n_targets)), 10))
    
    # 📊 可视化
    p_mirna_strategy <- mirna_therapeutic_strategy %>%
      head(20) %>%
      ggplot(aes(x = reorder(miRNA, n_targets),
                 y = n_targets,
                 fill = Strategy)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(
        values = c("Use_miRNA_Mimic" = "#2ecc71", 
                   "Use_AntagomiR" = "#e74c3c"),
        labels = c("Use miRNA mimic\n(suppress upregulated genes)", 
                   "Use antagomiR\n(release downregulated genes)")
      ) +
      coord_flip() +
      labs(
        title = "miRNA Therapeutic Strategy",
        subtitle = "Based on target gene expression patterns",
        x = "miRNA",
        y = "Number of Target Genes",
        fill = "Recommended Strategy"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "bottom",
        legend.text = element_text(size = 9)
      )
    
    ggsave(file.path(wgcna_dir, "Fig9_miRNA_Therapeutic_Strategy.pdf"),
           plot = p_mirna_strategy, width = 10, height = 10)
    
    cat("\n✅ Fig9: miRNA治疗策略图已生成\n")
  }
  
} else {
  cat("⚠️ 未找到关键miRNA，跳过Layer 3\n")
}

cat("\n✅ Layer 3 完成！\n\n")

# ============================================
# 9.4 综合药物推荐（整合三层结果）
# ============================================
cat("========================================\n")
cat("  综合药物推荐（三层整合）\n")
cat("========================================\n\n")

comprehensive_drugs <- data.frame()

# Layer 1: 直接靶向核心基因的药物
if (exists("dual_action_drugs") && !is.null(dual_action_drugs) && nrow(dual_action_drugs) > 0) {
  layer1 <- dual_action_drugs %>%
    head(15) %>%
    mutate(
      Drug_Name = Drug,
      Target_Level = "Direct_Gene",
      Mechanism = "Dual-action on genes",
      Priority_Score = Mean_Score * 1.5,  # 最高优先级
      Evidence_Strength = "Strong",
      Layer = "Layer 1"
    ) %>%
    dplyr::select(Drug_Name, Target_Level, Mechanism, Priority_Score, Evidence_Strength, Layer)
  
  comprehensive_drugs <- rbind(comprehensive_drugs, layer1)
  cat(paste0("✓ Layer 1贡献: ", nrow(layer1), " 个药物\n"))
}

# Layer 2: 靶向TF的药物
if (exists("tf_drugs") && !is.null(tf_drugs) && nrow(tf_drugs) > 0) {
  layer2 <- tf_drugs %>%
    head(15) %>%
    mutate(
      Drug_Name = Drug,
      Target_Level = "Transcription_Factor",
      Mechanism = paste0("Targeting TF: ", substr(TF_Targets, 1, 40), "..."),
      Priority_Score = Combined.Score * 1.2,  # 中等优先级
      Evidence_Strength = "Moderate",
      Layer = "Layer 2"
    ) %>%
    dplyr::select(Drug_Name, Target_Level, Mechanism, Priority_Score, Evidence_Strength, Layer)
  
  comprehensive_drugs <- rbind(comprehensive_drugs, layer2)
  cat(paste0("✓ Layer 2贡献: ", nrow(layer2), " 个药物\n"))
}

# Layer 3: miRNA治疗
if (exists("mirna_therapeutic_strategy") && !is.null(mirna_therapeutic_strategy) && 
    nrow(mirna_therapeutic_strategy) > 0) {
  layer3 <- mirna_therapeutic_strategy %>%
    head(15) %>%
    mutate(
      Drug_Name = paste0(miRNA, " ", 
                         ifelse(Strategy == "Use_miRNA_Mimic", 
                                "mimic", "antagomiR")),
      Target_Level = "microRNA",
      Mechanism = paste0(Strategy, " (targets: ", n_targets, " genes)"),
      Priority_Score = n_targets * abs(Mean_Target_Coef) * 15,
      Evidence_Strength = "Experimental",
      Layer = "Layer 3"
    ) %>%
    dplyr::select(Drug_Name, Target_Level, Mechanism, Priority_Score, Evidence_Strength, Layer)
  
  comprehensive_drugs <- rbind(comprehensive_drugs, layer3)
  cat(paste0("✓ Layer 3贡献: ", nrow(layer3), " 个策略\n\n"))
}

if (nrow(comprehensive_drugs) > 0) {
  # 去重并排序
  comprehensive_drugs <- comprehensive_drugs %>%
    distinct(Drug_Name, .keep_all = TRUE) %>%
    arrange(desc(Priority_Score))
  
  write.csv(comprehensive_drugs,
            file.path(wgcna_dir, "Comprehensive_Drug_Recommendations.csv"),
            row.names = FALSE)
  
  cat("🎯 最终推荐结果:\n")
  cat("========================================\n")
  cat(paste0("   总候选药物/策略数: ", nrow(comprehensive_drugs), "\n"))
  cat(paste0("      - 直接作用(Layer 1): ", 
             sum(comprehensive_drugs$Target_Level == "Direct_Gene"), "\n"))
  cat(paste0("      - TF层面(Layer 2): ", 
             sum(comprehensive_drugs$Target_Level == "Transcription_Factor"), "\n"))
  cat(paste0("      - miRNA层面(Layer 3): ", 
             sum(comprehensive_drugs$Target_Level == "microRNA"), "\n\n"))
  
  cat("Top 15 综合推荐:\n")
  cat("----------------------------------------\n")
  for (i in 1:min(15, nrow(comprehensive_drugs))) {
    cat(paste0(i, ". ", comprehensive_drugs$Drug_Name[i], "\n"))
    cat(paste0("   层级: ", comprehensive_drugs$Target_Level[i], 
               " (", comprehensive_drugs$Layer[i], ")\n"))
    cat(paste0("   机制: ", substr(comprehensive_drugs$Mechanism[i], 1, 60), "\n"))
    cat(paste0("   评分: ", round(comprehensive_drugs$Priority_Score[i], 2), "\n"))
    cat(paste0("   证据强度: ", comprehensive_drugs$Evidence_Strength[i], "\n\n"))
  }
  
  # 📊 【论文图10】综合推荐图
  p_comprehensive <- comprehensive_drugs %>%
    head(25) %>%
    ggplot(aes(x = reorder(Drug_Name, Priority_Score),
               y = Priority_Score,
               fill = Target_Level)) +
    geom_bar(stat = "identity", alpha = 0.85) +
    scale_fill_manual(
      values = c("Direct_Gene" = "#3498db",
                 "Transcription_Factor" = "#e74c3c",
                 "microRNA" = "#2ecc71"),
      name = "Target Level",
      labels = c("Direct Gene (Layer 1)",
                 "Transcription Factor (Layer 2)",
                 "microRNA (Layer 3)")
    ) +
    coord_flip() +
    labs(
      title = "Comprehensive Multi-Level Drug Recommendations",
      subtitle = "Integrating gene, TF, and miRNA targeting strategies",
      x = "Drug / Therapeutic Strategy",
      y = "Priority Score"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "bottom",
      axis.text.y = element_text(size = 8)
    )
  
  ggsave(file.path(wgcna_dir, "Fig10_Comprehensive_Drug_Recommendations.pdf"),
         plot = p_comprehensive, width = 12, height = 12)
  
  cat("✅ Fig10: 综合药物推荐图已生成\n\n")
  
  # 📊 【论文图11】三层药物对比气泡图
  if (nrow(comprehensive_drugs) >= 15) {
    p_bubble <- comprehensive_drugs %>%
      head(30) %>%
      mutate(rank = row_number()) %>%
      ggplot(aes(x = rank, y = Priority_Score, 
                 size = Priority_Score, color = Target_Level)) +
      geom_point(alpha = 0.7) +
      ggrepel::geom_text_repel(aes(label = Drug_Name), size = 3, 
                               max.overlaps = 25,
                               box.padding = 0.5) +
      scale_size_continuous(range = c(4, 15), name = "Priority Score") +
      scale_color_manual(
        values = c("Direct_Gene" = "#3498db",
                   "Transcription_Factor" = "#e74c3c",
                   "microRNA" = "#2ecc71"),
        name = "Target Level"
      ) +
      labs(
        title = "Multi-Level Drug Landscape",
        subtitle = "Three-layer therapeutic targeting strategy",
        x = "Rank",
        y = "Priority Score"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "right"
      )
    
    ggsave(file.path(wgcna_dir, "Fig11_Drug_Landscape_Bubble.pdf"),
           plot = p_bubble, width = 14, height = 10)
    
    cat("✅ Fig11: 药物景观气泡图已生成\n\n")
  }
  
} else {
  cat("⚠️ 未生成任何药物推荐\n")
}

# ============================================
# 9.5 生成详细药物筛选报告
# ============================================
cat("\n📝 生成药物筛选分析报告...\n")

report_lines <- c(
  "========================================",
  "   多层次药物筛选分析报告",
  "========================================",
  paste("\n分析日期:", Sys.Date()),
  paste("核心基因数:", length(lasso_genes)),
  paste("核心基因列表:", paste(lasso_genes, collapse = ", ")),
  "\n----------------------------------------",
  "\n[1] Layer 1: 直接靶向核心基因"
)

if (exists("dual_action_drugs") && nrow(dual_action_drugs) > 0) {
  report_lines <- c(report_lines,
                    paste("   - 双向调控药物:", nrow(dual_action_drugs), "个"),
                    "   - Top 5 推荐:")
  for (i in 1:min(5, nrow(dual_action_drugs))) {
    report_lines <- c(report_lines,
                      paste0("     ", i, ". ", dual_action_drugs$Drug[i],
                             " (Score: ", round(dual_action_drugs$Mean_Score[i], 2), ")"))
  }
}

if (exists("tf_drugs") && !is.null(tf_drugs) && nrow(tf_drugs) > 0) {
  report_lines <- c(report_lines,
                    "\n[2] Layer 2: 靶向关键转录因子",
                    paste("   - 关键TF数:", length(key_TFs), "个"),
                    paste("   - 靶向TF的药物:", nrow(tf_drugs), "个"),
                    "   - Top 5 推荐:")
  for (i in 1:min(5, nrow(tf_drugs))) {
    report_lines <- c(report_lines,
                      paste0("     ", i, ". ", tf_drugs$Drug[i],
                             " -> ", tf_drugs$TF_Targets[i]))
  }
}

if (exists("mirna_therapeutic_strategy") && !is.null(mirna_therapeutic_strategy) && 
    nrow(mirna_therapeutic_strategy) > 0) {
  report_lines <- c(report_lines,
                    "\n[3] Layer 3: miRNA治疗策略",
                    paste("   - 关键miRNA数:", nrow(mirna_therapeutic_strategy), "个"),
                    paste("   - miRNA mimic策略:", 
                          sum(mirna_therapeutic_strategy$Strategy == "Use_miRNA_Mimic"), "个"),
                    paste("   - AntagomiR策略:", 
                          sum(mirna_therapeutic_strategy$Strategy == "Use_AntagomiR"), "个"),
                    "   - Top 5 推荐:")
  for (i in 1:min(5, nrow(mirna_therapeutic_strategy))) {
    report_lines <- c(report_lines,
                      paste0("     ", i, ". ", mirna_therapeutic_strategy$miRNA[i],
                             " (", mirna_therapeutic_strategy$Strategy[i], 
                             ", targets: ", mirna_therapeutic_strategy$n_targets[i], ")"))
  }
}

report_lines <- c(report_lines,
                  "\n----------------------------------------",
                  "[4] 综合推荐",
                  paste("   总候选数:", nrow(comprehensive_drugs), "个"),
                  "\n   临床转化建议:",
                  "   • 短期策略: 使用Layer 1的成熟小分子药物",
                  "   • 中期策略: 开发Layer 2的TF靶向药物",
                  "   • 长期策略: 研发Layer 3的miRNA疗法",
                  "\n   联合用药建议:",
                  "   • Layer 1 + Layer 3: 小分子药物联合miRNA治疗",
                  "   • Layer 2 + Layer 3: TF调控剂联合miRNA治疗",
                  "\n----------------------------------------",
                  "[5] 输出文件列表",
                  "\n   Layer 1 (直接靶向基因):",
                  "   - L1_Upregulated_genes.txt",
                  "   - L1_Downregulated_genes.txt",
                  "   - L1_EnrichR_Upregulated_DSigDB.csv",
                  "   - L1_EnrichR_Downregulated_DSigDB.csv",
                  "   - L1_Dual_Action_Drugs_DSigDB.csv",
                  "\n   Layer 2 (靶向TF):",
                  "   - L2_Drugs_targeting_TFs.csv",
                  "   - L2_ChEA_TF_enrichment.csv",
                  "\n   Layer 3 (miRNA治疗):",
                  "   - L3_miRNA_therapeutic_strategy.csv",
                  "\n   综合推荐:",
                  "   - Comprehensive_Drug_Recommendations.csv",
                  "\n   可视化:",
                  "   - Fig9_miRNA_Therapeutic_Strategy.pdf",
                  "   - Fig10_Comprehensive_Drug_Recommendations.pdf",
                  "   - Fig11_Drug_Landscape_Bubble.pdf",
                  "\n========================================"
)

writeLines(report_lines, file.path(wgcna_dir, "Drug_Screening_Report.txt"))

cat("✅ 药物筛选报告已生成\n\n")

cat("✅✅✅ 多层次药物筛选完成！✅✅✅\n\n")

# ======================================================
# =========     最终总结     =========
# ======================================================
cat(rep("=", 70), "\n", sep="")
cat("🎉🎉🎉 所有分析流程已完成！🎉🎉🎉\n")
cat(rep("=", 70), "\n\n", sep="")

cat("📊 生成的主要图表:\n")
cat("   ✓ Fig1: 模块-性状相关性热图\n")
cat("   ✓ Fig2: DEG与模块基因Venn图\n")
cat("   ✓ Fig3: GO富集分析点图\n")
cat("   ✓ Fig4: KEGG富集分析柱状图\n")
cat("   ✓ Fig5: LASSO交叉验证图\n")
cat("   ✓ Fig6: 模型ROC曲线对比\n")
cat("   ✓ Fig7: 核心基因表达热图\n")
if (exists("tf_target_pairs")) cat("   ✓ Fig8A: TF调控网络图\n")
if (exists("mirna_target_pairs")) cat("   ✓ Fig8B: miRNA调控网络图\n")
if (exists("integrated_edges")) cat("   ✓ Fig8C: 整合调控网络图\n")
if (exists("mirna_therapeutic_strategy")) cat("   ✓ Fig9: miRNA治疗策略图\n")
if (nrow(comprehensive_drugs) > 0) {
  cat("   ✓ Fig10: 综合药物推荐图\n")
  cat("   ✓ Fig11: 药物景观气泡图\n")
}

cat("\n📁 所有结果保存在:\n")
cat("   主目录:", wgcna_dir, "\n")
cat("   调控网络:", regulation_dir, "\n")

cat("\n🔑 核心发现总结:\n")
cat("   - 核心基因数:", length(lasso_genes), "个\n")
cat("   - 核心基因:", paste(lasso_genes, collapse = ", "), "\n")
if (exists("key_TFs")) {
  cat("   - 关键TF数:", length(key_TFs), "个\n")
  cat("   - 关键TF:", paste(head(key_TFs, 10), collapse = ", "), "...\n")
}
if (exists("key_miRNAs")) {
  cat("   - 关键miRNA数:", length(key_miRNAs), "个\n")
  cat("   - 关键miRNA:", paste(head(key_miRNAs, 10), collapse = ", "), "...\n")
}
if (nrow(comprehensive_drugs) > 0) {
  cat("   - 候选药物/策略总数:", nrow(comprehensive_drugs), "个\n")
  cat("   - Top 3 推荐:\n")
  for (i in 1:min(3, nrow(comprehensive_drugs))) {
    cat(paste0("     ", i, ". ", comprehensive_drugs$Drug_Name[i], 
               " (", comprehensive_drugs$Target_Level[i], ")\n"))
  }
}

cat("\n💡 优化后的逻辑优势:\n")
cat("   ✅ 先识别调控网络，后筛选药物\n")
cat("   ✅ 三层药物筛选（基因 → TF → miRNA）\n")
cat("   ✅ 综合推荐，优先级明确\n")
cat("   ✅ 短中长期策略清晰\n")

cat("\n📖 推荐阅读报告:\n")
cat("   1. ", file.path(regulation_dir, "Regulatory_Network_Report.txt"), "\n")
cat("   2. ", file.path(wgcna_dir, "Drug_Screening_Report.txt"), "\n")

cat("\n💊 临床转化路径:\n")
cat("   短期(1-2年): Layer 1 双向调控小分子药物\n")
cat("   中期(3-5年): Layer 2 TF靶向药物开发\n")
cat("   长期(5-10年): Layer 3 miRNA基因治疗\n")

cat("\n✨ 分析完成！祝您科研顺利，论文发表成功！✨\n\n")
cat(rep("=", 70), "\n", sep="")
# ======================================================
# =========  第十部分: 药物安全性与成药性评估  =========
# ======================================================

cat("\n")
cat("========================================\n")
cat("  第十部分: 药物安全性与成药性评估\n")
cat("========================================\n\n")

# 明确加载包并解决命名空间冲突
library(httr)
library(jsonlite)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(tidyr)

# 创建药物安全性评估目录
safety_dir <- file.path(wgcna_dir, "Drug_Safety_Assessment")
if (!dir.exists(safety_dir)) {
  dir.create(safety_dir)
}

# ============================================
# 10.1 SIDER 副作用数据库分析
# ============================================
cat("========================================\n")
cat("  Step 1: SIDER 副作用筛选\n")
cat("========================================\n\n")

# SIDER数据库文件路径
sider_drug_se_file <- file.path(safety_dir, "meddra_all_se.tsv")
sider_freq_file <- file.path(safety_dir, "meddra_freq.tsv")
sider_names_file <- file.path(safety_dir, "drug_names.tsv")

# 提示用户手动下载
if (!file.exists(sider_drug_se_file) || !file.exists(sider_names_file)) {
  cat("⚠️ 请手动下载SIDER数据库文件:\n")
  cat("   访问: http://sideeffects.embl.de/download/\n")
  cat("   下载以下文件到:", safety_dir, "\n")
  cat("   1. meddra_all_se.tsv.gz (解压后重命名为 meddra_all_se.tsv)\n")
  cat("   2. meddra_freq.tsv.gz (解压后重命名为 meddra_freq.tsv)\n")
  cat("   3. drug_names.tsv\n\n")
  cat("   💡 如不使用SIDER，将采用备用评估方法\n\n")
}

# 读取SIDER数据
sider_data_available <- FALSE

if (file.exists(sider_drug_se_file) && file.exists(sider_names_file)) {
  
  cat("📖 读取SIDER数据库...\n")
  
  tryCatch({
    # 读取药物名称映射
    sider_names <- read.table(sider_names_file, sep = "\t", header = FALSE,
                              quote = "", stringsAsFactors = FALSE,
                              col.names = c("STITCH_flat", "STITCH_stereo", "Drug_Name"))
    
    # 读取药物-副作用关系
    sider_se <- read.table(sider_drug_se_file, sep = "\t", header = FALSE,
                           quote = "", stringsAsFactors = FALSE,
                           col.names = c("STITCH_flat", "STITCH_stereo",
                                         "UMLS_ID", "MedDRA_Type", 
                                         "UMLS_ID_Label", "Side_Effect"))
    
    # 合并药物名称（使用dplyr::）
    sider_se_named <- sider_se %>%
      dplyr::left_join(sider_names %>% 
                         dplyr::select(STITCH_flat, Drug_Name) %>%
                         dplyr::distinct(),
                       by = "STITCH_flat")
    
    cat(paste0("✅ 成功读取 ", nrow(sider_se_named), " 条副作用记录\n"))
    cat(paste0("✅ 涉及 ", length(unique(sider_se_named$Drug_Name)), " 个药物\n\n"))
    
    sider_data_available <- TRUE
    
  }, error = function(e) {
    cat("⚠️ SIDER数据读取失败:", e$message, "\n")
    cat("   将使用备用评估方法\n\n")
    sider_data_available <<- FALSE
  })
}

# ============================================
# 10.1.1 匹配候选药物的副作用
# ============================================

drug_safety_profile <- data.frame()

if (sider_data_available && exists("comprehensive_drugs")) {
  
  cat("🔍 匹配候选药物副作用信息...\n")
  
  # 提取药物名称
  candidate_drugs <- comprehensive_drugs %>%
    dplyr::filter(Target_Level %in% c("Direct_Gene", "Transcription_Factor")) %>%
    dplyr::mutate(Drug_Clean = tolower(trimws(gsub("-.*", "", Drug_Name))))
  
  for (i in 1:nrow(candidate_drugs)) {
    drug <- candidate_drugs$Drug_Clean[i]
    
    # 模糊匹配药物名称
    matches <- sider_se_named %>%
      dplyr::filter(grepl(drug, tolower(Drug_Name), fixed = TRUE))
    
    if (nrow(matches) > 0) {
      # 统计副作用
      se_summary <- matches %>%
        dplyr::group_by(Drug_Name) %>%
        dplyr::summarise(
          n_side_effects = n(),
          side_effects_list = paste(unique(Side_Effect), collapse = "; "),
          .groups = "drop"
        ) %>%
        head(1)
      
      # 严重副作用关键词
      severe_keywords <- c("death", "fatal", "cancer", "hepatotoxicity",
                           "cardiotoxicity", "neurotoxicity", "nephrotoxicity",
                           "anaphylaxis", "Stevens-Johnson")
      
      severe_se <- matches %>%
        dplyr::filter(grepl(paste(severe_keywords, collapse = "|"), 
                            Side_Effect, ignore.case = TRUE)) %>%
        nrow()
      
      drug_safety_profile <- rbind(
        drug_safety_profile,
        data.frame(
          Original_Drug = candidate_drugs$Drug_Name[i],
          Matched_Drug = se_summary$Drug_Name,
          Total_Side_Effects = se_summary$n_side_effects,
          Severe_Side_Effects = severe_se,
          Side_Effects_Summary = substr(se_summary$side_effects_list, 1, 200),
          Safety_Score = 100 - (se_summary$n_side_effects * 0.5 + severe_se * 5),
          Layer = candidate_drugs$Layer[i]
        )
      )
    }
  }
  
  if (nrow(drug_safety_profile) > 0) {
    # 确保安全评分在0-100之间
    drug_safety_profile$Safety_Score <- pmax(0, pmin(100, drug_safety_profile$Safety_Score))
    
    # 添加安全等级
    drug_safety_profile <- drug_safety_profile %>%
      dplyr::mutate(
        Safety_Grade = dplyr::case_when(
          Safety_Score >= 80 ~ "A (Excellent)",
          Safety_Score >= 60 ~ "B (Good)",
          Safety_Score >= 40 ~ "C (Moderate)",
          Safety_Score >= 20 ~ "D (Poor)",
          TRUE ~ "E (High Risk)"
        )
      ) %>%
      dplyr::arrange(desc(Safety_Score))
    
    write.csv(drug_safety_profile,
              file.path(safety_dir, "Drug_Safety_Profile_SIDER.csv"),
              row.names = FALSE)
    
    cat(paste0("✅ 成功匹配 ", nrow(drug_safety_profile), " 个药物的副作用信息\n"))
    cat("\nTop 10 最安全的候选药物:\n")
    print(drug_safety_profile[1:min(10, nrow(drug_safety_profile)), 
                              c("Original_Drug", "Total_Side_Effects", 
                                "Severe_Side_Effects", "Safety_Grade")])
    
    # 📊 【论文图12】副作用分析图
    p_safety <- drug_safety_profile %>%
      head(25) %>%
      ggplot(aes(x = reorder(Original_Drug, Safety_Score),
                 y = Safety_Score,
                 fill = Safety_Grade)) +
      geom_bar(stat = "identity", alpha = 0.8) +
      geom_text(aes(label = Total_Side_Effects), 
                hjust = -0.3, size = 3) +
      scale_fill_manual(
        values = c("A (Excellent)" = "#27ae60",
                   "B (Good)" = "#f39c12",
                   "C (Moderate)" = "#e67e22",
                   "D (Poor)" = "#e74c3c",
                   "E (High Risk)" = "#c0392b"),
        name = "Safety Grade"
      ) +
      coord_flip() +
      labs(
        title = "Drug Safety Profile Based on SIDER Database",
        subtitle = "Numbers indicate total reported side effects",
        x = "Drug Name",
        y = "Safety Score (0-100)"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "bottom"
      )
    
    ggsave(file.path(safety_dir, "Fig12_Drug_Safety_Profile.pdf"),
           plot = p_safety, width = 12, height = 10)
    
    cat("\n✅ Fig12: 药物安全性分析图已生成\n\n")
  }
  
} else if (!sider_data_available) {
  cat("⚠️ SIDER数据未加载，使用备用方法...\n\n")
  
  # 备用方法：基于文献的副作用风险评估
  if (exists("comprehensive_drugs")) {
    drug_safety_profile <- comprehensive_drugs %>%
      dplyr::filter(Target_Level %in% c("Direct_Gene", "Transcription_Factor")) %>%
      dplyr::mutate(
        # 基于药物类型的启发式安全评分
        Safety_Score = dplyr::case_when(
          grepl("metformin|aspirin|vitamin", Drug_Name, ignore.case = TRUE) ~ 85,
          grepl("statin|inhibitor", Drug_Name, ignore.case = TRUE) ~ 70,
          grepl("chemotherapy|cytotoxic", Drug_Name, ignore.case = TRUE) ~ 40,
          TRUE ~ 60
        ),
        Safety_Grade = dplyr::case_when(
          Safety_Score >= 80 ~ "A (Predicted Good)",
          Safety_Score >= 60 ~ "B (Predicted Moderate)",
          TRUE ~ "C (Needs Evaluation)"
        ),
        Data_Source = "Heuristic Estimation"
      )
    
    # 使用显式列选择
    drug_safety_profile <- data.frame(
      Original_Drug = drug_safety_profile$Drug_Name,
      Safety_Score = drug_safety_profile$Safety_Score,
      Safety_Grade = drug_safety_profile$Safety_Grade,
      Layer = drug_safety_profile$Layer,
      Data_Source = drug_safety_profile$Data_Source
    )
    
    write.csv(drug_safety_profile,
              file.path(safety_dir, "Drug_Safety_Profile_Estimated.csv"),
              row.names = FALSE)
    
    cat("✅ 生成基于启发式的安全评分\n")
    cat(paste0("   评估了 ", nrow(drug_safety_profile), " 个药物\n"))
    cat("   建议：使用PubChem或DrugBank API进一步验证\n\n")
  }
}

# ============================================
# 10.2 ADMET 性质预测
# ============================================
cat("========================================\n")
cat("  Step 2: ADMET 性质预测\n")
cat("========================================\n\n")

# 检查并安装必要的包
if (!require("webchem", quietly = TRUE)) {
  cat("📦 安装 webchem 包...\n")
  tryCatch({
    install.packages("webchem", dependencies = TRUE)
    library(webchem)
  }, error = function(e) {
    cat("⚠️ webchem 安装失败，将跳过PubChem查询\n")
  })
}

# ============================================
# 10.2.1 从PubChem获取化合物信息
# ============================================

get_pubchem_properties <- function(drug_name) {
  
  if (!requireNamespace("webchem", quietly = TRUE)) {
    return(NULL)
  }
  
  tryCatch({
    # 搜索化合物CID
    cid_result <- webchem::get_cid(drug_name, from = "name", match = "first")
    
    if (is.na(cid_result[[1]])) {
      return(NULL)
    }
    
    cid <- cid_result[[1]]
    
    # 获取化合物性质
    props <- webchem::pc_prop(cid, properties = c("MolecularWeight", 
                                                  "XLogP", 
                                                  "HBondDonorCount",
                                                  "HBondAcceptorCount",
                                                  "RotatableBondCount",
                                                  "TPSA"))
    
    if (!is.null(props) && nrow(props) > 0) {
      return(data.frame(
        Drug = drug_name,
        CID = cid,
        MW = as.numeric(props$MolecularWeight),
        LogP = as.numeric(props$XLogP),
        HBD = as.numeric(props$HBondDonorCount),
        HBA = as.numeric(props$HBondAcceptorCount),
        RotBonds = as.numeric(props$RotatableBondCount),
        TPSA = as.numeric(props$TPSA)
      ))
    }
    
    return(NULL)
    
  }, error = function(e) {
    return(NULL)
  })
}

# ============================================
# 10.2.2 计算 Lipinski's Rule of Five
# ============================================

calculate_lipinski <- function(mw, logp, hbd, hba) {
  violations <- 0
  
  if (!is.na(mw) && mw > 500) violations <- violations + 1
  if (!is.na(logp) && logp > 5) violations <- violations + 1
  if (!is.na(hbd) && hbd > 5) violations <- violations + 1
  if (!is.na(hba) && hba > 10) violations <- violations + 1
  
  return(violations)
}

# ============================================
# 10.2.3 批量获取候选药物的ADMET性质
# ============================================

admet_results <- data.frame()

if (exists("comprehensive_drugs") && requireNamespace("webchem", quietly = TRUE)) {
  
  # 提取小分子药物
  small_molecule_drugs <- comprehensive_drugs %>%
    dplyr::filter(Target_Level %in% c("Direct_Gene", "Transcription_Factor")) %>%
    dplyr::mutate(Drug_Clean = tolower(trimws(gsub("-.*", "", Drug_Name)))) %>%
    dplyr::distinct(Drug_Clean, .keep_all = TRUE) %>%
    head(20)  # 限制数量以避免API限流
  
  cat(paste0("📊 分析 ", nrow(small_molecule_drugs), " 个候选药物的ADMET性质...\n"))
  cat("   (这可能需要几分钟...)\n\n")
  
  pb <- txtProgressBar(min = 0, max = nrow(small_molecule_drugs), style = 3)
  
  for (i in 1:nrow(small_molecule_drugs)) {
    drug <- small_molecule_drugs$Drug_Clean[i]
    
    # 获取PubChem数据
    props <- get_pubchem_properties(drug)
    
    if (!is.null(props)) {
      # 计算Lipinski违规
      lipinski_violations <- calculate_lipinski(
        props$MW, props$LogP, props$HBD, props$HBA
      )
      
      # ADMET评分
      admet_score <- 100
      
      if (!is.na(props$MW)) {
        if (props$MW < 200 || props$MW > 500) admet_score <- admet_score - 15
      }
      
      if (!is.na(props$LogP)) {
        if (props$LogP < 0 || props$LogP > 3) admet_score <- admet_score - 10
      }
      
      if (!is.na(props$TPSA)) {
        if (props$TPSA > 140) admet_score <- admet_score - 15
      }
      
      if (!is.na(props$RotBonds)) {
        if (props$RotBonds > 10) admet_score <- admet_score - 10
      }
      
      admet_score <- admet_score - (lipinski_violations * 10)
      admet_score <- max(0, min(100, admet_score))
      
      # 药物等级
      drug_likeness <- dplyr::case_when(
        lipinski_violations == 0 ~ "Excellent",
        lipinski_violations == 1 ~ "Good",
        lipinski_violations == 2 ~ "Moderate",
        TRUE ~ "Poor"
      )
      
      # BBB渗透预测
      bbb_penetration <- "Unknown"
      if (!is.na(props$LogP) && !is.na(props$TPSA)) {
        if (props$LogP > 1 && props$LogP < 3 && props$TPSA < 90) {
          bbb_penetration <- "High"
        } else if (props$TPSA < 140) {
          bbb_penetration <- "Moderate"
        } else {
          bbb_penetration <- "Low"
        }
      }
      
      # 口服生物利用度预测
      oral_bioavailability <- "Unknown"
      if (!is.na(props$TPSA) && !is.na(props$RotBonds)) {
        if (props$TPSA < 140 && props$RotBonds < 10) {
          oral_bioavailability <- "Good"
        } else {
          oral_bioavailability <- "Poor"
        }
      }
      
      admet_results <- rbind(
        admet_results,
        data.frame(
          Original_Drug = small_molecule_drugs$Drug_Name[i],
          Drug_Clean = drug,
          CID = props$CID,
          MW = props$MW,
          LogP = props$LogP,
          HBD = props$HBD,
          HBA = props$HBA,
          RotBonds = props$RotBonds,
          TPSA = props$TPSA,
          Lipinski_Violations = lipinski_violations,
          Drug_Likeness = drug_likeness,
          ADMET_Score = admet_score,
          BBB_Penetration = bbb_penetration,
          Oral_Bioavailability = oral_bioavailability,
          Layer = small_molecule_drugs$Layer[i]
        )
      )
    }
    
    setTxtProgressBar(pb, i)
    Sys.sleep(0.5)  # 避免API限流
  }
  
  close(pb)
  
  if (nrow(admet_results) > 0) {
    admet_results <- admet_results %>%
      dplyr::arrange(desc(ADMET_Score))
    
    write.csv(admet_results,
              file.path(safety_dir, "ADMET_Analysis_Results.csv"),
              row.names = FALSE)
    
    cat(paste0("\n✅ 成功获取 ", nrow(admet_results), " 个药物的ADMET性质\n"))
    cat("\nTop 10 最佳成药性候选:\n")
    print(admet_results[1:min(10, nrow(admet_results)), 
                        c("Original_Drug", "MW", "LogP", "Lipinski_Violations",
                          "Drug_Likeness", "ADMET_Score")])
    
    # 📊 【论文图13】Lipinski规则散点图
    library(ggrepel)
    
    p_lipinski <- ggplot(admet_results, 
                         aes(x = MW, y = LogP, 
                             color = Drug_Likeness,
                             size = ADMET_Score)) +
      geom_point(alpha = 0.7) +
      geom_vline(xintercept = 500, linetype = "dashed", color = "red") +
      geom_hline(yintercept = 5, linetype = "dashed", color = "red") +
      geom_text_repel(aes(label = Original_Drug), size = 3, max.overlaps = 15) +
      scale_color_manual(
        values = c("Excellent" = "#27ae60",
                   "Good" = "#f39c12",
                   "Moderate" = "#e67e22",
                   "Poor" = "#e74c3c"),
        name = "Drug-likeness"
      ) +
      scale_size_continuous(range = c(3, 10), name = "ADMET Score") +
      labs(
        title = "Lipinski's Rule of Five Analysis",
        subtitle = "MW ≤ 500, LogP ≤ 5 (red lines)",
        x = "Molecular Weight (Da)",
        y = "LogP (Lipophilicity)"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "right"
      )
    
    ggsave(file.path(safety_dir, "Fig13_Lipinski_Rule_Scatter.pdf"),
           plot = p_lipinski, width = 12, height = 10)
    
    cat("\n✅ Fig13: Lipinski规则散点图已生成\n\n")
  }
  
} else {
  cat("⚠️ webchem包不可用或无候选药物，跳过ADMET分析\n\n")
}

# ============================================
# 10.3 整合安全性与成药性评分
# ============================================
cat("========================================\n")
cat("  Step 3: 整合安全性与成药性评估\n")
cat("========================================\n\n")

final_ranking <- data.frame()

if (exists("comprehensive_drugs")) {
  
  final_ranking <- comprehensive_drugs %>%
    dplyr::filter(Target_Level %in% c("Direct_Gene", "Transcription_Factor"))
  
  # 合并安全性评分
  if (exists("drug_safety_profile") && nrow(drug_safety_profile) > 0) {
    safety_subset <- data.frame(
      Original_Drug = drug_safety_profile$Original_Drug,
      Safety_Score = drug_safety_profile$Safety_Score,
      Safety_Grade = drug_safety_profile$Safety_Grade
    )
    
    final_ranking <- final_ranking %>%
      dplyr::left_join(safety_subset,
                       by = c("Drug_Name" = "Original_Drug"))
  } else {
    final_ranking$Safety_Score <- NA
    final_ranking$Safety_Grade <- NA
  }
  
  # 合并ADMET评分
  if (exists("admet_results") && nrow(admet_results) > 0) {
    admet_subset <- data.frame(
      Original_Drug = admet_results$Original_Drug,
      ADMET_Score = admet_results$ADMET_Score,
      Drug_Likeness = admet_results$Drug_Likeness,
      Lipinski_Violations = admet_results$Lipinski_Violations,
      Oral_Bioavailability = admet_results$Oral_Bioavailability
    )
    
    final_ranking <- final_ranking %>%
      dplyr::left_join(admet_subset,
                       by = c("Drug_Name" = "Original_Drug"))
  } else {
    final_ranking$ADMET_Score <- NA
    final_ranking$Drug_Likeness <- NA
    final_ranking$Lipinski_Violations <- NA
    final_ranking$Oral_Bioavailability <- NA
  }
  
  # 填充缺失值
  final_ranking$Safety_Score[is.na(final_ranking$Safety_Score)] <- 60
  final_ranking$ADMET_Score[is.na(final_ranking$ADMET_Score)] <- 60
  
  # 计算综合评分
  final_ranking <- final_ranking %>%
    dplyr::mutate(
      Efficacy_Score = Priority_Score / max(Priority_Score, na.rm = TRUE) * 100,
      Comprehensive_Score = 
        Efficacy_Score * 0.4 + 
        Safety_Score * 0.3 + 
        ADMET_Score * 0.3,
      
      Overall_Grade = dplyr::case_when(
        Comprehensive_Score >= 80 ~ "★★★★★ (Highly Recommended)",
        Comprehensive_Score >= 70 ~ "★★★★☆ (Recommended)",
        Comprehensive_Score >= 60 ~ "★★★☆☆ (Potential)",
        Comprehensive_Score >= 50 ~ "★★☆☆☆ (Limited)",
        TRUE ~ "★☆☆☆☆ (Not Recommended)"
      ),
      
      Clinical_Priority = dplyr::case_when(
        Comprehensive_Score >= 75 & !is.na(ADMET_Score) ~ "High Priority",
        Comprehensive_Score >= 60 ~ "Medium Priority",
        TRUE ~ "Low Priority"
      )
    ) %>%
    dplyr::arrange(desc(Comprehensive_Score))
  
  write.csv(final_ranking,
            file.path(safety_dir, "Final_Drug_Ranking_Integrated.csv"),
            row.names = FALSE)
  
  cat("✅ 综合评估完成！\n\n")
  cat("🏆 Top 15 最终推荐候选药物:\n")
  cat(rep("=", 80), "\n", sep = "")
  
  for (i in 1:min(15, nrow(final_ranking))) {
    drug_info <- final_ranking[i, ]
    cat(paste0(i, ". ", drug_info$Drug_Name, " ", drug_info$Overall_Grade, "\n"))
    cat(paste0("   疗效评分: ", round(drug_info$Efficacy_Score, 1), 
               " | 安全性: ", round(drug_info$Safety_Score, 1),
               " | 成药性: ", round(drug_info$ADMET_Score, 1), "\n"))
    cat(paste0("   综合评分: ", round(drug_info$Comprehensive_Score, 1), "\n"))
    cat(paste0("   临床优先级: ", drug_info$Clinical_Priority, "\n"))
    cat(paste0("   作用层级: ", drug_info$Target_Level, 
               " (", drug_info$Layer, ")\n\n"))
  }
  
  # 📊 【论文图14】综合评分对比图
  plot_data <- final_ranking %>%
    head(20) %>%
    dplyr::select(Drug_Name, Efficacy_Score, Safety_Score, ADMET_Score) %>%
    tidyr::pivot_longer(cols = c(Efficacy_Score, Safety_Score, ADMET_Score),
                        names_to = "Score_Type", values_to = "Score") %>%
    dplyr::mutate(Score_Type = factor(Score_Type, 
                                      levels = c("Efficacy_Score", "Safety_Score", "ADMET_Score"),
                                      labels = c("Efficacy", "Safety", "ADMET")))
  
  p_comprehensive <- ggplot(plot_data, 
                            aes(x = reorder(Drug_Name, Score),
                                y = Score,
                                fill = Score_Type)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    scale_fill_manual(
      values = c("Efficacy" = "#3498db",
                 "Safety" = "#2ecc71",
                 "ADMET" = "#e74c3c"),
      name = "Evaluation Dimension"
    ) +
    coord_flip() +
    labs(
      title = "Comprehensive Drug Evaluation (Top 20)",
      subtitle = "Three-dimensional assessment: Efficacy + Safety + ADMET",
      x = "Drug Name",
      y = "Score (0-100)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "bottom"
    )
  
  ggsave(file.path(safety_dir, "Fig14_Comprehensive_Drug_Evaluation.pdf"),
         plot = p_comprehensive, width = 12, height = 12)
  
  cat("✅ Fig14: 综合评估对比图已生成\n\n")
}

# ============================================
# 10.4 生成详细报告
# ============================================
cat("\n📝 生成药物安全性与成药性评估报告...\n")

safety_report_lines <- c(
  "========================================",
  "   药物安全性与成药性评估报告",
  "========================================",
  paste("\n评估日期:", Sys.Date()),
  paste("候选药物总数:", ifelse(exists("comprehensive_drugs"), nrow(comprehensive_drugs), 0)),
  paste("完成安全性评估:", ifelse(exists("drug_safety_profile"), nrow(drug_safety_profile), 0), "个"),
  paste("完成ADMET评估:", ifelse(exists("admet_results"), nrow(admet_results), 0), "个"),
  "\n----------------------------------------",
  "\n[1] 评估方法",
  "   - 安全性: SIDER数据库 / 启发式评估",
  "   - 成药性: PubChem + Lipinski规则",
  "   - 综合评分 = 疗效40% + 安全性30% + 成药性30%",
  "\n[2] Top 10 综合推荐"
)

if (exists("final_ranking") && nrow(final_ranking) > 0) {
  for (i in 1:min(10, nrow(final_ranking))) {
    safety_report_lines <- c(safety_report_lines,
                             paste0("\n   ", i, ". ", final_ranking$Drug_Name[i]),
                             paste0("      综合评分: ", round(final_ranking$Comprehensive_Score[i], 1)),
                             paste0("      临床优先级: ", final_ranking$Clinical_Priority[i]))
  }
}

safety_report_lines <- c(safety_report_lines,
                         "\n========================================")

writeLines(safety_report_lines, 
           file.path(safety_dir, "Drug_Safety_ADMET_Report.txt"))

cat("✅ 报告已生成\n\n")

cat("✅✅✅ 第十部分完成！✅✅✅\n\n")