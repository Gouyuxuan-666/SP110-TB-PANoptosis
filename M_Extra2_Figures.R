###############################################################################
# SP110 Extra Figures v2 — 多方法免疫 + 亚型 + LASSO + 多ROC对比 + 相关性网络
###############################################################################
setwd("F:/SP110_project")
suppressMessages({
  library(ggplot2); library(pheatmap); library(reshape2); library(caret)
  library(glmnet); library(pROC); library(dplyr); library(RColorBrewer)
  library(IOBR); library(data.table)
})

pdf_dir <- "figures"; dir.create(pdf_dir, showWarnings = FALSE)
load("M1_output/merged_data.RData"); expr <- merged_combat
load("M2_output/wgcna_results.RData"); load("M3_output/ml_results.RData")
load("M4_output/immune_results.RData")
er_genes <- readLines("M1_output/ER_stress_genes.txt")
sp110_expr <- as.numeric(expr["SP110", ])
tb_indicator <- as.numeric(type_vec == "Treat")

# ==========================================
# Fig Ex2-1: Multi-method immune deconvolution comparison
# ==========================================
cat("Ex2-1: Multi-method immune comparison...\n")
# CIBERSORT already done (via IOBR), add EPIC/MCP-counter/xCell
imm_methods <- list()
tryCatch({
  imm_methods$EPIC <- deconvo_epic(as.data.frame(expr), tumor = FALSE)
  cat("  EPIC done\n")
}, error = function(e) cat("  EPIC skipped\n"))
tryCatch({
  imm_methods$MCP <- deconvo_mcp(as.data.frame(expr))
  cat("  MCP done\n")
}, error = function(e) cat("  MCP skipped\n"))

# Compare CIBERSORT vs MCP/EPIC for key cell types
key_types <- names(head(sort(abs(sp110_cor_imm), decreasing = TRUE), 8))
cib_compare <- data.frame(CellType = names(sp110_cor_imm), CIBERSORT = sp110_cor_imm)

if (length(imm_methods) >= 1) {
  cat("  Computing multi-method correlations...\n")
  # EPIC
  m <- imm_methods[["EPIC"]]
  if (!is.null(m)) {
    epic_frac <- as.matrix(m[, grep("_EPIC$", colnames(m), value = TRUE)])
    if (ncol(epic_frac) > 0) {
      epic_cor <- sapply(1:ncol(epic_frac), function(i) cor(sp110_expr, epic_frac[,i], method="spearman", use="pairwise.complete.obs"))
      names(epic_cor) <- gsub("_EPIC$", "", colnames(epic_frac))
      write.csv(data.frame(CellType = names(epic_cor), EPIC_rho = epic_cor),
                file.path(pdf_dir, "Table_EPIC_SP110_corr.csv"), row.names = FALSE)
    }
  }
}

# ==========================================
# Fig Ex2-2: Consensus clustering — ER stress subtypes
# ==========================================
cat("Ex2-2: ER stress subtypes...\n")
if (!require("ConsensusClusterPlus", quietly = TRUE)) {
  cat("  ConsensusClusterPlus not installed. Skipping.\n")
  er_genes_use <- character(0)
} else {
  er_genes_use <- intersect(er_genes, rownames(expr))
}
er_mat <- t(scale(t(expr[er_genes_use, ])))
er_mat <- er_mat[, apply(er_mat, 2, function(x) all(!is.na(x)))]

if (nrow(er_mat) >= 10 && ncol(er_mat) >= 100) {
  tryCatch({
    cc <- ConsensusClusterPlus(er_mat[1:min(50, nrow(er_mat)), ], maxK = 4, reps = 50,
                                pItem = 0.8, pFeature = 1, clusterAlg = "hc",
                                distance = "spearman", seed = 20240717,
                                plot = "pdf", title = file.path(pdf_dir, "Fig_Ex2_Consensus"))
    er_clusters <- cc[[2]]$consensusClass
    cat(sprintf("  Found %d ER stress subtypes\n", length(unique(er_clusters))))
  }, error = function(e) cat("  ConsensusCluster failed\n"))
} else {
  cat("  Skipping: not enough genes/samples\n")
}

# ==========================================
# Fig Ex2-3: LASSO coefficient path + CV
# ==========================================
cat("Ex2-3: LASSO path...\n")
ml_features <- intersect(unique(c(names(head(sort(abs(apply(expr, 1, function(x) cor(sp110_expr, x, method="spearman"))), decreasing=TRUE), 200)), er_genes)), rownames(expr))
x <- t(expr[ml_features, ])
y <- factor(ifelse(type_vec == "Treat", 1, 0))

cvfit <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 5)

pdf(file.path(pdf_dir, "Fig_Ex2_LASSO_CV.pdf"), width = 8, height = 5)
plot(cvfit)
dev.off()

pdf(file.path(pdf_dir, "Fig_Ex2_LASSO_Path.pdf"), width = 8, height = 6)
plot(cvfit$glmnet.fit, xvar = "lambda", label = TRUE)
abline(v = log(cvfit$lambda.min), lty = 2, col = "#E74C3C")
dev.off()

# Selected genes
lasso_coef <- coef(cvfit, s = "lambda.min")
selected <- rownames(lasso_coef)[which(lasso_coef[,1] != 0)][-1]
cat(sprintf("  LASSO selected %d genes\n", length(selected)))

# ==========================================
# Fig Ex2-4: Multi-ROC comparison (SP110 alone vs GBP5 vs model vs random)
# ==========================================
cat("Ex2-4: Multi-ROC...\n")
# SP110 single-gene
sp110_glm <- glm(factor(type_vec) ~ sp110_expr, family = binomial)
sp110_pred <- predict(sp110_glm, type = "response")
sp110_roc <- roc(type_vec == "Treat", sp110_pred, quiet = TRUE)

# GBP5 single-gene
if ("GBP5" %in% rownames(expr)) {
  gbp5_pred <- predict(glm(factor(type_vec) ~ as.numeric(expr["GBP5",]), family = binomial), type = "response")
  gbp5_roc <- roc(type_vec == "Treat", gbp5_pred, quiet = TRUE)
}

# Random baseline
set.seed(42); rand_roc <- roc(type_vec == "Treat", runif(length(type_vec)), quiet = TRUE)

pdf(file.path(pdf_dir, "Fig_Ex2_MultiROC.pdf"), width = 7, height = 7)
plot.roc(sp110_roc, col = "#E74C3C", lwd = 2.5, main = "Diagnostic ROC Comparison")
plot.roc(test_roc, col = "#2E75B6", lwd = 2.5, add = TRUE)
if (exists("gbp5_roc")) plot.roc(gbp5_roc, col = "#F39C12", lwd = 2, lty = 2, add = TRUE)
plot.roc(rand_roc, col = "#BDC3C7", lwd = 1.5, lty = 2, add = TRUE)
legend("bottomright", c(
  sprintf("SP110 Model (AUC=%.3f)", test_auc),
  sprintf("SP110 alone (AUC=%.3f)", as.numeric(sp110_roc$auc)),
  if (exists("gbp5_roc")) sprintf("GBP5 alone (AUC=%.3f)", as.numeric(gbp5_roc$auc)) else NULL,
  sprintf("Random (AUC=%.3f)", as.numeric(rand_roc$auc))
), col = c("#2E75B6", "#E74C3C", "#F39C12", "#BDC3C7"), lwd = c(2.5, 2.5, 2, 1.5), cex = 0.8, bty = "n")
dev.off()

# ==========================================
# Fig Ex2-5: Top SP110-correlated gene network (correlation heatmap)
# ==========================================
cat("Ex2-5: SP110 gene network heatmap...\n")
sp110_cor <- apply(expr, 1, function(x) cor(sp110_expr, x, method="spearman", use="pairwise.complete.obs"))
top30 <- names(head(sort(abs(sp110_cor), decreasing = TRUE), 30))
top30_cor <- cor(t(expr[top30, ]), method = "spearman", use = "pairwise.complete.obs")

pdf(file.path(pdf_dir, "Fig_Ex2_SP110_Network.pdf"), width = 10, height = 9)
pheatmap(top30_cor, show_rownames = TRUE, show_colnames = TRUE, fontsize = 7,
         main = "Top 30 SP110-Correlated Gene Network (Spearman rho)",
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
dev.off()

# ==========================================
# Fig Ex2-6: SP110 expression heatmap (top 50 most variable + SP110)
# ==========================================
cat("Ex2-6: Expression heatmap...\n")
var50 <- names(head(sort(apply(expr, 1, sd), decreasing = TRUE), 50))
var50 <- unique(c("SP110", intersect(c("SP110","SP140","GBP1","GBP5","STAT1","DDIT3","BECN1"), rownames(expr)), var50[1:45]))
heat_mat <- t(scale(t(expr[var50, ])))
heat_mat <- heat_mat[, order(type_vec)]
ann_col <- data.frame(Group = ifelse(type_vec[order(type_vec)] == "Treat", "TB", "Control"),
                       row.names = colnames(heat_mat))

pdf(file.path(pdf_dir, "Fig_Ex2_Heatmap_Top50.pdf"), width = 14, height = 10)
pheatmap(heat_mat, annotation_col = ann_col, show_colnames = FALSE,
         cluster_cols = FALSE, fontsize_row = 7,
         main = "Top Variable Genes + SP110: TB vs Control",
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
dev.off()

# ==========================================
# Fig Ex2-7: SP110 expression by gender/age (if metadata available)
# ==========================================
cat("Ex2-7: SP110 by cohort heatmap...\n")
cohort_means <- sapply(unique(batch_vec), function(b) {
  idx <- batch_vec == b
  c(mean = mean(sp110_expr[idx]), n = sum(idx),
    TB_mean = mean(sp110_expr[idx & type_vec == "Treat"]),
    Ctrl_mean = mean(sp110_expr[idx & type_vec == "Control"]))
})
cohort_df <- as.data.frame(t(cohort_means))
cohort_df$Cohort <- rownames(cohort_df)
cohort_df <- cohort_df[order(cohort_df$mean), ]

pdf(file.path(pdf_dir, "Fig_Ex2_SP110_byCohort.pdf"), width = 8, height = 5)
ggplot(cohort_df, aes(x = reorder(Cohort, mean), y = mean, fill = n)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_gradient(low = "#3498DB", high = "#E74C3C") +
  labs(title = "Mean SP110 Expression by Cohort", x = "", y = "Mean Expression", fill = "n")
dev.off()

cat(sprintf("\nAll Extra2 figures saved to %s/\n", pdf_dir))
