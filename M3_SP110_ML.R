###############################################################################
# SP110 M3 — 机器学习诊断模型 (SP110-correlated genes as features)
###############################################################################
setwd("F:/SP110_project")
suppressMessages({
  library(caret); library(pROC); library(randomForest); library(glmnet)
  library(e1071); library(gbm); library(mboost); library(data.table)
})

load("M1_output/merged_data.RData")
er_genes <- readLines("M1_output/ER_stress_genes.txt")

# ---- Strategy: SP110-correlated genes + ER stress DEGs as features ----
expr <- merged_combat
cat("Computing SP110-correlated genes...\n")

sp110_expr <- as.numeric(expr["SP110", ])
sp110_cor <- apply(expr, 1, function(x) cor(sp110_expr, x, method = "spearman", use = "pairwise.complete.obs"))
sp110_cor <- sp110_cor[!is.na(sp110_cor)]

# Top 500 SP110-correlated genes
top_cor_genes <- names(sort(abs(sp110_cor), decreasing = TRUE)[1:500])

# Add ER stress DEGs
deg_tab <- fread("M1_output/DEG_full.txt", data.table = FALSE)
rownames(deg_tab) <- deg_tab[[1]]
er_degs <- intersect(er_genes, rownames(deg_tab)[abs(deg_tab$logFC) > 0.3 & deg_tab$adj.P.Val < 0.05])
cat(sprintf("ER stress DEGs (relaxed): %d\n", length(er_degs)))
cat("ER stress DEGs: "); cat(er_degs, sep = ", "); cat("\n")

# Union: top SP110-correlated + ER stress DEGs
ml_features <- union(top_cor_genes, er_degs)
ml_features <- intersect(ml_features, rownames(expr))
cat(sprintf("ML feature pool: %d genes\n", length(ml_features)))

# Get expression matrix
ml_expr <- t(expr[ml_features, ])
ml_labels <- factor(ifelse(type_vec == "Treat", "TB", "Control"))
cat(sprintf("ML data: %d samples x %d features\n", nrow(ml_expr), ncol(ml_expr)))

# ---- Train/test split (70/30) ----
set.seed(20240717)
train_idx <- createDataPartition(ml_labels, p = 0.7, list = FALSE)
train_x <- as.data.frame(ml_expr[train_idx, ])
test_x <- as.data.frame(ml_expr[-train_idx, ])
train_y <- ml_labels[train_idx]; test_y <- ml_labels[-train_idx]

# ---- LDA feature selection ----
cat("LDA feature selection...\n")
lda_subset <- train_x[, 1:min(50, ncol(train_x))]
lda_fit <- MASS::lda(x = lda_subset, grouping = train_y)
lda_scores <- abs(lda_fit$scaling[, 1])
sel_genes <- names(sort(lda_scores, decreasing = TRUE))[1:min(30, length(lda_scores))]
cat(sprintf("Selected %d genes by LDA\n", length(sel_genes)))

# ---- RF classifier ----
cat("Training RF...\n")
ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
rf_model <- caret::train(x = train_x[, sel_genes, drop = FALSE], y = train_y,
                         method = "rf", trControl = ctrl, metric = "ROC", tuneLength = 3)

train_auc <- max(rf_model$results$ROC, na.rm = TRUE)
cat(sprintf("Train AUC: %.4f\n", train_auc))

# ---- Internal test ----
test_pred <- predict(rf_model, newdata = test_x[, sel_genes, drop = FALSE], type = "prob")[, "TB"]
test_roc <- roc(test_y, test_pred, quiet = TRUE)
test_auc <- as.numeric(test_roc$auc)
cat(sprintf("Test AUC: %.4f\n", test_auc))

youden <- coords(test_roc, "best", ret = c("sensitivity", "specificity"))
cat(sprintf("Sens=%.3f Spec=%.3f (Youden)\n", youden[1], youden[2]))

# ---- Feature importance ----
cat("\nTop 10 features:\n")
imp <- varImp(rf_model)$importance
imp <- imp[order(imp$Overall, decreasing = TRUE), , drop = FALSE]
print(head(imp, 10))

# ---- SP110 importance ----
if ("SP110" %in% rownames(imp)) {
  cat(sprintf("SP110 rank: %d/%d\n", which(rownames(imp) == "SP110"), nrow(imp)))
}

# ---- PPV/NPV ----
sens <- as.numeric(youden[1]); spec <- as.numeric(youden[2])
cat("\nPPV/NPV scenarios:\n")
for (p in c(0.02, 0.05, 0.10, 0.20, 0.40)) {
  ppv <- (sens * p) / (sens * p + (1 - spec) * (1 - p))
  npv <- (spec * (1 - p)) / (spec * (1 - p) + (1 - sens) * p)
  cat(sprintf("  Prev=%3d%%  PPV=%.3f  NPV=%.3f\n", p*100, ppv, npv))
}

# ---- Save ----
save(rf_model, test_roc, train_auc, test_auc, sel_genes, youden, imp, ml_features,
     file = "M3_output/ml_results.RData")
write.csv(imp, "M3_output/feature_importance.csv")

cat("\nM3 Complete. Output: M3_output/\n")
