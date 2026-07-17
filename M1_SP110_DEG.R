###############################################################################
# SP110 新题 — M1: 数据整合与差异表达分析
# 在2号机运行: Rscript M1_SP110_DEG.R
# 数据路径: 复用 F:/GBP1_pipeline_2hao/5/*.normalize.txt
# 输出: F:/SP110_project/M1_output/
###############################################################################

options(repos = c(CRAN = "https://mirror.lzu.edu.cn/CRAN/"))
options(warn = 1)

# ---- Setup ----
proj_dir <- "F:/SP110_project"
data_dir <- "F:/GBP1_pipeline_2hao/5"
for (d in c("M1_output","M2_output","M3_output","M4_output","M5_output","M6_output","M7_output","M8_output")) {
  dir.create(file.path(proj_dir, d), showWarnings = FALSE, recursive = TRUE)
}
setwd(proj_dir)

suppressMessages({
  library(limma); library(sva); library(ggplot2); library(ggrepel)
  library(pheatmap); library(data.table); library(metafor)
})

# WGCNA cor wrapper for R 4.6
cor <- function(x, y = NULL, use = "everything",
                method = c("pearson", "kendall", "spearman"), ...) {
  if (use == "p") use <- "pairwise.complete.obs"
  if (use == "na.or.complete") use <- "na.or.complete"
  stats::cor(x, y, use = use, method = match.arg(method))
}

# ---- ER Stress Gene Set ----
er_stress_genes <- c(
  # UPR sensors
  "ERN1","EIF2AK3","ATF6","HSPA5","ERN2",
  # PERK pathway (apoptosis)
  "EIF2S1","ATF4","DDIT3","CASP12","TRIB3","PPP1R15A",
  # IRE1 pathway (autophagy)
  "XBP1","BECN1","MAP1LC3A","MAP1LC3B","SQSTM1","ATG5","ATG7","ATG12","ATG16L1",
  # ERAD
  "HERPUD1","SYVN1","SEL1L","EDEM1","DERL1",
  # Chaperones
  "CALR","PDIA3","PDIA4","PDIA6","HSP90B1","DNAJB9","DNAJC3",
  # Core: SP110 + downstream
  "SP110","SP140","MYBBP1A","RELA",
  # Apoptosis/autophagy crosstalk
  "BCL2","MCL1","BAX","BAK1","BID","CASP3","CASP7","CASP9","CYCS",
  "PARP1","RPS3A","NCL","GBP1","GBP2","GBP5"
)
er_genes_unique <- unique(er_stress_genes)
cat(sprintf("ER stress gene set: %d genes\n", length(er_genes_unique)))

# ---- 1. Merge 6 Discovery Cohorts ----
cat("\n========== Merging 6 cohorts ==========\n")
gse_all <- c("GSE83456","GSE34608","GSE19491","GSE37250","GSE28623","GSE42830")

expr_list <- list(); batch_vec <- c(); type_vec <- c()
for (gse in gse_all) {
  f <- file.path(data_dir, paste0(gse, ".normalize.txt"))
  if (!file.exists(f)) { cat(sprintf("SKIP %s\n", gse)); next }

  dat <- as.data.frame(fread(f, header = TRUE, sep = "\t", check.names = FALSE))
  rownames(dat) <- dat[[1]]; dat[[1]] <- NULL
  dat <- as.matrix(dat); mode(dat) <- "numeric"
  dat <- dat[rowSums(is.na(dat)) == 0, , drop = FALSE]
  cn <- colnames(dat)

  types <- rep("Control", length(cn))
  types[grepl("_Treat$|_treat$|_TB$", cn)] <- "Treat"
  types[grepl("tuberculosis", cn, ignore.case = TRUE)] <- "Treat"
  valid <- !is.na(types)

  dat <- dat[, valid, drop = FALSE]; types <- types[valid]
  cat(sprintf("  %s: %d genes x %d samples (T=%d C=%d)\n", gse, nrow(dat), ncol(dat), sum(types=="Treat"), sum(types=="Control")))

  if (ncol(dat) > 0 && sum(types == "Treat") > 2 && sum(types == "Control") > 2) {
    expr_list[[gse]] <- dat
    batch_vec <- c(batch_vec, rep(gse, ncol(dat)))
    type_vec <- c(type_vec, types)
  }
}

common_genes <- Reduce(intersect, lapply(expr_list, rownames))
cat(sprintf("Common genes: %d\n", length(common_genes)))

merged_expr <- do.call(cbind, lapply(expr_list, function(x) x[common_genes, , drop = FALSE]))
cat(sprintf("Merged: %d genes x %d samples\n", nrow(merged_expr), ncol(merged_expr)))

# ComBat
mod <- model.matrix(~ type_vec)
merged_combat <- ComBat(dat = merged_expr, batch = batch_vec, mod = mod)
save(merged_expr, merged_combat, batch_vec, type_vec, er_genes_unique, file = "M1_output/merged_data.RData")

# ---- 2. DEG Analysis ----
cat("\n========== DEG Analysis ==========\n")
expr <- merged_combat
meta <- data.frame(condition = type_vec, batch = batch_vec)
design <- model.matrix(~ condition + batch, data = meta)
fit <- lmFit(expr, design); fit <- eBayes(fit)
deg_tab <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")

logFCcut <- 0.585; pcut <- 0.05
degs <- deg_tab[abs(deg_tab$logFC) > logFCcut & deg_tab$adj.P.Val < pcut, ]
cat(sprintf("DEGs: %d total (%d up, %d down)\n", nrow(degs), sum(degs$logFC > 0), sum(degs$logFC < 0)))

# SP110
if ("SP110" %in% rownames(deg_tab)) {
  cat(sprintf("SP110: logFC=%.3f, adjP=%.2e\n", deg_tab["SP110","logFC"], deg_tab["SP110","adj.P.Val"]))
}

# ER stress genes in DEGs
er_in_degs <- intersect(er_genes_unique, rownames(degs))
cat(sprintf("ER stress genes in DEGs: %d\n", length(er_in_degs)))
cat("Top ER stress DEGs:\n")
print(deg_tab[er_in_degs, c("logFC","adj.P.Val")][order(abs(deg_tab[er_in_degs,"logFC"]), decreasing = TRUE)[1:min(20,length(er_in_degs))],])

# ---- 3. SP110 Expression per Cohort ----
cat("\n========== SP110 Expression ==========\n")
for (gse in names(expr_list)) {
  idx <- batch_vec == gse
  if (sum(idx) < 5) next
  sp110_level <- as.numeric(expr["SP110", idx])
  tb_idx <- type_vec[idx] == "Treat"
  t_test <- t.test(sp110_level[tb_idx], sp110_level[!tb_idx])
  cat(sprintf("  %-12s TB=%.3f Ctrl=%.3f logFC=%.3f P=%.4f\n",
              gse, mean(sp110_level[tb_idx]), mean(sp110_level[!tb_idx]),
              mean(sp110_level[tb_idx]) - mean(sp110_level[!tb_idx]), t_test$p.value))
}

# ---- 4. Meta Analysis (SP110) ----
cat("\n========== Meta Analysis: SP110 ==========\n")
meta_df <- data.frame()
for (gse in names(expr_list)) {
  idx <- batch_vec == gse
  if (sum(idx) < 6) next
  sub_expr <- expr[, idx]; sub_type <- type_vec[idx]
  sfit <- eBayes(lmFit(sub_expr, model.matrix(~ sub_type)))
  if ("SP110" %in% rownames(sfit$coefficients)) {
    stt <- topTable(sfit, coef = 2, number = Inf)
    sp <- stt["SP110", ]
    meta_df <- rbind(meta_df, data.frame(
      cohort = gse, n = sum(idx), n_TB = sum(sub_type=="Treat"), n_Control = sum(sub_type=="Control"),
      logFC = sp$logFC, SE = abs(sp$logFC)/abs(sp$t), P = sp$P.Value, stringsAsFactors = FALSE
    ))
  }
}

if (nrow(meta_df) >= 3) {
  ma <- rma(yi = logFC, sei = SE, data = meta_df, method = "REML")
  cat(sprintf("Pooled logFC: %.3f (95%% CI %.3f-%.3f), I2=%.1f%%, P=%.4f\n",
              ma$b, ma$ci.lb, ma$ci.ub, ma$I2, ma$pval))

  pdf("M1_output/Fig_SP110_Forest.pdf", width = 10, height = 5)
  forest(ma, slab = meta_df$cohort, xlab = "SP110 log2 Fold Change",
         header = "Cohort", mlab = "RE Model", cex = 0.9)
  dev.off()
}

# ---- 5. ER Stress Heatmap ----
cat("\n========== ER Stress Heatmap ==========\n")
er_expr <- expr[intersect(er_genes_unique, rownames(expr)), ]
er_expr <- er_expr[rowSums(is.na(er_expr)) == 0, ]
# Scale
er_scaled <- t(scale(t(er_expr)))

annotation <- data.frame(Group = ifelse(type_vec == "Treat", "TB", "Control"),
                          row.names = colnames(er_scaled))
pdf("M1_output/Fig_ER_Stress_Heatmap.pdf", width = 14, height = 10)
pheatmap(er_scaled[, order(type_vec)], annotation_col = annotation,
         show_colnames = FALSE, cluster_cols = FALSE,
         main = "ER Stress Gene Expression (TB vs Control)",
         fontsize_row = 7)
dev.off()

# ---- 6. Save results ----
write.table(deg_tab, "M1_output/DEG_full.txt", sep = "\t", quote = FALSE)
write.table(degs, "M1_output/DEG_filtered.txt", sep = "\t", quote = FALSE)
writeLines(er_genes_unique, "M1_output/ER_stress_genes.txt")

cat("\n========== M1 Complete ==========\n")
cat(sprintf("Output: %s/M1_output/\n", proj_dir))
