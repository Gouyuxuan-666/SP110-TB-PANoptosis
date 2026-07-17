###############################################################################
# SP110 M4 — 免疫微环境分析
###############################################################################
setwd("F:/SP110_project")
suppressMessages({
  library(IOBR); library(ggplot2); library(pheatmap); library(reshape2)
  library(data.table)
})

load("M1_output/merged_data.RData")
expr <- merged_combat

# ---- 1. CIBERSORT ----
cat("Running CIBERSORT...\n")
cib <- deconvo_cibersort(eset = as.data.frame(expr), arrays = TRUE, perm = 100)
cib_pass <- cib$`P-value_CIBERSORT` < 0.05
cat(sprintf("Samples passing QC: %d/%d (%.1f%%)\n", sum(cib_pass), length(cib_pass), 100*mean(cib_pass)))

cib_frac <- cib[, 2:23]  # 22 cell type fractions

# ---- 2. SP110 high vs low ----
sp110_expr <- as.numeric(expr["SP110", ])
sp110_high <- sp110_expr > median(sp110_expr)

cat("\nImmune differences (SP110 high vs low):\n")
for (i in 1:22) {
  ct <- colnames(cib_frac)[i]
  high_vals <- cib_frac[sp110_high & cib_pass, i]
  low_vals <- cib_frac[!sp110_high & cib_pass, i]
  if (length(high_vals) > 3 && length(low_vals) > 3) {
    p <- wilcox.test(high_vals, low_vals)$p.value
    d <- mean(high_vals) - mean(low_vals)
    if (abs(d) > 0.01) {
      cat(sprintf("  %-30s delta=%.3f P=%.4f\n", ct, d, p))
    }
  }
}

# ---- 3. SP110 vs 22 cell types correlation ----
cat("\nSP110 correlation with immune cells (Spearman):\n")
sp110_cor_imm <- sapply(1:22, function(i) {
  cor(sp110_expr[cib_pass], cib_frac[cib_pass, i], method = "spearman", use = "pairwise.complete.obs")
})
names(sp110_cor_imm) <- colnames(cib_frac)
sp110_cor_imm <- sort(sp110_cor_imm, decreasing = TRUE)
for (i in seq_along(sp110_cor_imm)) {
  ct <- names(sp110_cor_imm)[i]
  cat(sprintf("  %-30s rho=%.3f\n", ct, sp110_cor_imm[i]))
}

# ---- 4. ER stress pathway correlations ----
cat("\nSP110 vs immune checkpoints:\n")
icb_genes <- c("CD274","PDCD1","CTLA4","HAVCR2","LAG3","TIGIT")
for (g in intersect(icb_genes, rownames(expr))) {
  ct <- cor.test(sp110_expr, as.numeric(expr[g, ]), method = "spearman")
  cat(sprintf("  SP110 vs %-10s rho=%.3f P=%.2e\n", g, ct$estimate, ct$p.value))
}

# ---- 5. Save ----
save(cib, cib_frac, sp110_cor_imm, sp110_high, cib_pass,
     file = "M4_output/immune_results.RData")
write.csv(data.frame(CellType = names(sp110_cor_imm), Rho = sp110_cor_imm),
          "M4_output/SP110_immune_corr.csv", row.names = FALSE)

cat("\nM4 Complete.\n")
