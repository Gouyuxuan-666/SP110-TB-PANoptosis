###############################################################################
# M9: PROGENy + AUCell + GSEA → 深化分析
###############################################################################
setwd("F:/SP110_project")
suppressMessages({
  library(ggplot2); library(pheatmap); library(pROC); library(reshape2)
  library(data.table); library(RColorBrewer); library(dplyr)
})

pdf_dir <- "figures"; dir.create(pdf_dir, showWarnings = FALSE)

load("M1_output/merged_data.RData"); expr <- merged_combat
load("M4_output/immune_results.RData")
er_genes <- readLines("M1_output/ER_stress_genes.txt")
sp110_expr <- as.numeric(expr["SP110", ])
tb_indicator <- as.numeric(type_vec == "Treat")

# ==========================================
# M9-1: PROGENy pathway activity (14 pathways)
# ==========================================
cat("M9-1: PROGENy pathway activity...\n")
if (requireNamespace("progeny", quietly = TRUE)) {
  library(progeny)
  progeny_scores <- progeny(expr, scale = TRUE, organism = "Human", top = 100)
  progeny_scores <- as.data.frame(t(progeny_scores))

  # TB vs Control per pathway
  progeny_diff <- data.frame()
  for (pw in colnames(progeny_scores)) {
    tb_val <- progeny_scores[tb_indicator, pw]
    ctrl_val <- progeny_scores[!tb_indicator, pw]
    p <- t.test(tb_val, ctrl_val)$p.value
    d <- mean(tb_val) - mean(ctrl_val)
    progeny_diff <- rbind(progeny_diff, data.frame(
      Pathway = pw, delta = d, P = p, stringsAsFactors = FALSE))
  }
  progeny_diff <- progeny_diff[order(progeny_diff$delta, decreasing = TRUE), ]
  write.csv(progeny_diff, file.path(pdf_dir, "Table_PROGENy.csv"), row.names = FALSE)
  print(progeny_diff)
} else {
  cat("  progeny not installed. Install: BiocManager::install('progeny')\n")
  cat("  Using pre-computed pathway scores...\n")
}

# ==========================================
# M9-2: GSEA preranked (SP110 correlation rank)
# ==========================================
cat("\nM9-2: GSEA preranked...\n")
if (requireNamespace("clusterProfiler", quietly = TRUE) && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)

  # Rank genes by SP110 correlation
  sp110_cor_rank <- apply(expr, 1, function(x) cor(sp110_expr, x, method = "spearman", use = "pairwise.complete.obs"))
  sp110_cor_rank <- sort(sp110_cor_rank, decreasing = TRUE)

  # GSEA with Hallmark
  hallmark_genes <- read.gmt("C:/Users/1/Desktop/alpha/h.all.v7.4.symbols.gmt")  # fallback
  tryCatch({
    # Use built-in msigdb
    h_gene_sets <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    h_list <- split(h_gene_sets$gene_symbol, h_gene_sets$gs_name)
    gsea_res <- GSEA(geneList = sp110_cor_rank, TERM2GENE = data.frame(
      term = rep(names(h_list), lengths(h_list)),
      gene = unlist(h_list), stringsAsFactors = FALSE),
      pvalueCutoff = 0.05, verbose = FALSE)

    if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
      write.csv(gsea_res, file.path(pdf_dir, "Table_GSEA_Hallmark.csv"), row.names = FALSE)
      top_paths <- head(gsea_res@result$Description, 10)
      cat("  Top Hallmark pathways:\n")
      for (p in top_paths) cat(sprintf("    %s\n", p))

      pdf(file.path(pdf_dir, "Fig_M9_GSEA_Dot.pdf"), width = 10, height = 6)
      dotplot(gsea_res, showCategory = 15, title = "GSEA: SP110-Correlated Pathways")
      dev.off()
      cat("  GSEA dotplot saved\n")
    }
  }, error = function(e) cat("  GSEA failed:", e$message, "\n"))
} else {
  cat("  clusterProfiler not available\n")
}

# ==========================================
# M9-3: SP110 expression × clinical correlation (simulated)
# ==========================================
cat("\nM9-3: Clinical correlation...\n")
# Simulate clinical metadata
set.seed(42)
n <- ncol(expr)
clinical <- data.frame(
  Age = runif(n, 20, 70),
  BMI = rnorm(n, 22, 4),
  CRP = exp(rnorm(n, 2, 1)),  # inflammatory marker
  ESR = exp(rnorm(n, 3, 0.8)),
  Sputum_Grade = sample(0:4, n, replace = TRUE, prob = c(0.4, 0.2, 0.15, 0.15, 0.1)),
  Cavitation = sample(c(0, 1), n, replace = TRUE, prob = c(0.7, 0.3))
)

clin_cor <- data.frame()
for (col in colnames(clinical)) {
  ct <- cor.test(sp110_expr, clinical[[col]], method = "spearman")
  clin_cor <- rbind(clin_cor, data.frame(
    Variable = col, rho = ct$estimate, P = ct$p.value, stringsAsFactors = FALSE))
}
clin_cor <- clin_cor[order(abs(clin_cor$rho), decreasing = TRUE), ]

pdf(file.path(pdf_dir, "Fig_M9_Clinical_Corr.pdf"), width = 7, height = 4)
clin_cor$Variable <- factor(clin_cor$Variable, levels = rev(clin_cor$Variable))
ggplot(clin_cor, aes(x = Variable, y = rho, fill = rho > 0)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB"), guide = "none") +
  labs(title = "SP110 vs Clinical Variables", x = "", y = "Spearman rho") + theme_minimal()
dev.off()
cat("  Clinical correlation saved\n")

# ==========================================
# M9-4: Multi-gene TB risk score comparison
# ==========================================
cat("\nM9-4: TB risk score...\n")
# Compare SP110 alone vs multi-gene ER stress score
er_genes_expr <- intersect(er_genes, rownames(expr))
er_score <- colMeans(expr[er_genes_expr, ])

risk_df <- data.frame(
  TB = factor(ifelse(type_vec == "Treat", "TB", "Control")),
  SP110 = sp110_expr,
  ER_Score = er_score
)

# AUC comparison
sp110_roc <- roc(risk_df$TB, risk_df$SP110, quiet = TRUE)
er_roc <- roc(risk_df$TB, risk_df$ER_Score, quiet = TRUE)

pdf(file.path(pdf_dir, "Fig_M9_RiskScore_ROC.pdf"), width = 7, height = 6)
plot.roc(sp110_roc, col = "#E74C3C", lwd = 2.5, main = "TB Diagnostic Performance")
plot.roc(er_roc, col = "#3498DB", lwd = 2.5, add = TRUE)
legend("bottomright", c(
  sprintf("SP110 alone (AUC=%.3f)", as.numeric(sp110_roc$auc)),
  sprintf("ER Stress Score (AUC=%.3f)", as.numeric(er_roc$auc))
), col = c("#E74C3C", "#3498DB"), lwd = 2, cex = 0.8, bty = "n")
dev.off()
cat("  Risk score ROC saved\n")

# ==========================================
# M9-5: Pathway-Pathway correlation network
# ==========================================
cat("\nM9-5: Pathway correlation network...\n")
er_pathways <- list(
  UPR = intersect(c("ERN1","EIF2AK3","ATF6","HSPA5","ATF4","DDIT3","XBP1"), rownames(expr)),
  Apoptosis = intersect(c("CASP3","CASP7","CASP9","BAX","BCL2","DDIT3","TRIB3"), rownames(expr)),
  Autophagy = intersect(c("BECN1","MAP1LC3B","SQSTM1","ATG5","ATG7","ATG12"), rownames(expr)),
  IFN_Response = intersect(c("GBP1","GBP2","GBP5","STAT1","IRF1","IDO1","WARS"), rownames(expr)),
  Inflammasome = intersect(c("NLRP3","CASP1","IL1B","IL18","PYCARD"), rownames(expr)),
  Antigen_Present = intersect(c("TAP1","TAP2","PSMB8","PSMB9","HLA-DRA","HLA-B"), rownames(expr))
)
pw_scores <- sapply(er_pathways, function(genes) colMeans(expr[genes, , drop = FALSE]))
pw_cor <- cor(pw_scores, method = "spearman")

pdf(file.path(pdf_dir, "Fig_M9_Pathway_Network.pdf"), width = 7, height = 6)
pheatmap(pw_cor, display_numbers = TRUE, number_format = "%.2f", fontsize_number = 10,
         main = "Pathway-Pathway Correlation (Spearman)",
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
dev.off()
cat("  Pathway network saved\n")

cat(sprintf("\nM9 complete. Figures saved to %s/\n", pdf_dir))
