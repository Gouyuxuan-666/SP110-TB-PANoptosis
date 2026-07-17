###############################################################################
# SP110 M5 — 调控网络 + PPI + 与GBP1网络交叉比较
###############################################################################
setwd("F:/SP110_project")
suppressMessages({
  library(data.table); library(ggplot2)
})

load("M1_output/merged_data.RData")
er_genes <- readLines("M1_output/ER_stress_genes.txt")

# ---- 1. SP110 co-expression network ----
cat("SP110 co-expression network...\n")
sp110_expr <- as.numeric(expr["SP110", ])
sp110_cor <- apply(expr, 1, function(x) cor(sp110_expr, x, method = "spearman", use = "pairwise.complete.obs"))
sp110_cor <- sp110_cor[!is.na(sp110_cor)]

# Top 100 co-expressed genes
top100_cor <- names(sort(abs(sp110_cor), decreasing = TRUE)[1:100])
cor_vals <- sp110_cor[top100_cor]

# ER stress overlap
er_overlap <- intersect(top100_cor, er_genes)
cat(sprintf("ER stress genes in top100 SP110-coexpressed: %d\n", length(er_overlap)))
if (length(er_overlap) > 0) cat(paste(er_overlap, collapse = ", "), "\n")

# Print top genes
cat("\nTop 20 SP110-coexpressed genes:\n")
top20 <- head(sort(cor_vals, decreasing = TRUE), 20)
for (i in seq_along(top20)) {
  cat(sprintf("  %-12s rho=%.3f\n", names(top20)[i], top20[i]))
}

# ---- 2. Pathway correlation ----
cat("\nSP110 vs key pathway genes:\n")
pathways <- list(
  ER_stress_sensors = c("ERN1","EIF2AK3","ATF6","HSPA5"),
  Apoptosis = c("DDIT3","CASP12","CASP3","CASP7","CASP9","BAX","BCL2"),
  Autophagy = c("BECN1","MAP1LC3B","SQSTM1","ATG5","ATG7"),
  Inflammasome = c("NLRP3","CASP1","IL1B","IL18"),
  IFN_response = c("GBP1","GBP2","GBP5","STAT1","IRF1","IDO1"),
  Antigen_present = c("TAP1","TAP2","PSMB8","PSMB9","HLA-DRA")
)
for (pn in names(pathways)) {
  genes <- intersect(pathways[[pn]], rownames(expr))
  cors <- sp110_cor[genes]
  cat(sprintf("  %-20s mean_rho=%.3f (n=%d)\n", pn, mean(abs(cors), na.rm = TRUE), length(genes)))
  if (length(genes) > 0) {
    top_g <- names(sort(abs(cors), decreasing = TRUE))[1:min(3, length(genes))]
    for (g in top_g) cat(sprintf("    %s rho=%.3f\n", g, cors[g]))
  }
}

# ---- 3. Comparison: SP110 vs GBP1 co-expression patterns ----
cat("\n=== SP110 vs GBP1 network comparison ===\n")
if ("GBP1" %in% rownames(expr)) {
  gbp1_expr <- as.numeric(expr["GBP1", ])
  gbp1_cor <- apply(expr, 1, function(x) cor(gbp1_expr, x, method = "spearman", use = "pairwise.complete.obs"))
  gbp1_cor <- gbp1_cor[!is.na(gbp1_cor)]

  top100_gbp1 <- names(sort(abs(gbp1_cor), decreasing = TRUE)[1:100])
  shared <- intersect(top100_cor, top100_gbp1)

  cat(sprintf("Shared top100 genes (SP110 ∩ GBP1): %d\n", length(shared)))
  cat("Shared genes: "); cat(head(shared, 20), sep = ", "); cat("\n")

  # Correlation of correlations
  common <- intersect(names(sp110_cor), names(gbp1_cor))
  r_cor <- cor(sp110_cor[common], gbp1_cor[common], method = "spearman")
  cat(sprintf("SP110 vs GBP1 co-expression profile correlation: rho=%.3f (n=%d genes)\n", r_cor, length(common)))
}

# ---- 4. SP110 vs GBP1 co-expression correlation ----
if ("GBP1" %in% rownames(expr)) {
  ct <- cor.test(sp110_expr, gbp1_expr, method = "spearman")
  cat(sprintf("\nSP110 vs GBP1 expression: rho=%.3f P=%.2e\n", ct$estimate, ct$p.value))
}

# ---- Save ----
save(top100_cor, cor_vals, er_overlap, file = "M5_output/network_results.RData")
writeLines(top100_cor, "M5_output/SP110_top100_coexpressed.txt")

cat("\nM5 Complete.\n")
