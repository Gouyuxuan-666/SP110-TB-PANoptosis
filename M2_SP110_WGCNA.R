###############################################################################
# SP110 M2 — WGCNA + ER Stress 共表达模块
# 在2号机运行: C:\R\R-4.6.1\bin\Rscript.exe M2_SP110_WGCNA.R
###############################################################################

setwd("F:/SP110_project")
suppressMessages({
  library(WGCNA); library(clusterProfiler); library(org.Hs.eg.db)
  library(enrichplot); library(ggplot2); library(pheatmap)
  library(GSVA); library(data.table)
})

allowWGCNAThreads(nThreads = 4)
enableWGCNAThreads()

# ---- Load merged data ----
cat("Loading M1 data...\n")
load("M1_output/merged_data.RData")
expr <- merged_combat

# ---- ER stress gene set ----
er_genes <- readLines("M1_output/ER_stress_genes.txt")
cat(sprintf("ER stress genes: %d\n", length(er_genes)))

# ---- WGCNA Input: top5000 MAD ∪ all DEGs ∪ ER stress genes ----
deg_tab <- fread("M1_output/DEG_full.txt", data.table = FALSE)
rownames(deg_tab) <- deg_tab[[1]]
deg_genes <- rownames(deg_tab)[abs(deg_tab$logFC) > 0.585 & deg_tab$adj.P.Val < 0.05]

mad_rank <- apply(expr, 1, mad)
top5000 <- names(sort(mad_rank, decreasing = TRUE)[1:5000])

wgcna_input <- union(union(top5000, deg_genes), intersect(er_genes, rownames(expr)))
wgcna_expr <- t(expr[wgcna_input, ])
cat(sprintf("WGCNA input: %d genes x %d samples\n", ncol(wgcna_expr), nrow(wgcna_expr)))

# ---- Soft threshold ----
cat("Selecting soft threshold...\n")
powers <- c(1:30)
sft <- pickSoftThreshold(wgcna_expr, powerVector = powers, verbose = 0,
                          networkType = "signed", corFnc = "bicor")
soft_power <- sft$powerEstimate
if (is.na(soft_power)) soft_power <- 4
cat(sprintf("Soft power: %d (R2=%.2f)\n", soft_power, sft$fitIndices[soft_power, "SFT.R.sq"]))

# ---- Network construction ----
cat("Constructing network...\n")
net <- blockwiseModules(wgcna_expr, power = soft_power,
                         TOMType = "signed", minModuleSize = 30,
                         deepSplit = 2, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = FALSE, verbose = 0, corType = "bicor",
                         maxBlockSize = 30000)

module_colors <- labels2colors(net$colors)
n_modules <- length(unique(module_colors))
cat(sprintf("Modules: %d\n", n_modules))

# ---- Module-trait correlation ----
tb_indicator <- as.numeric(type_vec == "Treat")
MEs <- net$MEs
module_trait_cor <- cor(MEs, tb_indicator, use = "pairwise.complete.obs")
module_trait_p <- apply(MEs, 2, function(x) cor.test(x, tb_indicator)$p.value)

# Find TB-associated module
cat("\nModule-TB correlations:\n")
for (i in seq_along(unique(module_colors))) {
  col <- unique(module_colors)[i]
  cat(sprintf("  %-15s r=%.3f P=%.2e\n", col, module_trait_cor[i], module_trait_p[i]))
}

top_module <- unique(module_colors)[which.max(abs(module_trait_cor))]
top_genes <- wgcna_input[module_colors == top_module]
cat(sprintf("\nTop module: %s (%d genes, r=%.3f)\n", top_module, length(top_genes), max(abs(module_trait_cor))))

# ---- SP110 in modules ----
sp110_module <- module_colors[which(wgcna_input == "SP110")[1]]
kME_all <- cor(t(expr[wgcna_input, ]), MEs, use = "pairwise.complete.obs")
sp110_kME <- kME_all[which(wgcna_input == "SP110"), which(unique(module_colors) == sp110_module)]
sp110_GS <- cor(expr["SP110", ], tb_indicator, use = "pairwise.complete.obs")
cat(sprintf("SP110: module=%s, kME=%.3f, GS=%.3f\n", sp110_module, sp110_kME, sp110_GS))

# ---- ER stress gene overlap with top module ----
er_in_top <- intersect(er_genes, top_genes)
cat(sprintf("ER stress genes in '%s': %d\n", top_module, length(er_in_top)))
cat("Top ER genes in module:\n")
er_kme <- kME_all[match(er_in_top, wgcna_input), which(unique(module_colors) == top_module)]
names(er_kme) <- er_in_top
print(sort(er_kme, decreasing = TRUE)[1:min(15, length(er_kme))])

# ---- Enrichment ----
cat("\nGO/KEGG enrichment...\n")
entrez_ids <- tryCatch({
  bitr(top_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
}, error = function(e) NULL)

if (!is.null(entrez_ids)) {
  ego <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                  ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.2)
  if (!is.null(ego) && nrow(ego) > 0) {
    ego <- simplify(ego, cutoff = 0.7)
    cat("Top GO terms:\n")
    print(head(ego[, c("Description","p.adjust")], 10))
  }
}

# ---- ER stress gene distribution across modules ----
cat("\nER stress gene distribution:\n")
for (g in intersect(er_genes, wgcna_input)) {
  mod <- module_colors[which(wgcna_input == g)]
  cat(sprintf("  %-12s -> %s\n", g, mod))
}

# ---- ssGSEA ----
cat("\nssGSEA pathway activity...\n")
er_pathway <- list(
  UPR = c("ERN1","EIF2AK3","ATF6","HSPA5","ATF4","DDIT3","XBP1"),
  ER_stress_apoptosis = c("DDIT3","CASP12","TRIB3","PPP1R15A","ATF4","BAX","BCL2"),
  Autophagy = c("BECN1","MAP1LC3B","SQSTM1","ATG5","ATG7","ATG12"),
  IFN_gamma = c("GBP1","GBP2","GBP5","STAT1","IRF1","IDO1","WARS")
)
ssgsea <- tryCatch(gsva(expr, er_pathway, method = "ssgsea", verbose = FALSE), error = function(e) NULL)
if (!is.null(ssgsea)) {
  for (pw in rownames(ssgsea)) {
    ct <- cor.test(as.numeric(expr["SP110", ]), ssgsea[pw, ], method = "spearman")
    cat(sprintf("  SP110 vs %-20s rho=%.3f P=%.2e\n", pw, ct$estimate, ct$p.value))
  }
} else {
  cat("ssGSEA failed (GSVA API issue). Using mean z-score instead...\n")
  for (pn in names(er_pathway)) {
    genes <- intersect(er_pathway[[pn]], rownames(expr))
    if (length(genes) >= 3) {
      score <- colMeans(expr[genes, , drop = FALSE])
      ct <- cor.test(as.numeric(expr["SP110", ]), score, method = "spearman")
      cat(sprintf("  SP110 vs %-20s rho=%.3f P=%.2e\n", pn, ct$estimate, ct$p.value))
    }
  }
}

# ---- Save ----
save(wgcna_input, module_colors, top_module, top_genes, MEs, kME_all, sp110_kME, sp110_GS, sp110_module, file = "M2_output/wgcna_results.RData")
writeLines(top_genes, "M2_output/top_module_genes.txt")
writeLines(er_in_top, "M2_output/er_genes_in_module.txt")

cat("\nM2 Complete. Output: M2_output/\n")
