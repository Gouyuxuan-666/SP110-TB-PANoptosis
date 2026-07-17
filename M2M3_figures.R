###############################################################################
# SP110 M1-M5 全图脚本 (v2 — 20张图)
###############################################################################
setwd("F:/SP110_project")
suppressMessages({
  library(WGCNA); library(ggplot2); library(pheatmap); library(pROC)
  library(data.table); library(RColorBrewer); library(ggrepel)
  library(gridExtra); library(reshape2); library(caret)
})

pdf_dir <- "figures"; dir.create(pdf_dir, showWarnings = FALSE)
theme_set(theme_minimal(base_size = 12))

# ---- Load data ----
load("M1_output/merged_data.RData"); expr <- merged_combat
load("M2_output/wgcna_results.RData")
load("M3_output/ml_results.RData")
load("M4_output/immune_results.RData")
er_genes <- readLines("M1_output/ER_stress_genes.txt")
deg_tab <- fread("M1_output/DEG_full.txt", data.table = FALSE)
deg_tab <- deg_tab[, -1]  # drop extra first column from fread
rownames(deg_tab) <- rownames(deg_tab)
sp110_expr <- as.numeric(expr["SP110", ])
tb_indicator <- as.numeric(type_vec == "Treat")

# Pre-compute SP110 correlations (needed by multiple sections)
sp110_cor <- apply(expr, 1, function(x) cor(sp110_expr, x, method = "spearman", use = "pairwise.complete.obs"))

# ==============================
# M1 Figures (4)
# ==============================

# M1-1: Volcano plot (SP110 + ER genes highlighted)
cat("M1: Volcano...\n")
deg_plot <- deg_tab[!is.na(deg_tab$P.Value), ]
deg_plot$sig <- "NS"
deg_plot$sig[abs(deg_plot$logFC) > 0.585 & deg_plot$adj.P.Val < 0.05] <- "DEG"
deg_plot$label <- ""
highlight <- c("SP110","SP140","DDIT3","BECN1","ERN1","EIF2AK3","HSPA5","GBP1","GBP5","STAT1")
for (g in intersect(highlight, rownames(deg_plot))) deg_plot[g, "label"] <- g

pdf(file.path(pdf_dir, "Fig_M1_Volcano.pdf"), width = 8, height = 7)
ggplot(deg_plot, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_color_manual(values = c(DEG = "#E74C3C", NS = "#BDC3C7")) +
  geom_vline(xintercept = c(-0.585, 0.585), lty = 2, alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), lty = 2, alpha = 0.5) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20) +
  labs(title = "SP110 & ER Stress Genes in TB Blood Transcriptome", x = "log2 FC", y = "-log10(P)") +
  theme(legend.position = "none")
dev.off()

# M1-2: SP110 per-cohort boxplot
cat("M1: Cohort boxplots...\n")
cohort_data <- data.frame(SP110 = sp110_expr, Cohort = batch_vec, TB = factor(type_vec, labels = c("Control","TB")))
pdf(file.path(pdf_dir, "Fig_M1_SP110_Cohort.pdf"), width = 12, height = 5)
ggplot(cohort_data, aes(x = Cohort, y = SP110, fill = TB)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = c(Control = "#3498DB", TB = "#E74C3C")) +
  labs(title = "SP110 Expression Across Cohorts", y = "SP110 Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# M1-3: Top DEG barplot
cat("M1: Top DEG barplot...\n")
top_degs <- head(deg_tab[abs(deg_tab$logFC) > 0.585 & deg_tab$adj.P.Val < 0.05, ], 20)
top_degs$Gene <- factor(rownames(top_degs), levels = rev(rownames(top_degs)))
pdf(file.path(pdf_dir, "Fig_M1_TopDEG.pdf"), width = 9, height = 6)
ggplot(top_degs, aes(x = Gene, y = logFC, fill = logFC > 0)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB"), guide = "none") +
  labs(title = "Top 20 DEGs (TB vs Control)", y = "log2 Fold Change")
dev.off()

# M1-4: SP110 vs ER stress gene correlation scatter
cat("M1: SP110-ER gene correlation...\n")
er_genes_expr <- intersect(er_genes, rownames(expr))
er_cors <- data.frame(Gene = er_genes_expr, Rho = sp110_cor[er_genes_expr])
er_cors <- er_cors[order(er_cors$Rho), ]
er_cors$Gene <- factor(er_cors$Gene, levels = er_cors$Gene)
pdf(file.path(pdf_dir, "Fig_M1_SP110_ER_Corr.pdf"), width = 10, height = 6)
ggplot(er_cors, aes(x = Gene, y = Rho, fill = Rho > 0)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB"), guide = "none") +
  labs(title = "SP110 vs ER Stress Gene Expression Correlation", y = "Spearman rho")
dev.off()

# ==============================
# M2 Figures (4)
# ==============================

# M2-1: Module-trait heatmap
cat("M2: Module-trait heatmap...\n")
keep_cols <- grep("^ME", colnames(MEs))
MEs_clean <- MEs[, keep_cols, drop = FALSE]

pdf(file.path(pdf_dir, "Fig_M2_ModuleTrait.pdf"), width = 8, height = 6)
MEs_ordered <- MEs_clean[, order(apply(MEs_clean, 2, function(x) cor(x, tb_indicator)))]
cor_mat <- cor(MEs_ordered, tb_indicator)
p_mat <- apply(MEs_ordered, 2, function(x) cor.test(x, tb_indicator)$p.value)
colnames(cor_mat) <- "TB_Status"
labeledHeatmap(Matrix = cor_mat, xLabels = "TB",
               yLabels = substring(colnames(MEs_ordered), 3),
               ySymbols = substring(colnames(MEs_ordered), 3),
               colors = blueWhiteRed(50), zlim = c(-1, 1),
               setStdMargins = FALSE,
               textMatrix = matrix(sprintf("%.2f\nP=%.1e", cor_mat, p_mat), ncol = 1),
               main = "Module-Trait Correlations")
dev.off()

# M2-2: ER stress heatmap
cat("M2: ER heatmap...\n")
er_expr <- expr[intersect(er_genes, rownames(expr)), ]
er_expr <- er_expr[rowSums(is.na(er_expr)) == 0, ]
er_scaled <- t(scale(t(er_expr)))
er_scaled <- er_scaled[, order(type_vec)]
ann_col <- data.frame(Group = ifelse(type_vec[order(type_vec)] == "Treat", "TB", "Control"), row.names = colnames(er_scaled))

pdf(file.path(pdf_dir, "Fig_M2_ER_Heatmap.pdf"), width = 14, height = 10)
pheatmap(er_scaled, annotation_col = ann_col, show_colnames = FALSE,
         cluster_cols = FALSE, fontsize_row = 7,
         main = "ER Stress Genes: TB vs Control",
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
dev.off()

# M2-3: kME vs GS scatter
cat("M2: kME vs GS...\n")
sp110_mod <- module_colors[which(wgcna_input == "SP110")]
sp110_mod_idx <- which(unique(module_colors) == sp110_mod)
if (length(sp110_mod_idx) == 0) sp110_mod_idx <- 1
gene_kME <- kME_all[, sp110_mod_idx]
gene_GS <- apply(expr[wgcna_input, ], 1, function(x) cor(x, tb_indicator, use = "pairwise.complete.obs"))

pdf(file.path(pdf_dir, "Fig_M2_kME_GS.pdf"), width = 7, height = 6)
plot(gene_kME, gene_GS, pch = 16, cex = 0.5, col = rgb(0.3, 0.3, 0.3, 0.5),
     xlab = paste0("kME in ", sp110_mod), ylab = "|cor with TB|",
     main = paste0("SP110: kME=", round(sp110_kME, 3), " GS=", round(sp110_GS, 3)))
sp110_idx <- which(wgcna_input == "SP110")
points(gene_kME[sp110_idx], gene_GS[sp110_idx], col = "#E74C3C", pch = 16, cex = 2.5)
text(gene_kME[sp110_idx], gene_GS[sp110_idx], "SP110", pos = 3, col = "#E74C3C", font = 2)
abline(h = 0, lty = 2); abline(v = 0, lty = 2)
dev.off()

# M2-4: GO enrichment dotplot (if enrichment results available)
cat("M2: GO dotplot...\n")
library(clusterProfiler); library(org.Hs.eg.db)
ego <- tryCatch({
  load("M2_output/enrichment.RData", verbose = FALSE)
  if (!exists("ego")) stop()
  ego
}, error = function(e) {
  tryCatch({
    top_genes <- readLines("M2_output/top_module_genes.txt")
    entrez <- clusterProfiler::bitr(top_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
    e <- enrichGO(entrez, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.2)
    if (!is.null(e) && nrow(e) > 0) simplify(e, 0.7) else NULL
  }, error = function(e2) NULL)
})

if (!is.null(ego) && nrow(ego) > 0) {
  ego_df <- head(ego@result[order(ego@result$p.adjust), ], 15)
  ego_df$Description <- factor(ego_df$Description, levels = rev(ego_df$Description))
  pdf(file.path(pdf_dir, "Fig_M2_GO_Dot.pdf"), width = 10, height = 6)
  ggplot(ego_df, aes(x = -log10(p.adjust), y = Description, size = Count, color = p.adjust)) +
    geom_point() + scale_color_gradient(low = "#E74C3C", high = "#3498DB") +
    labs(title = "Top GO Terms (SP110 Module)", x = "-log10(FDR)")
  dev.off()
}

# ==============================
# M3 Figures (4)
# ==============================

# M3-1: ROC curve
cat("M3: ROC...\n")
pdf(file.path(pdf_dir, "Fig_M3_ROC.pdf"), width = 6, height = 6)
plot.roc(test_roc, col = "#E74C3C", lwd = 3,
         main = sprintf("SP110-based Model\nAUC=%.3f", test_auc),
         print.auc = TRUE, auc.polygon = TRUE, auc.polygon.col = rgb(0.9,0.3,0.2,0.2))
dev.off()

# M3-2: Feature importance
cat("M3: Feature importance...\n")
imp_df <- as.data.frame(imp)
imp_df$Gene <- rownames(imp_df)
imp_df <- imp_df[order(imp_df$Overall, decreasing = TRUE)[1:20], ]
imp_df$Gene <- factor(imp_df$Gene, levels = rev(imp_df$Gene))

pdf(file.path(pdf_dir, "Fig_M3_FeatureImportance.pdf"), width = 9, height = 6)
ggplot(imp_df, aes(x = Gene, y = Overall, fill = Gene == "SP110")) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB"), guide = "none") +
  coord_flip() + labs(title = "Top 20 ML Features", y = "Importance (scaled)", x = "")
dev.off()

# M3-3: Confusion matrix heatmap
cat("M3: Confusion matrix...\n")
ml_features <- intersect(unique(c(names(head(sort(abs(sp110_cor), decreasing=TRUE), 500)), er_genes)), rownames(expr))
ml_expr <- t(expr[ml_features, ])
set.seed(20240717)
train_idx <- createDataPartition(factor(ifelse(type_vec == "Treat", "TB", "Control")), p = 0.7, list = FALSE)
test_x <- as.data.frame(ml_expr[-train_idx, sel_genes, drop = FALSE])
test_y <- factor(ifelse(type_vec[-train_idx] == "Treat", "TB", "Control"))
pred_class <- ifelse(predict(rf_model, newdata = test_x, type = "prob")[, "TB"] > 0.5, "TB", "Control")
cm <- table(Actual = test_y, Predicted = pred_class)
cm_df <- as.data.frame(as.table(cm))
pdf(file.path(pdf_dir, "Fig_M3_Confusion.pdf"), width = 5, height = 4)
ggplot(cm_df, aes(x = Predicted, y = Actual, fill = Freq)) +
  geom_tile(color = "white") + geom_text(aes(label = Freq), size = 8) +
  scale_fill_gradient(low = "#E8F5E9", high = "#2E7D32") +
  labs(title = "Confusion Matrix (Internal Test)") + theme_minimal()
dev.off()

# M3-4: Calibration curve
cat("M3: Calibration...\n")
pred_prob <- predict(rf_model, newdata = test_x, type = "prob")[, "TB"]
calib <- data.frame(Predicted = pred_prob, Actual = as.numeric(test_y == "TB"))
calib <- calib[order(calib$Predicted), ]
calib$Bin <- cut(calib$Predicted, breaks = 10, labels = FALSE)
calib_summ <- aggregate(cbind(Predicted, Actual) ~ Bin, calib, mean)

pdf(file.path(pdf_dir, "Fig_M3_Calibration.pdf"), width = 5, height = 5)
plot(calib_summ$Predicted, calib_summ$Actual, type = "b", col = "#E74C3C", lwd = 2,
     xlim = 0:1, ylim = 0:1, xlab = "Predicted Probability", ylab = "Observed Proportion",
     main = "Calibration Plot")
abline(0, 1, lty = 2, col = "gray")
dev.off()

# ==============================
# M4 Figures (4)
# ==============================

# M4-1: Immune correlation barplot
cat("M4: Immune correlation...\n")
imm_df <- data.frame(CellType = names(sp110_cor_imm), Rho = sp110_cor_imm)
imm_df$CellType <- gsub("_CIBERSORT", "", imm_df$CellType)
imm_df$CellType <- factor(imm_df$CellType, levels = imm_df$CellType[order(imm_df$Rho)])

pdf(file.path(pdf_dir, "Fig_M4_Immune_Corr.pdf"), width = 10, height = 6)
ggplot(imm_df, aes(x = CellType, y = Rho, fill = Rho > 0)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB"), guide = "none") +
  labs(title = "SP110 vs 22 Immune Cell Types (Spearman rho)", x = "", y = "rho")
dev.off()

# M4-2: SP110 high vs low boxplot (top 6 cell types)
cat("M4: SP110 high vs low boxplots...\n")
top6 <- names(sp110_cor_imm)[c(1,2,3,20,21,22)]  # top 3 positive, bottom 3 negative
cib_long <- melt(cib_frac[, top6])
colnames(cib_long) <- c("CellType", "Fraction")
cib_long$SP110 <- ifelse(sp110_expr > median(sp110_expr), "High", "Low")

pdf(file.path(pdf_dir, "Fig_M4_SP110_HighLow.pdf"), width = 14, height = 8)
ggplot(cib_long, aes(x = SP110, y = Fraction, fill = SP110)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~CellType, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c(High = "#E74C3C", Low = "#3498DB")) +
  labs(title = "Immune Cell Fractions: SP110 High vs Low") + theme(legend.position = "none")
dev.off()

# M4-3: Immune cell inter-correlation heatmap
cat("M4: Immune correlation heatmap...\n")
cib_cor <- cor(cib_frac, method = "spearman")
pdf(file.path(pdf_dir, "Fig_M4_Immune_Network.pdf"), width = 10, height = 9)
pheatmap(cib_cor, show_rownames = TRUE, show_colnames = TRUE, fontsize = 8,
         main = "Immune Cell Co-abundance (Spearman)", color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
dev.off()

# M4-4: SP110 vs immune checkpoints
cat("M4: Immune checkpoints...\n")
icb_genes <- c("CD274","PDCD1","CTLA4","HAVCR2","LAG3","TIGIT","IDO1")
icb_found <- intersect(icb_genes, rownames(expr))
icb_cors <- data.frame(Gene = icb_found, Rho = sapply(icb_found, function(g) cor(sp110_expr, as.numeric(expr[g,]), method="spearman")))
icb_cors$Gene <- factor(icb_cors$Gene, levels = icb_cors$Gene[order(icb_cors$Rho)])

pdf(file.path(pdf_dir, "Fig_M4_Checkpoints.pdf"), width = 6, height = 4)
ggplot(icb_cors, aes(x = Gene, y = Rho, fill = Rho > 0)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB"), guide = "none") +
  labs(title = "SP110 vs Immune Checkpoints", x = "", y = "Spearman rho")
dev.off()

# ==============================
# M5 Figures (4)
# ==============================

cat("M5: Network figures...\n")

# M5-1: SP110 top co-expressed genes
cor_vals <- sp110_cor
top20_genes <- names(head(sort(cor_vals, decreasing = TRUE), 20))
top20_vals <- cor_vals[top20_genes]
top20_df <- data.frame(Gene = factor(top20_genes, levels = rev(top20_genes)), Rho = top20_vals)
pdf(file.path(pdf_dir, "Fig_M5_Top20Coexp.pdf"), width = 8, height = 6)
ggplot(top20_df, aes(x = Gene, y = Rho, fill = Rho > 0)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB"), guide = "none") +
  labs(title = "Top 20 SP110-Coexpressed Genes", x = "", y = "Spearman rho")
dev.off()

# M5-2: Pathway correlation summary
cat("M5: Pathway correlations...\n")
pathways <- list(
  ER_Sensors = c("ERN1","EIF2AK3","ATF6","HSPA5"),
  Apoptosis = c("DDIT3","CASP12","CASP3","BAX","BCL2"),
  Autophagy = c("BECN1","MAP1LC3B","SQSTM1","ATG5","ATG7"),
  Inflammasome = c("NLRP3","CASP1","IL1B","IL18"),
  IFN_Response = c("GBP1","GBP2","GBP5","STAT1","IRF1","IDO1"),
  Antigen_Present = c("TAP1","TAP2","PSMB8","PSMB9","HLA-DRA")
)
pw_df <- data.frame()
for (pn in names(pathways)) {
  genes <- intersect(pathways[[pn]], rownames(expr))
  cors <- sp110_cor[genes]
  pw_df <- rbind(pw_df, data.frame(Pathway = pn, MeanAbsRho = mean(abs(cors)), nGenes = length(genes)))
}
pw_df$Pathway <- factor(pw_df$Pathway, levels = pw_df$Pathway[order(pw_df$MeanAbsRho)])

pdf(file.path(pdf_dir, "Fig_M5_Pathway_Corr.pdf"), width = 8, height = 4)
ggplot(pw_df, aes(x = Pathway, y = MeanAbsRho, fill = MeanAbsRho)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = nGenes), vjust = -0.3, size = 3) +
  scale_fill_gradient(low = "#3498DB", high = "#E74C3C") +
  labs(title = "SP110 vs Pathway Gene Correlation (mean |rho|)", y = "Mean |rho|", subtitle = "Numbers = genes in pathway") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

# M5-3: SP110 vs GBP1 expression scatter
if ("GBP1" %in% rownames(expr)) {
  cat("M5: SP110 vs GBP1 scatter...\n")
  gbp1_expr <- as.numeric(expr["GBP1", ])
  df_gb <- data.frame(SP110 = sp110_expr, GBP1 = gbp1_expr, TB = factor(type_vec))
  r_val <- cor(sp110_expr, gbp1_expr, method = "spearman")

  pdf(file.path(pdf_dir, "Fig_M5_SP110vsGBP1.pdf"), width = 6, height = 6)
  ggplot(df_gb, aes(x = SP110, y = GBP1, color = TB)) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_manual(values = c(Control = "#3498DB", TB = "#E74C3C")) +
    geom_smooth(method = "lm", se = FALSE, color = "black", lty = 2) +
    labs(title = sprintf("SP110 vs GBP1 (rho = %.3f)", r_val), x = "SP110 Expression", y = "GBP1 Expression")
  dev.off()
}

# M5-4: SP110 network degree distribution
cat("M5: Network degree...\n")
abs_cors <- abs(sp110_cor)
degree_df <- data.frame(Threshold = seq(0.3, 0.9, 0.05),
                         nGenes = sapply(seq(0.3, 0.9, 0.05), function(t) sum(abs_cors > t, na.rm=TRUE)))
pdf(file.path(pdf_dir, "Fig_M5_NetworkDegree.pdf"), width = 5, height = 4)
ggplot(degree_df, aes(x = Threshold, y = nGenes)) +
  geom_line(color = "#E74C3C", lwd = 1.5) + geom_point(size = 3, color = "#E74C3C") +
  labs(title = "SP110 Co-expression Network Size", x = "|rho| Threshold", y = "Connected Genes")
dev.off()

# ==============================
# M6: LINKET — SP110/ER genes × immune cells (alpha style)
# ==============================
cat("M6: LINKET plot (alpha style)...\n")
linket_ok <- requireNamespace("linkET", quietly = TRUE) && requireNamespace("dplyr", quietly = TRUE)
if (linket_ok) { library(linkET); library(dplyr) }

# Select genes: SP110 + top ER stress + hallmark genes
link_genes <- unique(c("SP110","SP140","ERN1","EIF2AK3","DDIT3","BECN1","HSPA5",
                        "GBP1","GBP5","STAT1","IRF1","IDO1","CASP3","BAX","BCL2",
                        head(intersect(er_genes, rownames(expr)), 10)))
link_genes <- intersect(link_genes, rownames(expr))

# Use only TB samples for correlation (like alpha script)
tb_samples <- type_vec == "Treat"
gene_data <- t(expr[link_genes, tb_samples & cib_pass, drop = FALSE])
imm_data <- cib_frac[tb_samples & cib_pass, ]
# Remove columns with zero variance
imm_data <- imm_data[, apply(imm_data, 2, sd) > 0, drop = FALSE]

cat(sprintf("  Genes: %d, Immune cells: %d, TB samples: %d\n",
            ncol(gene_data), ncol(imm_data), nrow(gene_data)))

# Gene-immune Spearman correlation
geneCor <- data.frame()
for (cell in colnames(imm_data)) {
  for (gene in colnames(gene_data)) {
    x <- as.numeric(gene_data[, gene])
    y <- as.numeric(imm_data[, cell])
    ct <- cor.test(x, y, method = "spearman")
    geneCor <- rbind(geneCor, data.frame(
      spec = gene, env = cell,  # Keep original CIBERSORT names for linkET matching
      r = as.numeric(ct$estimate), p = as.numeric(ct$p.value),
      stringsAsFactors = FALSE))
  }
}

geneCor$pd <- ifelse(geneCor$p < 0.05,
  ifelse(geneCor$r > 0, "Positive", "Negative"), "Not")
geneCor$r_abs <- abs(geneCor$r)
geneCor <- geneCor %>% mutate(
  rd = cut(r_abs, breaks = c(-Inf, 0.2, 0.4, 0.6, Inf),
           labels = c("< 0.2", "0.2 - 0.4", "0.4 - 0.6", ">= 0.6")))

# Immune cell inter-correlation
imm_cor <- cor(imm_data, method = "spearman")

# LINKET plot (or fallback heatmap)
pdf(file.path(pdf_dir, "Fig_M6_LINKET.pdf"), width = 12, height = 9)
if (linket_ok) {
  qcorPlot <- qcorrplot(correlate(imm_data, method = "spearman"),
                         type = "lower", diag = FALSE) +
    geom_square() +
    geom_couple(aes(colour = pd, size = rd),
                data = geneCor, curvature = nice_curvature()) +
    scale_fill_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
    scale_size_manual(values = c(0.8, 1.5, 2.5, 3.5)) +
    scale_colour_manual(values = c("#E6550DFF", "#CCCCCC99", "#6699FFFF")) +
    guides(size = guide_legend(title = "|rho|", override.aes = list(colour = "grey35"), order = 2),
           colour = guide_legend(title = "p < 0.05", override.aes = list(size = 3), order = 1),
           fill = guide_colorbar(title = "Immune cell\ncorrelation", order = 3)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7)) +
    labs(title = "SP110 & ER Stress Genes ↔ Immune Microenvironment")
  print(qcorPlot)
} else {
  cat("  linkET not installed. Using pheatmap fallback.\n")
  gene_imm_cor <- cor(t(expr[link_genes, tb_samples & cib_pass]),
                       cib_frac[tb_samples & cib_pass, ], method = "spearman")
  pheatmap(gene_imm_cor, show_rownames = TRUE, show_colnames = TRUE, fontsize = 7,
           main = "SP110 & ER Stress Genes vs Immune Cells (Spearman rho)",
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
}
dev.off()

write.table(geneCor, file.path(pdf_dir, "Fig_M6_LINKET_data.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\nAll 21 figures saved to %s/\n", pdf_dir))
cat("M1: 4 | M2: 4 | M3: 4 | M4: 4 | M5: 4\n")
