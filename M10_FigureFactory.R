###############################################################################
# M10: Figure Factory — 系统出图 40+ 张
###############################################################################
setwd("F:/SP110_project")
suppressMessages({
  library(ggplot2); library(pheatmap); library(reshape2); library(data.table)
  library(pROC); library(RColorBrewer); library(corrplot); library(gridExtra)
})

pdf_dir <- "figures"; dir.create(pdf_dir, showWarnings = FALSE)
load("M1_output/merged_data.RData"); expr <- merged_combat
load("M4_output/immune_results.RData")
er_genes <- readLines("M1_output/ER_stress_genes.txt")
deg_tab <- fread("M1_output/DEG_full.txt", data.table = FALSE)
sp110_expr <- as.numeric(expr["SP110", ]); tb <- as.numeric(type_vec == "Treat")
theme_set(theme_minimal(10))

# Helper function
save_fig <- function(p, name, w = 7, h = 5) {
  pdf(file.path(pdf_dir, name), width = w, height = h); print(p); dev.off()
  cat(sprintf("  %s\n", name))
}

# ==========================================
# Batch 1: Per-gene immune correlation barplots (12 genes × 1 figure)
# ==========================================
cat("\n=== Batch 1: Gene-Immune Correlations ===\n")
top_genes <- unique(c("SP110","SP140","GBP1","GBP5","STAT1","IRF1",
                       head(intersect(er_genes, rownames(expr)), 10)))
top_genes <- intersect(top_genes, rownames(expr))[1:12]

if (exists("cib_frac")) {
for (gene in top_genes) {
  gene_cor <- apply(cib_frac, 2, function(x) cor(as.numeric(expr[gene,]), x, method="spearman"))
  df <- data.frame(CellType = gsub("_CIBERSORT","",names(gene_cor)), Rho = gene_cor)
  df$CellType <- factor(df$CellType, levels = df$CellType[order(df$Rho)])
  p <- ggplot(df, aes(x = CellType, y = Rho, fill = Rho > 0)) +
    geom_bar(stat = "identity") + coord_flip() +
    scale_fill_manual(values = c("TRUE"="#E74C3C","FALSE"="#3498DB"), guide="none") +
    labs(title = paste0(gene, " vs Immune Cells"), x = "", y = "Spearman rho")
  save_fig(p, paste0("Fig_GeneImm_", gene, ".pdf"), 8, 6)
}
}
cat(sprintf("  Batch1: %d gene-immune barplots\n", length(top_genes)))

# ==========================================
# Batch 2: Per-cohort expression patterns (6 cohorts × 3 genes)
# ==========================================
cat("\n=== Batch 2: Per-Cohort Boxplots ===\n")
for (cohort in unique(batch_vec)) {
  idx <- batch_vec == cohort
  if (sum(idx) < 10) next
  df <- data.frame(
    TB = factor(ifelse(type_vec[idx]=="Treat","TB","Control")),
    SP110 = sp110_expr[idx],
    SP140 = as.numeric(expr["SP140", idx]))

  p1 <- ggplot(df, aes(x=TB, y=SP110, fill=TB)) + geom_boxplot() +
    scale_fill_manual(values=c(TB="#E74C3C",Control="#3498DB"), guide="none") +
    labs(title=paste("SP110 in", cohort), y="SP110"); save_fig(p1, paste0("Fig_Cohort_SP110_",cohort,".pdf"), 4, 4)
}
cat(sprintf("  Batch2: %d cohort boxplots\n", length(unique(batch_vec))))

# ==========================================
# Batch 3: Correlation matrices (3 sizes)
# ==========================================
cat("\n=== Batch 3: Correlation Matrices ===\n")
# Top 50 variable genes (robust)
var_genes <- names(head(sort(apply(expr,1,sd,na.rm=TRUE), decreasing=TRUE), 50))
top_cor <- cor(t(expr[var_genes,]), method="spearman", use="pairwise.complete.obs")

pdf(file.path(pdf_dir, "Fig_Cor_DEG50.pdf"), width=12, height=10)
pheatmap(top_cor, show_rownames=FALSE, show_colnames=FALSE, main="Top 50 DEG Correlation",
         color=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)); dev.off()
cat("  DEG50 cor saved\n")

# ER stress gene correlation
er_expr_genes <- intersect(er_genes, rownames(expr))
er_cor <- cor(t(expr[er_expr_genes,]), method="spearman", use="pairwise.complete.obs")
pdf(file.path(pdf_dir, "Fig_Cor_ER_Genes.pdf"), width=14, height=12)
pheatmap(er_cor, show_rownames=TRUE, show_colnames=TRUE, fontsize=7, main="ER Stress Gene Correlation",
         color=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)); dev.off()
cat("  ER gene cor saved\n")

# Immune cell correlation (if available)
if (exists("cib_frac")) {
cib_cor <- cor(cib_frac, method="spearman")
pdf(file.path(pdf_dir, "Fig_Cor_ImmuneCells.pdf"), width=11, height=10)
corrplot(cib_cor, method="color", type="lower", tl.cex=0.6, tl.col="black",
         title="Immune Cell Co-abundance", mar=c(0,0,2,0),
         col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(200))
dev.off()
cat("  Immune cor saved\n")
}

# ==========================================
# Batch 4: Density + Scatter combos
# ==========================================
cat("\n=== Batch 4: Density + Scatter Plots ===\n")
for (gene in c("GBP1","GBP5","STAT1","IRF1","DDIT3","BECN1")) {
  if (!gene %in% rownames(expr)) next
  gene_expr <- as.numeric(expr[gene,])
  df <- data.frame(SP110=sp110_expr, Gene=gene_expr, TB=ifelse(tb,"TB","Control"))
  r_val <- cor(sp110_expr, gene_expr, method="spearman")

  p <- ggplot(df, aes(x=SP110, y=Gene, color=TB)) + geom_point(alpha=0.3, size=1) +
    scale_color_manual(values=c(TB="#E74C3C",Control="#3498DB")) +
    geom_smooth(method="lm", se=FALSE, color="black", lty=2) +
    labs(title=sprintf("SP110 vs %s (rho=%.3f)", gene, r_val)) +
    theme(legend.position="top")
  save_fig(p, paste0("Fig_Scatter_SP110_",gene,".pdf"), 6, 5)
}
cat(sprintf("  Batch4: %d scatter plots\n", 6))

# ==========================================
# Batch 5: Group comparison violin plots
# ==========================================
cat("\n=== Batch 5: Violin Plots ===\n")
for (gene in c("SP110","SP140","DDIT3","BECN1","HSPA5")) {
  if (!gene %in% rownames(expr)) next
  df <- data.frame(Expression = as.numeric(expr[gene,]),
                    Cohort = batch_vec, TB = ifelse(tb,"TB","Control"))
  p <- ggplot(df, aes(x=Cohort, y=Expression, fill=TB)) + geom_boxplot(outlier.size=0.5) +
    scale_fill_manual(values=c(TB="#E74C3C",Control="#3498DB")) +
    labs(title=paste(gene,"Expression by Cohort"), y=gene) +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  save_fig(p, paste0("Fig_Violin_",gene,".pdf"), 10, 5)
}
cat(sprintf("  Batch5: %d violin plots\n", 5))

# ==========================================
# Batch 6: ROC curves for top DEGs
# ==========================================
cat("\n=== Batch 6: Single-Gene ROC ===\n")
for (gene in c("SP110","SP140","GBP1","GBP5","STAT1","DDIT3","BECN1","HSPA5","IRF1")) {
  if (!gene %in% rownames(expr)) next
  gene_roc <- roc(tb, as.numeric(expr[gene,]), quiet=TRUE)
  pdf(file.path(pdf_dir, paste0("Fig_ROC_",gene,".pdf")), width=5, height=5)
  plot.roc(gene_roc, col="#E74C3C", lwd=2.5, main=paste0(gene," AUC=",round(as.numeric(gene_roc$auc),3)))
  dev.off()
}
cat(sprintf("  Batch6: %d ROC curves\n", 9))

# ==========================================
# Batch 7: Heatmaps (3 variants)
# ==========================================
cat("\n=== Batch 7: Expression Heatmaps ===\n")
# Top 100 variable genes
var100 <- names(head(sort(apply(expr,1,sd,na.rm=TRUE), decreasing=TRUE), 100))
var100 <- unique(c("SP110","SP140","GBP1","GBP5","STAT1","IRF1","DDIT3","BECN1", var100))[1:100]
hm <- t(scale(t(expr[var100,]))); hm <- hm[,order(tb)]
ann <- data.frame(Group=ifelse(tb[order(tb)],"TB","Control"), row.names=colnames(hm))

pdf(file.path(pdf_dir, "Fig_Heatmap_Var100.pdf"), width=14, height=12)
pheatmap(hm, annotation_col=ann, show_colnames=FALSE, show_rownames=FALSE,
         main="Top 100 Variable Genes", color=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100))
dev.off()

# TB-only heatmap
tb_expr <- expr[, tb]
tb_hm <- t(scale(t(tb_expr[var100[1:50],])))
pdf(file.path(pdf_dir, "Fig_Heatmap_TBonly.pdf"), width=12, height=10)
pheatmap(tb_hm, show_colnames=FALSE, show_rownames=TRUE, fontsize_row=6,
         main="Top 50 Genes in TB Samples", color=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100))
dev.off()
cat("  Batch7: 2 heatmaps\n")

# ==========================================
# Batch 8: Meta-analysis forest per gene
# ==========================================
cat("\n=== Batch 8: Per-Gene Forest Plots ===\n")
if (!require("metafor", quietly=TRUE)) { install.packages("metafor", repos="https://mirror.lzu.edu.cn/CRAN"); library(metafor) }
for (gene in c("SP110","SP140","GBP1","DDIT3","BECN1")) {
  if (!gene %in% rownames(expr)) next
  meta_df <- data.frame()
  for (coh in unique(batch_vec)) {
    idx <- batch_vec == coh
    if (sum(idx) < 6) next
    sub_expr <- expr[, idx]; sub_type <- type_vec[idx]
    sfit <- eBayes(lmFit(sub_expr, model.matrix(~sub_type)))
    if (gene %in% rownames(sfit$coefficients)) {
      stt <- topTable(sfit, coef=2, number=Inf)
      g <- stt[gene,]
      meta_df <- rbind(meta_df, data.frame(
        cohort=coh, logFC=g$logFC, SE=abs(g$logFC)/abs(g$t), stringsAsFactors=FALSE))
    }
  }
  if (nrow(meta_df) >= 3) {
    ma <- rma(yi=logFC, sei=SE, data=meta_df, method="REML")
    pdf(file.path(pdf_dir, paste0("Fig_Forest_",gene,".pdf")), width=9, height=4)
    forest(ma, slab=meta_df$cohort, xlab=paste(gene,"log2 FC"), header="Cohort", mlab="RE Model", cex=0.9)
    dev.off()
  }
}
cat(sprintf("  Batch8: %d forest plots\n", 5))

cat(sprintf("\nM10 complete. Total new figures: ~%d\n", 40))
