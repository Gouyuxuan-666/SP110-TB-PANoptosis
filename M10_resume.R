###############################################################################
# M10 Resume — 从 Batch 7 续跑
###############################################################################
setwd("F:/SP110_project")
suppressMessages({
  library(ggplot2); library(pheatmap); library(data.table)
  library(pROC); library(RColorBrewer); library(limma); library(metafor)
})

pdf_dir <- "figures"
load("M1_output/merged_data.RData"); expr <- merged_combat
sp110_expr <- as.numeric(expr["SP110", ])
load("M4_output/immune_results.RData")
load("M2_output/wgcna_results.RData")
tb <- as.numeric(type_vec == "Treat")
theme_set(theme_minimal(10))

save_fig <- function(p, name, w = 7, h = 5) {
  pdf(file.path(pdf_dir, name), width = w, height = h); print(p); dev.off()
  cat(sprintf("  %s\n", name))
}

cat("=== Batch 7: Heatmaps (fixed) ===\n")
var100 <- names(head(sort(apply(expr,1,sd,na.rm=TRUE), decreasing=TRUE), 100))
var100 <- unique(c("SP110","SP140","GBP1","GBP5","STAT1","IRF1","DDIT3","BECN1", var100))[1:100]
hm <- t(scale(t(expr[var100,])))
hm[!is.finite(hm)] <- 0
hm <- hm[, order(tb)]
ann <- data.frame(Group=ifelse(tb[order(tb)],"TB","Control"), row.names=colnames(hm))

pdf(file.path(pdf_dir, "Fig_Heatmap_Var100.pdf"), width=14, height=12)
pheatmap(hm, annotation_col=ann, show_colnames=FALSE, show_rownames=FALSE,
         main="Top 100 Variable Genes", color=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100))
dev.off(); cat("  Var100 heatmap saved\n")

cat("  TB-only heatmap skipped (scaling issue)\n")

cat("\n=== Batch 8: Forest Plots ===\n")
for (gene in c("SP110","SP140","GBP1","DDIT3","BECN1")) {
  if (!gene %in% rownames(expr)) next
  meta_df <- data.frame()
  for (coh in unique(batch_vec)) {
    idx <- batch_vec == coh
    if (sum(idx) < 6 || sum(type_vec[idx]=="Treat") < 2 || sum(type_vec[idx]=="Control") < 2) next
    sub_expr <- expr[, idx]; sub_type <- type_vec[idx]
    sfit <- tryCatch(eBayes(lmFit(sub_expr, model.matrix(~sub_type))), error=function(e) NULL)
    if (is.null(sfit) || !gene %in% rownames(sfit$coefficients)) next
    stt <- topTable(sfit, coef=2, number=Inf)
    g <- stt[gene,]; if (is.na(g$t) || g$t == 0) next
    meta_df <- rbind(meta_df, data.frame(cohort=coh, logFC=g$logFC, SE=abs(g$logFC)/abs(g$t)))
  }
  if (nrow(meta_df) >= 3) {
    ma <- rma(yi=logFC, sei=SE, data=meta_df, method="REML")
    pdf(file.path(pdf_dir, paste0("Fig_Forest_",gene,".pdf")), width=9, height=4)
    forest(ma, slab=meta_df$cohort, xlab=paste(gene,"log2 FC"), header="Cohort", mlab="RE Model", cex=0.9)
    dev.off(); cat(sprintf("  Forest: %s (I2=%.1f%%)\n", gene, ma$I2))
  } else { cat(sprintf("  Forest: %s skipped (n=%d)\n", gene, nrow(meta_df))) }
}

cat("\n=== Batch 9: TB vs Control Barplot per ER gene ===\n")
er_genes <- readLines("M1_output/ER_stress_genes.txt")
er_expr <- intersect(er_genes, rownames(expr))
for (gene in er_expr[1:15]) {
  df <- data.frame(TB=ifelse(tb,"TB","Control"), Expression=as.numeric(expr[gene,]))
  p <- ggplot(df, aes(x=TB, y=Expression, fill=TB)) + geom_boxplot(outlier.size=0.5) +
    scale_fill_manual(values=c(TB="#E74C3C",Control="#3498DB"), guide="none") +
    labs(title=paste(gene,"in TB vs Control"), y=gene)
  save_fig(p, paste0("Fig_ER_",gene,".pdf"), 4, 4)
}

cat("\n=== Batch 10: SP110 high vs low pathway heatmap ===\n")
sp110_hi <- sp110_expr > median(sp110_expr)
er_paths <- list(
  UPR=intersect(c("ERN1","EIF2AK3","ATF6","HSPA5","ATF4","DDIT3","XBP1"),rownames(expr)),
  Apop=intersect(c("CASP3","CASP7","CASP9","BAX","BCL2","DDIT3","TRIB3"),rownames(expr)),
  Auto=intersect(c("BECN1","MAP1LC3B","SQSTM1","ATG5","ATG7"),rownames(expr)),
  Inflam=intersect(c("NLRP3","CASP1","IL1B","IL18"),rownames(expr)),
  IFN=intersect(c("GBP1","GBP2","GBP5","STAT1","IRF1","IDO1"),rownames(expr)))
pw_hi <- sapply(er_paths, function(g) colMeans(expr[g, sp110_hi, drop=FALSE]))
pw_lo <- sapply(er_paths, function(g) colMeans(expr[g, !sp110_hi, drop=FALSE]))
pw_diff <- data.frame(Pathway=names(er_paths), Delta=colMeans(pw_hi)-colMeans(pw_lo))
pw_diff$Pathway <- factor(pw_diff$Pathway, levels=pw_diff$Pathway[order(pw_diff$Delta)])

pdf(file.path(pdf_dir, "Fig_SP110_HiLo_Pathway.pdf"), width=7, height=4)
ggplot(pw_diff, aes(x=Pathway, y=Delta, fill=Delta>0)) + geom_bar(stat="identity") +
  scale_fill_manual(values=c("TRUE"="#E74C3C","FALSE"="#3498DB"), guide="none") +
  labs(title="Pathway Activity: SP110 High vs Low", y="Delta (High-Low)") + coord_flip()
dev.off(); cat("  Pathway hi/lo saved\n")

cat("\n=== Done ===\n")
