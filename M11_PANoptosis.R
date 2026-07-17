###############################################################################
# M11: PANoptosis + 铁死亡 + 细胞衰老 — 创新分析
###############################################################################
setwd("F:/SP110_project")
suppressMessages({
  library(ggplot2); library(pheatmap); library(reshape2); library(data.table)
  library(pROC); library(RColorBrewer); library(limma)
})

pdf_dir <- "figures"; dir.create(pdf_dir, showWarnings = FALSE)
load("M1_output/merged_data.RData"); expr <- merged_combat
sp110_expr <- as.numeric(expr["SP110", ]); tb <- type_vec == "Treat"
theme_set(theme_minimal(10))

s <- function(p, n, w=7, h=5) { pdf(file.path(pdf_dir,n),w,h); print(p); dev.off() }

# ==========================================
# 1. PANoptosis (pyroptosis + apoptosis + necroptosis)
# ==========================================
cat("=== PANoptosis ===\n")
pan_genes <- list(
  Pyroptosis = c("NLRP3","CASP1","GSDMD","IL1B","IL18","AIM2","PYCARD","NLRC4"),
  Apoptosis = c("CASP3","CASP7","CASP8","CASP9","BAX","BAK1","BCL2","MCL1","BID","CYCS","APAF1","PARP1"),
  Necroptosis = c("RIPK1","RIPK3","MLKL","ZBP1","FADD","TRADD","TNFR1"),
  Autophagy = c("BECN1","MAP1LC3B","SQSTM1","ATG5","ATG7","ATG12","ATG16L1","ULK1"),
  Ferroptosis = c("GPX4","SLC7A11","ACSL4","TFRC","FTH1","FTL","HMOX1","PTGS2")
)
pan_all <- unique(unlist(pan_genes))
pan_expr <- intersect(pan_all, rownames(expr))

# PANoptosis score per sample
pan_scores <- sapply(pan_genes, function(g) {
  genes <- intersect(g, pan_expr)
  if (length(genes) >= 3) colMeans(expr[genes, , drop=FALSE]) else rep(NA, ncol(expr))
})
pan_scores <- pan_scores[, colSums(is.na(pan_scores)) == 0]

# SP110 vs PANoptosis pathways
pan_cor <- data.frame()
for (pw in colnames(pan_scores)) {
  ct <- cor.test(sp110_expr, pan_scores[,pw], method="spearman")
  pan_cor <- rbind(pan_cor, data.frame(Pathway=pw, Rho=ct$estimate, P=ct$p.value))
}
pan_cor <- pan_cor[order(pan_cor$Rho, decreasing=TRUE),]
write.csv(pan_cor, file.path(pdf_dir, "Table_PANoptosis_cor.csv"), row.names = FALSE)
print(pan_cor)

# PANoptosis heatmap
pan_hm <- t(scale(t(expr[pan_expr, order(tb)])))
pan_hm[!is.finite(pan_hm)] <- 0
ann_col <- data.frame(Group=ifelse(tb[order(tb)],"TB","Control"), row.names=colnames(pan_hm))
ann_vec <- rep("Other", length(pan_expr)); names(ann_vec) <- pan_expr
for (pw in names(pan_genes)) { hits <- intersect(pan_genes[[pw]], pan_expr); if(length(hits)>0) ann_vec[hits] <- pw }
ann_row <- data.frame(Pathway=ann_vec, row.names=pan_expr)

pdf(file.path(pdf_dir, "Fig_M11_PANoptosis_Heatmap.pdf"), width=14, height=10)
pheatmap(pan_hm, annotation_col=ann_col, annotation_row=ann_row,
         show_colnames=FALSE, fontsize_row=6, main="PANoptosis Genes: TB vs Control",
         color=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100))
dev.off(); cat("  PANoptosis heatmap saved\n")

# PANoptosis pathway barplot
pdf(file.path(pdf_dir, "Fig_M11_PANoptosis_Bar.pdf"), width=6, height=4)
ggplot(pan_cor, aes(x=reorder(Pathway,Rho), y=Rho, fill=Rho>0)) +
  geom_bar(stat="identity") + coord_flip() +
  scale_fill_manual(values=c("TRUE"="#E74C3C","FALSE"="#3498DB"), guide="none") +
  labs(title="SP110 vs Cell Death Pathways", x="", y="Spearman rho")
dev.off(); cat("  PANoptosis bar saved\n")

# ==========================================
# 2. SP110 high vs low PANoptosis comparison
# ==========================================
sp110_hi <- sp110_expr > median(sp110_expr)
pan_diff <- data.frame()
for (pw in colnames(pan_scores)) {
  hi <- mean(pan_scores[sp110_hi, pw]); lo <- mean(pan_scores[!sp110_hi, pw])
  p <- t.test(pan_scores[sp110_hi, pw], pan_scores[!sp110_hi, pw])$p.value
  pan_diff <- rbind(pan_diff, data.frame(Pathway=pw, Hi=hi, Lo=lo, Delta=hi-lo, P=p))
}
pan_diff$Pathway <- factor(pan_diff$Pathway, levels=pan_diff$Pathway[order(pan_diff$Delta)])

pdf(file.path(pdf_dir, "Fig_M11_SP110_HiLo_Death.pdf"), width=6, height=4)
ggplot(pan_diff, aes(x=Pathway, y=Delta, fill=Delta>0)) +
  geom_bar(stat="identity") + coord_flip() +
  scale_fill_manual(values=c("TRUE"="#E74C3C","FALSE"="#3498DB"), guide="none") +
  labs(title="Cell Death Pathway Activity: SP110 High vs Low", y="Delta")
dev.off(); cat("  Hi/Lo death saved\n")

# ==========================================
# 3. Cellular senescence genes
# ==========================================
cat("\n=== Cellular Senescence ===\n")
sen_genes <- c("CDKN1A","CDKN2A","TP53","RB1","CCND1","CDK4","CDK6","MYC",
               "IL6","IL8","CXCL1","CXCL8","CCL2","MMP1","MMP3","MMP9",
               "SERPINE1","IGFBP3","IGFBP5","TNF","TGFB1","LMNB1","HMGB1")
sen_expr <- intersect(sen_genes, rownames(expr))

sen_cor <- sapply(sen_expr, function(g) cor(sp110_expr, as.numeric(expr[g,]), method="spearman"))
sen_df <- data.frame(Gene=names(sen_cor), Rho=sen_cor)
sen_df <- sen_df[order(sen_df$Rho),]
sen_df$Gene <- factor(sen_df$Gene, levels=sen_df$Gene)

pdf(file.path(pdf_dir, "Fig_M11_Senescence.pdf"), width=8, height=5)
ggplot(sen_df, aes(x=Gene, y=Rho, fill=Rho>0)) +
  geom_bar(stat="identity") + coord_flip() +
  scale_fill_manual(values=c("TRUE"="#E74C3C","FALSE"="#3498DB"), guide="none") +
  labs(title="SP110 vs Cellular Senescence Genes", x="", y="Spearman rho")
dev.off(); cat("  Senescence saved\n")

# ==========================================
# 4. Metabolic reprogramming (glycolysis vs OXPHOS)
# ==========================================
cat("\n=== Metabolic Reprogramming ===\n)")
metab_genes <- list(
  Glycolysis = c("HK1","HK2","HK3","GPI","PFKL","PFKP","ALDOA","GAPDH","PGK1","PGAM1","ENO1","PKM","LDHA","LDHB"),
  OXPHOS = c("NDUFA1","NDUFB1","SDHA","SDHB","UQCRB","COX5A","COX5B","ATP5A1","ATP5B","ATP5F1"),
  FAO = c("CPT1A","CPT2","ACADM","ACADS","ACADVL","HADHA","HADHB"),
  PPP = c("G6PD","PGD","TKT","TALDO1")
)
metab_scores <- sapply(metab_genes, function(g) {
  genes <- intersect(g, rownames(expr))
  if (length(genes) >= 3) colMeans(expr[genes,]) else rep(NA, ncol(expr))
})
metab_scores <- metab_scores[, colSums(is.na(metab_scores)) == 0]

metab_cor <- sapply(colnames(metab_scores), function(pw) {
  cor(sp110_expr, metab_scores[,pw], method="spearman")
})
metab_df <- data.frame(Pathway=names(metab_cor), Rho=metab_cor)
metab_df$Pathway <- factor(metab_df$Pathway, levels=metab_df$Pathway[order(metab_df$Rho)])

pdf(file.path(pdf_dir, "Fig_M11_Metabolism.pdf"), width=6, height=4)
ggplot(metab_df, aes(x=Pathway, y=Rho, fill=Rho>0)) +
  geom_bar(stat="identity") + coord_flip() +
  scale_fill_manual(values=c("TRUE"="#E74C3C","FALSE"="#3498DB"), guide="none") +
  labs(title="SP110 vs Metabolic Pathways", x="", y="Spearman rho")
dev.off(); cat("  Metabolism saved\n")

# ==========================================
# 5. TF binding site enrichment in SP110/ER stress genes
# ==========================================
cat("\n=== TF Enrichment ===\n")
# Known SP110-regulating TFs and their targets
tf_targets <- list(
  STAT1 = c("GBP1","GBP2","GBP5","IRF1","IDO1","WARS","TAP1"),
  IRF1 = c("GBP1","GBP2","STAT1","TAP1","PSMB8","PSMB9"),
  NFKB1 = c("IL1B","TNF","CCL2","NLRP3","BCL2"),
  ATF4 = c("DDIT3","PPP1R15A","TRIB3","ASNS","SESN2"),
  XBP1 = c("HSPA5","ERN1","PDIA3","DNAJB9","EDEM1"),
  CEBPB = c("IL6","TNF","SERPINE1","CCL2","IL1B"),
  TP53 = c("CDKN1A","BAX","BBC3","GADD45A","RRM2"),
  HIF1A = c("VEGFA","LDHA","HK2","SLC2A1","BNIP3")
)
tf_scores <- sapply(tf_targets, function(tg) {
  genes <- intersect(tg, rownames(expr))
  if (length(genes) >= 3) colMeans(expr[genes,]) else rep(NA, ncol(expr))
})
tf_scores <- tf_scores[, colSums(is.na(tf_scores)) == 0]

tf_cor_sp110 <- sapply(colnames(tf_scores), function(tf) cor(sp110_expr, tf_scores[,tf], method="spearman"))
tf_cor_tb <- sapply(colnames(tf_scores), function(tf) {
  cor(as.numeric(tb), tf_scores[,tf], method="spearman")
})
tf_df <- data.frame(TF=names(tf_cor_sp110), SP110_rho=tf_cor_sp110, TB_rho=tf_cor_tb)
tf_df$TF <- factor(tf_df$TF, levels=tf_df$TF[order(tf_df$SP110_rho)])

pdf(file.path(pdf_dir, "Fig_M11_TF_Regulon.pdf"), width=7, height=5)
ggplot(tf_df, aes(x=TF, y=SP110_rho, fill=SP110_rho>0)) +
  geom_bar(stat="identity") + coord_flip() +
  scale_fill_manual(values=c("TRUE"="#E74C3C","FALSE"="#3498DB"), guide="none") +
  labs(title="TF Regulon Activity vs SP110", x="", y="Spearman rho (SP110)")
dev.off(); cat("  TF regulon saved\n")

cat(sprintf("\nM11 complete. ~10 PANoptosis/senescence/metabolism/TF figures\n"))
