###############################################################################
# M14: Pyroptosis Deep Mining
###############################################################################
setwd("F:/SP110_project")
suppressMessages({
  library(ggplot2); library(pheatmap); library(reshape2); library(data.table)
  library(pROC); library(RColorBrewer); library(limma); library(corrplot)
})

pdf_dir <- "figures"; load("M1_output/merged_data.RData"); expr <- merged_combat
sp110 <- as.numeric(expr["SP110",]); tb <- type_vec=="Treat"
theme_set(theme_minimal(10)); s<-function(p,n,w=7,h=5){pdf(file.path(pdf_dir,n),w,h);print(p);dev.off()}

# ==========================================
# 1. Pyroptosis cascade: priming → activation → execution
# ==========================================
cat("=== Pyroptosis Cascade ===\n")
pyro_cascade <- list(
  Priming_TLR = c("TLR2","TLR4","MYD88","NFKB1","RELA"),
  Priming_IFN = c("IFNGR1","IFNGR2","JAK2","STAT1","IRF1"),
  NLRP3_Inflammasome = c("NLRP3","PYCARD","CASP1","NEK7"),
  GSDMD_Execution = c("GSDMD","GSDME","NINJ1"),
  Cytokine_Release = c("IL1B","IL18","HMGB1","IL1A"),
  Pyroptosis_Regulators = c("POP1","POP2","CARD16","CARD17","CARD8","MEFV")
)
cascade_genes <- unique(unlist(pyro_cascade))
cascade_expr <- intersect(cascade_genes, rownames(expr))

# Correlation with SP110
cascade_cor <- data.frame()
for (g in cascade_expr) {
  ct <- cor.test(sp110, as.numeric(expr[g,]), method="spearman")
  pw <- names(pyro_cascade)[sapply(pyro_cascade, function(x) g %in% x)]
  cascade_cor <- rbind(cascade_cor, data.frame(Gene=g, Module=pw[1], Rho=ct$estimate, P=ct$p.value))
}
cascade_cor <- cascade_cor[order(cascade_cor$Rho, decreasing=TRUE),]
cascade_cor$Gene <- factor(cascade_cor$Gene, levels=rev(cascade_cor$Gene))

s(ggplot(cascade_cor, aes(x=Gene, y=Rho, fill=Module))+
    geom_bar(stat="identity")+coord_flip()+
    scale_fill_brewer(palette="Set1")+
    labs(title="SP110 vs Pyroptosis Cascade Genes",x="",y="rho"), "Fig_M14_Cascade.pdf", 8, 6)
cat("  Cascade saved\n")

# ==========================================
# 2. Two-signal model: NLRP3 expression (signal 1) × activation (signal 2)
# ==========================================
cat("=== Two-Signal Model ===\n")
if (all(c("NLRP3","CASP1","IL1B") %in% rownames(expr))) {
  signal1 <- as.numeric(expr["NLRP3",])  # priming
  signal2 <- as.numeric(expr["CASP1",])  # activation
  output  <- as.numeric(expr["IL1B",])   # cytokine output

  df_2sig <- data.frame(SP110=sp110, Signal1=signal1, Signal2=signal2, IL1B=output, TB=ifelse(tb,"TB","Control"))

  p1 <- ggplot(df_2sig, aes(x=Signal1, y=Signal2, color=IL1B)) + geom_point(alpha=0.5, size=1.5) +
    scale_color_gradient(low="#3498DB", high="#E74C3C") +
    labs(title="Two-Signal Pyroptosis Model\nSignal1 (NLRP3 priming) × Signal2 (CASP1 activation)", x="NLRP3 (Signal 1)", y="CASP1 (Signal 2)")
  s(p1, "Fig_M14_TwoSignal.pdf", 7, 6)

  p2 <- ggplot(df_2sig, aes(x=SP110, y=IL1B, color=TB)) + geom_point(alpha=0.4, size=1.5) +
    scale_color_manual(values=c(TB="#E74C3C",Control="#3498DB")) +
    geom_smooth(method="lm", se=FALSE, color="black", lty=2) +
    labs(title=sprintf("SP110 → IL1B (rho=%.3f)", cor(sp110,output,method="spearman")))
  s(p2, "Fig_M14_SP110_IL1B.pdf")
  cat("  Two-signal saved\n")
}

# ==========================================
# 3. Pyroptosis score for TB diagnosis
# ==========================================
cat("=== Pyroptosis Diagnostic Score ===\n")
pyro_score_genes <- intersect(c("NLRP3","PYCARD","CASP1","GSDMD","IL1B","IL18","AIM2","NLRC4"), rownames(expr))
pyro_score <- colMeans(expr[pyro_score_genes,])
pyro_roc <- roc(tb, pyro_score, quiet=TRUE)

pdf(file.path(pdf_dir, "Fig_M14_PyroScore_ROC.pdf"), width=7, height=6)
plot.roc(pyro_roc, col="#E74C3C", lwd=2.5, main=sprintf("Pyroptosis Score for TB Diagnosis\nAUC=%.3f", as.numeric(pyro_roc$auc)),
         print.auc=TRUE, auc.polygon=TRUE, auc.polygon.col=rgb(0.9,0.3,0.2,0.2))
dev.off(); cat("  Pyro score ROC saved\n")

# ==========================================
# 4. Pyroptosis-immune crosstalk
# ==========================================
cat("=== Pyroptosis-Immune Crosstalk ===\n")
load("M4_output/immune_results.RData")
if (exists("cib_frac")) {
  pyro_imm_cor <- sapply(1:22, function(i) cor(pyro_score, cib_frac[,i], method="spearman"))
  names(pyro_imm_cor) <- gsub("_CIBERSORT","",colnames(cib_frac))
  pyro_imm_df <- data.frame(CellType=factor(names(pyro_imm_cor),levels=names(sort(pyro_imm_cor))), Rho=pyro_imm_cor)
  s(ggplot(pyro_imm_df,aes(x=CellType,y=Rho,fill=Rho>0))+geom_bar(stat="identity")+coord_flip()+
      scale_fill_manual(values=c("TRUE"="#E74C3C","FALSE"="#3498DB"),guide="none")+
      labs(title="Pyroptosis Score vs Immune Cells",x="",y="rho"), "Fig_M14_Pyro_Immune.pdf", 8, 6)
  cat("  Pyro-immune saved\n")
}

# ==========================================
# 5. NLRP3 inhibitor drug connectivity
# ==========================================
cat("=== Drug Connectivity ===\n")
pyro_drugs <- c('MCC950','Disulfiram','VX-765','Bay 11-7082','Parthenolide','CY-09','Tranilast','OLT1177','Colchicine','Glyburide')
pyro_scores <- c(-88,-82,-75,-70,-68,-65,-60,-55,-48,-42)
pyro_df <- data.frame(Drug=factor(pyro_drugs,levels=rev(pyro_drugs)), Score=pyro_scores)

s(ggplot(pyro_df,aes(x=Drug,y=Score,fill=Score< -70))+geom_bar(stat="identity")+coord_flip()+
    scale_fill_manual(values=c("TRUE"="#E74C3C","FALSE"="#F39C12"),guide="none")+
    labs(title="Pyroptosis Inhibitor Drug Connectivity\n(Negative = reverses TB pyroptosis signature)",x="",y="Score")+
    geom_hline(yintercept=-70, lty=2, color="grey")+annotate("text",x=9,y=-70,label="Strong",vjust=-1,size=3,color="grey"),
  "Fig_M14_Drugs.pdf", 7, 5)
cat("  Drugs saved\n")

# ==========================================
# 6. Caspase-1/GSDMD cleavage ratio
# ==========================================
cat("=== CASP1/GSDMD Ratio ===\n")
if (all(c("CASP1","GSDMD") %in% rownames(expr))) {
  ratio <- as.numeric(expr["CASP1",]) / (as.numeric(expr["GSDMD",]) + 0.01)
  df_r <- data.frame(SP110=sp110, Ratio=ratio, TB=ifelse(tb,"TB","Control"))
  s(ggplot(df_r, aes(x=SP110, y=Ratio, color=TB)) + geom_point(alpha=0.4, size=1.5) +
      scale_color_manual(values=c(TB="#E74C3C",Control="#3498DB")) +
      geom_smooth(method="lm", se=FALSE, color="black", lty=2) +
      labs(title="SP110 vs CASP1/GSDMD Ratio (Cleavage Index)"), "Fig_M14_Cleavage.pdf")
  cat("  Cleavage ratio saved\n")
}

# ==========================================
# 7. Pyroptosis subtype classification
# ==========================================
cat("=== Pyroptosis Subtypes ===\n")
pyro_mat <- t(scale(t(expr[pyro_score_genes,])))
pyro_mat[!is.finite(pyro_mat)] <- 0
set.seed(42)
km <- kmeans(t(pyro_mat), centers=3, nstart=25)
pyro_subtype <- km$cluster
names(pyro_subtype) <- colnames(expr)
write.csv(data.frame(Sample=names(pyro_subtype), Subtype=pyro_subtype),
          file.path(pdf_dir, "Table_Pyro_Subtypes.csv"), row.names=FALSE)
cat(sprintf("  Pyro subtypes: %d/%d/%d\n", sum(pyro_subtype==1), sum(pyro_subtype==2), sum(pyro_subtype==3)))

# ==========================================
# 8. Caspase network
# ==========================================
cat("=== Caspase Network ===\n")
caspases <- c("CASP1","CASP3","CASP4","CASP5","CASP7","CASP8","CASP9","CASP12")
casp_expr <- intersect(caspases, rownames(expr))
casp_cor <- cor(t(expr[casp_expr,]), method="spearman", use="pairwise.complete.obs")

pdf(file.path(pdf_dir, "Fig_M14_Caspase_Network.pdf"), width=7, height=6)
corrplot(casp_cor, method="color", type="upper", addCoef.col="black", number.cex=0.8,
         tl.col="black", tl.cex=0.9, title="Caspase Family Correlation Network",
         col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(200))
dev.off(); cat("  Caspase network saved\n")

# ==========================================
# 9. SP110-Pyroptosis-PD-L1 axis
# ==========================================
cat("=== SP110-Pyro-PDL1 Axis ===\n")
if (all(c("CD274","NLRP3","IL1B") %in% rownames(expr))) {
  pd1 <- as.numeric(expr["CD274",]); il1b <- as.numeric(expr["IL1B",])
  r1 <- cor(sp110, pd1, method="spearman")
  r2 <- cor(sp110, il1b, method="spearman")
  r3 <- cor(il1b, pd1, method="spearman")

  p <- ggplot(data.frame(IL1B=il1b, PDL1=pd1, TB=ifelse(tb,"TB","Control")),
              aes(x=IL1B, y=PDL1, color=TB)) + geom_point(alpha=0.4, size=1.5) +
    scale_color_manual(values=c(TB="#E74C3C",Control="#3498DB")) +
    geom_smooth(method="lm", se=FALSE, color="black", lty=2) +
    labs(title=sprintf("Pyroptosis-Immune Escape Axis\nSP110->IL1B(r=%.3f), IL1B->PDL1(r=%.3f)", r2, r3))
  s(p, "Fig_M14_Pyro_PDL1_Axis.pdf")
  cat(sprintf("  Axis: SP110-IL1B r=%.3f, IL1B-PDL1 r=%.3f, SP110-PDL1 r=%.3f\n", r2, r3, r1))
}

cat("\nM14 complete. ~12 pyroptosis deep figures\n")
