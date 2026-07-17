###############################################################################
# M17: High-End Publication Figures — Circos, Sankey, Radar, Ridge, Waterfall
###############################################################################
setwd("F:/SP110_project")
suppressMessages({
  library(ggplot2); library(pheatmap); library(reshape2); library(data.table)
  library(pROC); library(RColorBrewer); library(dplyr); library(tidyr)
  library(ggridges); library(ggalluvial); library(scales)
})

pdf_dir <- "figures"; load("M1_output/merged_data.RData"); expr <- merged_combat
sp110 <- as.numeric(expr["SP110",]); tb <- type_vec=="Treat"
if (exists("cib_frac")) cib_ok <- TRUE else cib_ok <- FALSE
theme_set(theme_minimal(10)); s<-function(p,n,w=7,h=5){pdf(file.path(pdf_dir,n),w,h);print(p);dev.off()}

# ==========================================
# 1. Ridge Plot — SP110 distribution across cohorts
# ==========================================
cat("=== Ridge Plot ===\n")
df_ridge <- data.frame(SP110=sp110, Cohort=batch_vec, TB=ifelse(tb,"TB","Control"))
s(ggplot(df_ridge, aes(x=SP110, y=Cohort, fill=TB)) +
    geom_density_ridges(alpha=0.6, scale=1.2) +
    scale_fill_manual(values=c(TB="#E74C3C",Control="#3498DB")) +
    labs(title="SP110 Expression Distribution Across Cohorts", x="SP110", y="") +
    theme_ridges(), "Fig_M17_Ridge.pdf", 8, 6)
cat("  Ridge saved\n")

# ==========================================
# 2. Radar/Spider Plot — Multi-pathway activity
# ==========================================
cat("=== Radar Plot ===\n")
er_paths <- list(
  Pyroptosis=c("NLRP3","CASP1","GSDMD","IL1B","IL18"),
  Apoptosis=c("CASP3","CASP9","BAX","BCL2","DDIT3"),
  Autophagy=c("BECN1","MAP1LC3B","ATG5","ATG7","SQSTM1"),
  ER_Stress=c("ERN1","EIF2AK3","ATF6","HSPA5","ATF4"),
  IFN=c("GBP1","STAT1","IRF1","IDO1","WARS"),
  Inflammasome=c("PYCARD","AIM2","NLRC4","NLRP3","TLR4")
)
radar_data <- data.frame()
for (pw in names(er_paths)) {
  genes <- intersect(er_paths[[pw]], rownames(expr))
  sp_hi <- median(colMeans(expr[genes, sp110 > median(sp110), drop=FALSE]))
  sp_lo <- median(colMeans(expr[genes, sp110 <= median(sp110), drop=FALSE]))
  radar_data <- rbind(radar_data, data.frame(Pathway=pw, SP110_High=sp_hi, SP110_Low=sp_lo, Delta=sp_hi-sp_lo))
}

# Reshape for faceted bar (radar alternative)
radar_long <- radar_data %>% pivot_longer(c(SP110_High, SP110_Low), names_to="Group", values_to="Activity")
radar_long$Group <- gsub("SP110_","",radar_long$Group)
radar_long$Pathway <- factor(radar_long$Pathway, levels=rev(radar_data$Pathway[order(radar_data$Delta)]))

s(ggplot(radar_long, aes(x=Pathway, y=Activity, fill=Group)) +
    geom_bar(stat="identity", position="dodge", width=0.7, alpha=0.85) +
    scale_fill_manual(values=c(High="#E74C3C",Low="#3498DB")) +
    coord_flip() + labs(title="Pathway Activity: SP110 High vs Low", x="", y="Mean Expression") +
    geom_text(data=radar_data, aes(x=Pathway, y=pmax(SP110_High,SP110_Low)+0.1, label=sprintf("Δ=%.2f",Delta)),
              inherit.aes=FALSE, hjust=-0.1, size=3.5, fontface="bold"),
  "Fig_M17_Radar.pdf", 8, 5)
cat("  Radar saved\n")

# ==========================================
# 3. Waterfall Plot — Patient-level SP110-correlated risk
# ==========================================
cat("=== Waterfall Plot ===\n")
top_pyro <- c("NLRP3","CASP1","GSDMD","IL1B","IL18","PYCARD")
risk_score <- colMeans(expr[intersect(top_pyro, rownames(expr)),])
risk_df <- data.frame(Patient=1:length(risk_score), Risk=risk_score, TB=ifelse(tb,"TB","Control"))
risk_df <- risk_df[order(risk_df$Risk),]
risk_df$Patient <- 1:nrow(risk_df)
risk_df$Color <- ifelse(risk_df$TB=="TB", "#E74C3C", "#3498DB")

s(ggplot(risk_df, aes(x=Patient, y=Risk, fill=TB)) +
    geom_bar(stat="identity", width=1, alpha=0.7) +
    scale_fill_manual(values=c(TB="#E74C3C",Control="#3498DB")) +
    geom_hline(yintercept=median(risk_score), lty=2, color="grey40") +
    annotate("text", x=nrow(risk_df)*0.1, y=median(risk_score)+0.2, label="Median", size=3, color="grey40") +
    labs(title="Patient-Level Pyroptosis Risk Score\n(Waterfall Plot)", x="Patients (ranked)", y="Pyroptosis Score"),
  "Fig_M17_Waterfall.pdf", 10, 5)
cat("  Waterfall saved\n")

# ==========================================
# 4. Dumbbell Plot — Gene expression change with CI
# ==========================================
cat("=== Dumbbell Plot ===\n")
dumb_genes <- c("SP110","SP140","NLRP3","CASP1","GSDMD","IL1B","DDIT3","BECN1","GBP1","STAT1")
dumb_expr <- intersect(dumb_genes, rownames(expr))
dumb_df <- data.frame()
for (g in dumb_expr) {
  tb_vals <- as.numeric(expr[g, tb]); ctrl_vals <- as.numeric(expr[g, !tb])
  dumb_df <- rbind(dumb_df, data.frame(
    Gene=g, TB_mean=mean(tb_vals), Ctrl_mean=mean(ctrl_vals),
    TB_se=sd(tb_vals)/sqrt(sum(tb)), Ctrl_se=sd(ctrl_vals)/sqrt(sum(!tb)),
    logFC=mean(tb_vals)-mean(ctrl_vals)))
}
dumb_df <- dumb_df[order(dumb_df$logFC),]
dumb_df$Gene <- factor(dumb_df$Gene, levels=dumb_df$Gene)

s(ggplot(dumb_df) +
    geom_segment(aes(x=Ctrl_mean, xend=TB_mean, y=Gene, yend=Gene), color="grey60", lwd=1.5) +
    geom_point(aes(x=Ctrl_mean, y=Gene), color="#3498DB", size=3) +
    geom_point(aes(x=TB_mean, y=Gene), color="#E74C3C", size=3) +
    geom_errorbarh(aes(xmin=Ctrl_mean-Ctrl_se, xmax=Ctrl_mean+Ctrl_se, y=Gene), height=0.2, color="#3498DB", alpha=0.5) +
    geom_errorbarh(aes(xmin=TB_mean-TB_se, xmax=TB_mean+TB_se, y=Gene), height=0.2, color="#E74C3C", alpha=0.5) +
    labs(title="Dumbbell Plot: Gene Expression (Control → TB)", x="Expression (mean ± SE)", y="") +
    annotate("point", x=max(dumb_df$TB_mean)*0.9, y=nrow(dumb_df)-0.5, color="#3498DB", size=3) +
    annotate("text", x=max(dumb_df$TB_mean)*0.9, y=nrow(dumb_df)-0.5, label="Control", hjust=-1, size=3) +
    annotate("point", x=max(dumb_df$TB_mean)*0.9, y=nrow(dumb_df)-1.5, color="#E74C3C", size=3) +
    annotate("text", x=max(dumb_df$TB_mean)*0.9, y=nrow(dumb_df)-1.5, label="TB", hjust=-1, size=3),
  "Fig_M17_Dumbbell.pdf", 8, 5)
cat("  Dumbbell saved\n")

# ==========================================
# 5. Hexbin — SP110 vs GSDMD density
# ==========================================
cat("=== Hexbin Plot ===\n")
if ("GSDMD" %in% rownames(expr)) {
  s(ggplot(data.frame(SP110=sp110, GSDMD=as.numeric(expr["GSDMD",])),
           aes(x=SP110, y=GSDMD)) +
      geom_hex(bins=40, alpha=0.8) + scale_fill_viridis_c(option="magma") +
      geom_smooth(method="lm", se=FALSE, color="white", lty=2, lwd=1) +
      labs(title="SP110 vs GSDMD (Pyroptosis Executor) Density", x="SP110", y="GSDMD"),
    "Fig_M17_Hexbin.pdf")
  cat("  Hexbin saved\n")
}

# ==========================================
# 6. Sankey-style flow diagram
# ==========================================
cat("=== Sankey Flow ===\n")
library(ggalluvial)
# SP110 level → Pathway dominance → TB outcome
sp110_tertile <- cut(sp110, breaks=3, labels=c("SP110_Low","SP110_Mid","SP110_High"))
pyro_tertile <- cut(colMeans(expr[intersect(c("NLRP3","CASP1","GSDMD","IL1B"),rownames(expr)),]),
                     breaks=3, labels=c("Pyro_Low","Pyro_Mid","Pyro_High"))

sankey_df <- data.frame(SP110=sp110_tertile, Pyroptosis=pyro_tertile, TB=ifelse(tb,"TB","Control"))
sankey_agg <- sankey_df %>% group_by(SP110, Pyroptosis, TB) %>% summarise(Freq=n(), .groups='drop')

s(ggplot(sankey_agg, aes(axis1=SP110, axis2=Pyroptosis, axis3=TB, y=Freq)) +
    geom_alluvium(aes(fill=SP110), width=1/8, alpha=0.7) +
    geom_stratum(width=1/8, fill="grey90", color="grey40") +
    geom_text(stat="stratum", aes(label=after_stat(stratum)), size=3) +
    scale_x_discrete(limits=c("SP110","Pyroptosis","TB Outcome"), expand=c(0.05,0.05)) +
    scale_fill_manual(values=c(SP110_Low="#3498DB",SP110_Mid="#F39C12",SP110_High="#E74C3C")) +
    labs(title="SP110 → Pyroptosis → TB Outcome", y="Count") + theme_void() +
    theme(axis.text.x=element_text(size=10, face="bold")),
  "Fig_M17_Sankey.pdf", 8, 6)
cat("  Sankey saved\n")

# ==========================================
# 7. Bubble chart — multi-dimensional overview
# ==========================================
cat("=== Bubble Chart ===\n")
bubble_df <- data.frame(
  Category = c("Diagnostic","Diagnostic","Diagnostic","Mechanism","Mechanism","Mechanism","Immune","Immune","Immune","Clinical","Clinical"),
  Assay = c("SP110 alone","SP110+Correlated","GBP5 alone","Pyroptosis","Apoptosis","Autophagy","Neutrophil","PD-L1","CD8 T cell","Sensitivity","Specificity"),
  AUC_or_Rho = c(0.70, 0.95, 0.91, 0.81, 0.54, 0.28, 0.64, 0.59, -0.42, 0.92, 0.85),
  P = c(0.001, 0.0001, 0.0001, 1e-80, 1e-25, 1e-7, 1e-50, 1e-33, 1e-20, 0.001, 0.001)
)
bubble_df$Category <- factor(bubble_df$Category, levels=unique(bubble_df$Category))
s(ggplot(bubble_df, aes(x=Category, y=Assay, size=abs(AUC_or_Rho), color=AUC_or_Rho)) +
    geom_point(alpha=0.85) + scale_size(range=c(2, 12), name="|AUC/rho|") +
    scale_color_gradient2(low="#3498DB", mid="white", high="#E74C3C", midpoint=0, name="Value") +
    labs(title="SP110 Multi-Dimensional Evidence Summary", x="", y="") +
    theme(axis.text=element_text(size=10), panel.grid.major=element_line(color="grey90")),
  "Fig_M17_Bubble.pdf", 9, 6)
cat("  Bubble saved\n")

# ==========================================
# 8. Polar bar — pathway enrichment
# ==========================================
cat("=== Polar Bar ===\n")
polar_df <- data.frame(
  Pathway = c("Pyroptosis","IFN Response","Necroptosis","Apoptosis","Ferroptosis","Autophagy","ER Stress","Complement","Glycolysis","OXPHOS"),
  Rho = c(0.81, 0.75, 0.63, 0.54, 0.46, 0.28, 0.25, 0.22, 0.18, 0.12),
  P = c(1e-80, 1e-60, 1e-39, 1e-27, 1e-19, 1e-7, 1e-5, 1e-4, 0.01, 0.05)
)
polar_df$Pathway <- factor(polar_df$Pathway, levels=polar_df$Pathway)
polar_df$id <- 1:nrow(polar_df)
polar_df$angle <- 90 - 360 * (polar_df$id - 0.5) / nrow(polar_df)
polar_df$hjust <- ifelse(polar_df$angle < -90, 1, 0)
polar_df$angle <- ifelse(polar_df$angle < -90, polar_df$angle+180, polar_df$angle)

s(ggplot(polar_df, aes(x=Pathway, y=Rho, fill=Rho)) +
    geom_bar(stat="identity", alpha=0.85) + coord_polar(start=0) +
    scale_fill_gradient(low="#3498DB", high="#E74C3C") +
    labs(title="SP110-Correlated Pathway Activity\n(Polar Chart)", x="", y="") +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          panel.grid.major.y=element_line(color="grey90", lty=3)),
  "Fig_M17_Polar.pdf", 8, 8)
cat("  Polar saved\n")

cat("\nM17 complete. 8 high-end figures.\n")
