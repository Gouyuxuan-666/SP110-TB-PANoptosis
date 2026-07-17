# SP110-TB-PANoptosis

**SP110 orchestrates PERK/ERN1 double-branch ER stress and a pyroptosis–PD-L1 immune-evasion axis in human tuberculosis**

Computational pipeline for the SP110-TB project. 17 analysis modules covering differential expression, WGCNA, ensemble machine learning, immune deconvolution, PANoptosis pathway scoring, spatial transcriptomics, virtual gene knockout, cell-cell communication, SCENIC TF regulon analysis, molecular docking, and high-end publication figures.

## Repository Structure

| Module | Script | Description |
|--------|--------|-------------|
| M1 | `M1_SP110_DEG.R` | Data integration, ComBat batch correction, differential expression (limma), random-effects meta-analysis |
| M2 | `M2_SP110_WGCNA.R` | Weighted gene co-expression network analysis (WGCNA), module–trait correlation, kME/GS quantification |
| M2/3 Figs | `M2M3_figures.R` | Generation of WGCNA and ML diagnostic figures |
| M3 | `M3_SP110_ML.R` | Ensemble ML: 5 feature selectors (Lasso/Ridge/ElasticNet/StepGLM/SVM-RFE) × 7 classifiers (LDA/RF/GBM/XGBoost/PLS/glmBoost/NB), DeLong test, DCA, Nomogram |
| M4 | `M4_SP110_Immune.R` | CIBERSORT deconvolution (LM22), immune checkpoint correlation, LINKET network |
| M5 | `M5_SP110_Network.R` | Regulatory network, SP110 vs GBP1 comparison, STRING PPI |
| M6 | `M6_SP110_Spatial.py` | Single-cell RNA-seq analysis (Scanpy): QC, normalization, HVG, PCA, UMAP, Leiden clustering, cell type annotation |
| M7 | `M7_KO.py` | Virtual gene knockout simulation |
| M7b | `M7_scTenifoldKnk.R` | scTenifoldKnk virtual knockout with single-cell resolution |
| M8 | `M8_CellChat.py` | Cell–cell communication inference (CellChat) |
| M8b | `M8B_SCENIC.py` | TF regulon inference and kinase activity scoring (pySCENIC) |
| M9 | `M9_DeepAnalysis.R` | PROGENy pathway activity, GSEA preranked, clinical correlation |
| M10 | `M10_FigureFactory.R` | Batch figure generation (40+ publication-ready figures) |
| M10b | `M10_resume.R` | Figure factory resume mode for incremental updates |
| M11 | `M11_PANoptosis.R` | Five-pathway PANoptosis scoring (pyroptosis/apoptosis/autophagy/necroptosis/ferroptosis) via ssGSEA |
| M12 | `M12_PANoptosis_Deep.R` | PANoptosome complex analysis, death subtype classification, GSDMD cleavage index |
| M13 | `M13_Spatial_Deep.py` | Spatial transcriptomics deep analysis: metabolic zonation, 3-layer granuloma architecture |
| M14 | `M14_Pyroptosis_Deep.R` | Pyroptosis cascade dissection: NLRP3→CASP1→GSDMD→IL1B, Caspase co-expression network |
| M15 | `M15_Wetlab_Integration.R` | PERK/ERN1 branch scores, drug perturbation transcriptional signatures, mouse–human ortholog validation, 7/7 closed-loop verification |
| M16 | `M16_HotTopics.R` | Hot-topic analyses: LLPS phase separation, m6A epitranscriptomics (writers/erasers/readers), MALAT1–TLR4 axis, ubiquitination/ISGylation, SP110 domain architecture |
| M17 | `M17_FancyFigs.R` | High-end figures: Ridge, Radar, Waterfall, Dumbbell, Hexbin, Sankey, Bubble, Polar bar |
| Extra | `M_Extra_Figures.R` | Supplementary figures batch |
| Extra | `M_Extra2_Figures.R` | Additional supplementary figures |
| Tools | `tools/run_all_tools.py` | Integrated spatial visualization (pure Scanpy + matplotlib) |

## Dependencies

### R (≥4.3.2)
- **Core**: limma, sva (ComBat), WGCNA, clusterProfiler, GSVA, CIBERSORT
- **ML**: glmnet, ranger, gbm, xgboost, pls, mboost, klaR, caret, pROC
- **Figures**: ggplot2, pheatmap, reshape2, ggalluvial, ggridges, RColorBrewer
- **Systems**: metafor, PROGENy, dorothea, scTenifoldKnk

### Python (≥3.11)
- **Core**: scanpy, anndata, pandas, numpy, scipy
- **Figures**: matplotlib, seaborn
- **Spatial**: squidpy, cellchat (R/Python)

## Data

All human transcriptomic datasets are publicly available at NCBI GEO:
GSE19491, GSE28623, GSE34608, GSE37250, GSE42830, GSE83456, GSE19444, GSE39940, GSE296400

## License

MIT License. See [LICENSE](LICENSE) for details.

## Citation

If you use this code, please cite our manuscript (in preparation) and the individual method papers referenced in each module.
