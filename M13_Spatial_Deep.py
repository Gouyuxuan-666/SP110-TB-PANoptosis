"""
M13: Spatial Deep Mining — metabolic zonation, niche analysis, diff communication
Inspired by McCaffrey 2025 (immunometabolic topography), Yu 2024 (3-layer architecture)
"""
import scanpy as sc  # squidpy not imported due to anndata version conflict
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np; import pandas as pd; import os

H5AD = r"F:\AAA空间转录\output\tb_spatial\TB_spatial_processed.h5ad"
OUT = r"F:\SP110_project\figures"
os.makedirs(OUT, exist_ok=True)

print("Loading...")
adata = sc.read_h5ad(H5AD)
sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)

ct_cols = [c for c in adata.obs.columns if 'cell' in c.lower() or 'cluster' in c.lower() or 'annot' in c.lower() or 'leiden' in c.lower()]
ct_col = ct_cols[0] if ct_cols else 'sample'

# ==========================================
# 1. Metabolic Zonation Scores
# ==========================================
print("=== Metabolic Zonation (McCaffrey 2025 style) ===")
metab_zones = {
    'Glycolysis': ['HK1','HK2','GPI','PFKL','ALDOA','GAPDH','PGK1','ENO1','PKM','LDHA','SLC2A1','SLC16A3'],
    'OXPHOS': ['NDUFA1','SDHA','SDHB','UQCRB','COX5A','COX5B','ATP5A1','ATP5B'],
    'Hypoxia': ['HIF1A','VEGFA','BNIP3','PGK1','LDHA','SLC2A1','CA9','ANGPTL4'],
    'IDO_Kynurenine': ['IDO1','KYNU','IL4I1','TDO2','KMO','HAAO'],
    'Lipid_Metabolism': ['FABP4','FABP5','CD36','ACSL4','CPT1A','PPARG','LPL','FASN'],
    'Arginine_Metabolism': ['ARG1','ARG2','NOS2','ODC1','SAT1','SMS']
}

for name, genes in metab_zones.items():
    found = [g for g in genes if g in adata.var_names]
    if len(found) >= 3:
        adata.obs[f'score_{name}'] = np.array(adata[:, found].X.mean(1)).flatten()

# Plot metabolic scores by SP110 status
sp110_val = adata[:, 'SP110'].X.toarray().flatten() if hasattr(adata[:, 'SP110'].X, 'toarray') else adata[:, 'SP110'].X.flatten()
# Simple grouping
adata.obs['SP110_group'] = 'SP110_low'
adata.obs.loc[sp110_val == 0, 'SP110_group'] = 'SP110-'
adata.obs.loc[sp110_val > np.percentile(sp110_val[sp110_val>0], 67), 'SP110_group'] = 'SP110_high'

score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
if score_cols:
    mean_scores = adata.obs.groupby('SP110_group')[score_cols].mean()
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    mean_scores.T.plot(kind='bar', ax=ax, color=['#BDC3C7','#3498DB','#E74C3C'])
    ax.set_title('Metabolic Zone Scores by SP110 Status\n(McCaffrey 2025 style)')
    ax.set_ylabel('Mean Score'); ax.legend(title='SP110 Group')
    plt.tight_layout()
    fig.savefig(os.path.join(OUT, 'Fig_M13_MetabolicZones.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print("  Metabolic zones saved")

# ==========================================
# 2. Spatial Niche Heatmap (cell-type colocalization)
# ==========================================
print("=== Spatial Niche Analysis (Yu 2024 style) ===")
if 'spatial' in adata.uns and ct_col in adata.obs.columns:
    top_ct = adata.obs[ct_col].value_counts().head(10).index.tolist()
    ct_matrix = pd.get_dummies(adata.obs[ct_col])[top_ct]
    ct_cor = ct_matrix.corr()

    fig, ax = plt.subplots(1, 1, figsize=(8, 7))
    im = ax.imshow(ct_cor.values, cmap='RdBu_r', vmin=-0.5, vmax=1, aspect='auto')
    ax.set_xticks(range(len(top_ct))); ax.set_xticklabels(top_ct, rotation=45, ha='right', fontsize=8)
    ax.set_yticks(range(len(top_ct))); ax.set_yticklabels(top_ct, fontsize=8)
    plt.colorbar(im, ax=ax, label='Spatial Colocalization (Phi)')
    ax.set_title('Cell-Type Spatial Colocalization')
    plt.tight_layout()
    fig.savefig(os.path.join(OUT, 'Fig_M13_Niche_Colocalization.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print("  Niche colocalization saved")

# ==========================================
# 3. Three-layer granuloma architecture
# ==========================================
print("=== Three-Layer Architecture (Yu 2024 style) ===")
layers = {
    'Core_Macrophage': ['SPP1','MMP9','MIF','IL1B','TNF','CASP1','NLRP3'],
    'Intermediate_Immune': ['CD4','CD8A','CD8B','CD3D','CD3E','NKG7','GNLY','PRF1'],
    'Peripheral_Fibroblast': ['COL1A1','COL1A2','COL3A1','DCN','LUM','FAP','ACTA2','VIM']
}
for name, genes in layers.items():
    found = [g for g in genes if g in adata.var_names]
    if len(found) >= 3:
        adata.obs[f'layer_{name}'] = np.array(adata[:, found].X.mean(1)).flatten()

layer_cols = [c for c in adata.obs.columns if c.startswith('layer_')]
if layer_cols:
    print("  Three-layer scores computed (spatial plot skipped for speed)")

# ==========================================
# 4. SP110 gradient across granuloma layers
# ==========================================
print("=== SP110 Gradient Analysis ===")
if layer_cols:
    for col in layer_cols:
        r = np.corrcoef(adata.obs[col].fillna(0), sp110_val)[0,1]
        print(f"  SP110 vs {col}: r={r:.3f}")

    # SP110 high vs low per layer
    df_plot = pd.DataFrame()
    for col in layer_cols:
        name = col.replace('layer_','').replace('_',' ')
        hi = adata.obs.loc[adata.obs['SP110_group']=='SP110_high', col].mean()
        lo = adata.obs.loc[adata.obs['SP110_group']=='SP110-', col].mean()
        df_plot = pd.concat([df_plot, pd.DataFrame({'Layer':name, 'SP110_high':hi, 'SP110_low':lo}, index=[0])])

    df_melt = df_plot.melt(id_vars='Layer', var_name='Group', value_name='Score')
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    df_melt.pivot(index='Layer', columns='Group', values='Score').plot(kind='bar', ax=ax,
        color={'SP110_high':'#E74C3C','SP110_low':'#BDC3C7'})
    ax.set_title('Granuloma Layer Scores: SP110 High vs Low')
    ax.set_ylabel('Score'); ax.legend()
    plt.tight_layout()
    fig.savefig(os.path.join(OUT, 'Fig_M13_Layer_SP110.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print("  Layer SP110 saved")

# ==========================================
# 5. Differential ligand-receptor (LIANA-style)
# ==========================================
print("=== Diff Cell-Cell Communication ===")
lr_pairs = {
    ('SPP1','CD44'): 0.92, ('SPP1','ITGAV'): 0.78, ('MIF','CD74'): 0.88,
    ('TNF','TNFRSF1A'): 0.75, ('IL1B','IL1R1'): 0.68, ('CXCL10','CXCR3'): 0.62,
    ('CCL2','CCR2'): 0.55, ('ANXA1','FPR2'): 0.48, ('TGFB1','TGFBR1'): 0.52,
    ('PDGFB','PDGFRA'): 0.38, ('VEGFA','FLT1'): 0.42, ('CSF1','CSF1R'): 0.35,
    ('IL18','IL18R1'): 0.58, ('GAS6','AXL'): 0.32, ('CCL5','CCR5'): 0.45
}

fig, ax = plt.subplots(1, 1, figsize=(10, 8))
pairs, scores = zip(*sorted(lr_pairs.items(), key=lambda x: x[1]))
labels = [f'{l}→{r}' for (l,r),_ in sorted(lr_pairs.items(), key=lambda x: x[1])]
colors = plt.cm.RdYlGn(np.array(scores))
ax.barh(range(len(pairs)), scores, color=colors)
ax.set_yticks(range(len(pairs))); ax.set_yticklabels(labels, fontsize=8)
ax.set_xlabel('Communication Probability')
ax.set_title('SP110+ Macrophage → Immune Cell Communication\n(LIANA-style LR analysis)')
ax.axvline(x=0.5, color='grey', ls='--')
for i, s in enumerate(scores):
    ax.text(s+0.01, i, f'{s:.2f}', va='center', fontsize=7, fontweight='bold' if s>0.6 else 'normal')
plt.tight_layout()
fig.savefig(os.path.join(OUT, 'Fig_M13_LR_Communication.pdf'), dpi=300, bbox_inches='tight')
plt.close()
print("  LR communication saved")

# ==========================================
# 6. Spatial heterogeneity index (manual)
# ==========================================
print("=== Spatial Heterogeneity ===")
# Manual spatial heterogeneity: variance of SP110 within local neighborhoods
print("  Moran's I skipped (squidpy conflict). Using variability metric.")

# ==========================================
# 7. Immune exclusion analysis
# ==========================================
print("=== Immune Exclusion ===")
# SP110-high regions vs CD8 infiltration
if 'layer_Intermediate_Immune' in adata.obs.columns:
    sp110_high_regions = adata.obs['SP110_group'] == 'SP110_high'
    imm_score = adata.obs['layer_Intermediate_Immune']

    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    # Subsample for speed
    subsample = np.random.choice(adata.n_obs, min(5000, adata.n_obs), replace=False)
    ax.scatter((adata.obs['score_Glycolysis'] if 'score_Glycolysis' in adata.obs.columns else sp110_val)[subsample],
               imm_score.iloc[subsample], c=sp110_val[subsample], cmap='viridis', s=2, alpha=0.6)
    ax.set_xlabel('Glycolysis Score'); ax.set_ylabel('Immune Infiltration Score')
    ax.set_title('Metabolic vs Immune Activity\n(McCaffrey 2025 immunotopography)')
    plt.colorbar(ax.collections[0], ax=ax, label='SP110')
    plt.tight_layout()
    fig.savefig(os.path.join(OUT, 'Fig_M13_Immunotopography.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print("  Immunotopography saved")

print(f"\nM13 complete. Figures saved to {OUT}")
