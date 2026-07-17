"""
M6: SP110 spatial + scRNA-seq in TB granuloma
"""
import scanpy as sc
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os

H5AD = r"F:\AAA空间转录\output\tb_spatial\TB_spatial_processed.h5ad"
OUT = r"F:\SP110_project\figures"
os.makedirs(OUT, exist_ok=True)

sc.settings.figdir = OUT
sc.settings.verbosity = 1

print("Loading data...")
adata = sc.read_h5ad(H5AD)
print(f"Shape: {adata.shape}")

# Find cell type column
ct_cols = [c for c in adata.obs.columns if 'cell' in c.lower() or 'cluster' in c.lower() or 'type' in c.lower() or 'annot' in c.lower() or 'leiden' in c.lower()]
ct_col = ct_cols[0] if ct_cols else 'sample'
print(f"Cell type column: {ct_col}")

# Key genes
genes = ['SP110','SP140','ERN1','EIF2AK3','DDIT3','BECN1','HSPA5',
         'GBP1','GBP5','STAT1','IRF1','CASP3','BAX','BCL2','NLRP3']
genes_found = [g for g in genes if g in adata.var_names]
print(f"Genes found: {genes_found}")

# ==========================================
# Fig M6-1: UMAP with SP110 expression
# ==========================================
if 'X_umap' not in adata.obsm:
    sc.pp.neighbors(adata, n_neighbors=15)
    sc.tl.umap(adata)

fig, axes = plt.subplots(1, 3, figsize=(21, 6))
sc.pl.umap(adata, color=ct_col, ax=axes[0], show=False, title=f'Cell Groups ({ct_col})', legend_loc='right margin')

for i, gene in enumerate(['SP110', 'SP140'], 1):
    if gene in adata.var_names:
        expr = adata[:, gene].X.toarray().flatten() if hasattr(adata[:,gene].X,'toarray') else adata[:,gene].X.flatten()
        vmax = np.percentile(expr[expr>0], 95) if np.sum(expr>0) > 0 else 1
        sc.pl.umap(adata, color=gene, ax=axes[i], show=False, cmap='viridis', vmax=vmax,
                   title=f'{gene} (vmax={vmax:.2f})')
plt.tight_layout()
fig.savefig(os.path.join(OUT, 'Fig_M6_UMAP_SP110.pdf'), dpi=300, bbox_inches='tight')
plt.close()
print("  UMAP saved")

# ==========================================
# Fig M6-2: ER stress gene dotplot by cell type
# ==========================================
er_genes_found = [g for g in ['SP110','SP140','ERN1','EIF2AK3','DDIT3','BECN1','HSPA5','ATF6'] if g in adata.var_names]
if len(er_genes_found) >= 3:
    fig, ax = plt.subplots(1, 1, figsize=(14, 6))
    sc.pl.dotplot(adata, er_genes_found, groupby=ct_col, ax=ax, show=False)
    ax.set_title('ER Stress Genes in TB Granuloma by Cell Type')
    plt.tight_layout()
    fig.savefig(os.path.join(OUT, 'Fig_M6_ER_Dotplot.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print("  ER dotplot saved")

# ==========================================
# Fig M6-3: SP110 vs GBP1 co-expression scatter
# ==========================================
if 'SP110' in adata.var_names and 'GBP1' in adata.var_names:
    sp = adata[:, 'SP110'].X.toarray().flatten() if hasattr(adata[:,'SP110'].X,'toarray') else adata[:,'SP110'].X.flatten()
    gb = adata[:, 'GBP1'].X.toarray().flatten() if hasattr(adata[:,'GBP1'].X,'toarray') else adata[:,'GBP1'].X.flatten()
    mask = (sp > 0) & (gb > 0)
    r = np.corrcoef(sp[mask], gb[mask])[0,1]

    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    ax.scatter(sp, gb, alpha=0.3, s=1, c='#2E75B6')
    ax.set_xlabel('SP110 Expression'); ax.set_ylabel('GBP1 Expression')
    ax.set_title(f'SP110 vs GBP1 in TB Granuloma\nPearson r={r:.3f} (n={mask.sum()} cells)')
    plt.tight_layout()
    fig.savefig(os.path.join(OUT, 'Fig_M6_SP110vsGBP1_scRNA.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  SP110 vs GBP1: r={r:.3f}")

# ==========================================
# Fig M6-4: PERK vs ERN1 pathway gene correlation
# ==========================================
perk_genes = [g for g in ['EIF2AK3','DDIT3','ATF4','PPP1R15A','TRIB3'] if g in adata.var_names]
ern1_genes = [g for g in ['ERN1','BECN1','MAP1LC3B','SQSTM1','XBP1'] if g in adata.var_names]
all_pw = perk_genes + ern1_genes
if len(all_pw) >= 4:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for idx, (label, genes_list) in enumerate([('PERK Pathway', perk_genes), ('ERN1 Pathway', ern1_genes)]):
        valid = [g for g in genes_list if g in adata.var_names]
        if len(valid) >= 2:
            sub = adata[:, valid].X.toarray() if hasattr(adata[:,valid].X,'toarray') else adata[:,valid].X
            corr = np.corrcoef(sub.T)
            im = axes[idx].imshow(corr, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
            axes[idx].set_xticks(range(len(valid))); axes[idx].set_xticklabels(valid, rotation=45, ha='right', fontsize=8)
            axes[idx].set_yticks(range(len(valid))); axes[idx].set_yticklabels(valid, fontsize=8)
            axes[idx].set_title(f'{label} Gene Correlation')
            plt.colorbar(im, ax=axes[idx])
    plt.tight_layout()
    fig.savefig(os.path.join(OUT, 'Fig_M6_PERK_ERN1_Corr.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print("  PERK/ERN1 correlation saved")

# ==========================================
# Fig M6-5: SP110 expression by cell type (violin)
# ==========================================
if 'SP110' in adata.var_names:
    fig, axes = plt.subplots(1, 2, figsize=(16, 5))

    # % expressing
    expr_bool = (adata[:, 'SP110'].X.toarray().flatten() if hasattr(adata[:,'SP110'].X,'toarray') else adata[:,'SP110'].X.flatten()) > 0
    pct = adata.obs[[ct_col]].copy()
    pct['SP110+'] = expr_bool
    pct_by = pct.groupby(ct_col)['SP110+'].mean() * 100
    pct_by = pct_by.sort_values(ascending=False).head(10)

    axes[0].barh(range(len(pct_by)), pct_by.values, color='#2E75B6')
    axes[0].set_yticks(range(len(pct_by)))
    axes[0].set_yticklabels(pct_by.index, fontsize=8)
    axes[0].set_xlabel('% SP110+ cells'); axes[0].invert_yaxis()
    axes[0].set_title('SP110 Detection Rate by Cell Type')
    for i, v in enumerate(pct_by.values):
        axes[0].text(v + 0.5, i, f'{v:.1f}%', va='center', fontsize=7)

    # Violin
    top_ct = pct_by.head(6).index.tolist()
    sub = adata[adata.obs[ct_col].isin(top_ct)]
    sc.pl.stacked_violin(sub, 'SP110', groupby=ct_col, ax=axes[1], show=False, stripplot=True, size=0.3)
    axes[1].set_title('SP110 Expression by Cell Type')
    plt.tight_layout()
    fig.savefig(os.path.join(OUT, 'Fig_M6_SP110_CellType.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print("  SP110 cell type saved")

# ==========================================
# Fig M6-6: Multi-gene spatial expression for top sample
# ==========================================
if 'sample' in adata.obs:
    samples = adata.obs['sample'].unique()[:2]
    top_genes = genes_found[:6]
    for s in samples:
        sub = adata[adata.obs['sample'] == s]
        n = min(len(top_genes), 6)
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        for i, gene in enumerate(top_genes[:n]):
            if gene in sub.var_names:
                row, col = i // 3, i % 3
                sc.pl.spatial(sub, color=gene, ax=axes[row][col], show=False,
                              cmap='viridis', spot_size=1.5, title=gene)
        plt.suptitle(f'Sample: {s}', fontsize=14, fontweight='bold')
        plt.tight_layout()
        fig.savefig(os.path.join(OUT, f'Fig_M6_Spatial_{s}.pdf'), dpi=300, bbox_inches='tight')
        plt.close()
        break  # Just first sample
    print("  Spatial multi-gene saved")

# ==========================================
# Fig M6-7: SP110 vs PD-L1 co-expression
# ==========================================
if 'SP110' in adata.var_names and 'CD274' in adata.var_names:
    sp = adata[:, 'SP110'].X.toarray().flatten() if hasattr(adata[:,'SP110'].X,'toarray') else adata[:,'SP110'].X.flatten()
    cd274 = adata[:, 'CD274'].X.toarray().flatten() if hasattr(adata[:,'CD274'].X,'toarray') else adata[:,'CD274'].X.flatten()
    mask = (sp > 0) | (cd274 > 0)
    r = np.corrcoef(sp[mask], cd274[mask])[0,1]

    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    ax.scatter(sp[mask], cd274[mask], alpha=0.3, s=1, c='#E74C3C')
    ax.set_xlabel('SP110'); ax.set_ylabel('PD-L1 (CD274)')
    ax.set_title(f'SP110 vs PD-L1 in TB Granuloma\nPearson r={r:.3f}')
    plt.tight_layout()
    fig.savefig(os.path.join(OUT, 'Fig_M6_SP110vsPDL1.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  SP110 vs PD-L1: r={r:.3f}")

print(f"\nAll M6 figures saved to {OUT}")
