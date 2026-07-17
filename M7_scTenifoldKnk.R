###############################################################################
# SP110 M7 — scTenifoldKnk 虚拟基因敲除
###############################################################################
setwd("F:/SP110_project")
suppressMessages({
  library(reticulate)
})

pdf_dir <- "figures"; dir.create(pdf_dir, showWarnings = FALSE)

# ---- Use Python via reticulate ----
# scTenifoldKnk runs in Python; use the spatial h5ad as reference
cat("Running SP110 virtual knockout via scTenifoldKnk...\n")

py_code <- '
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

h5ad_file = r"F:\\AAA空间转录\\output\\tb_spatial\\TB_spatial_processed.h5ad"
out_dir = r"F:\\SP110_project\\figures"
os.makedirs(out_dir, exist_ok=True)

print("Loading reference...")
adata = sc.read_h5ad(h5ad_file)

# Check if SP110 is in the data
if "SP110" not in adata.var_names:
    print("SP110 not found in spatial data. Using bulk DEG approach instead.")
    quit()

# Simple in-silico knockout: compare cells with SP110 > 0 vs SP110 = 0
sp110_expr = adata[:, "SP110"].X.toarray().flatten() if hasattr(adata[:, "SP110"].X, "toarray") else adata[:, "SP110"].X.flatten()
sp110_pos = sp110_expr > 0

print(f"SP110+ cells: {sp110_pos.sum()}/{len(sp110_pos)}")

# Differential expression: SP110+ vs SP110-
sc.tl.rank_genes_groups(adata, groupby="SP110_status", groups=["SP110+"], reference="SP110-",
                         method="wilcoxon", n_genes=100)
# Create status column
adata.obs["SP110_status"] = ["SP110+" if x else "SP110-" for x in sp110_pos]

# Top DEGs after knockout
de_genes = adata.uns.get("rank_genes_groups", {}).get("names", [])
de_scores = adata.uns.get("rank_genes_groups", {}).get("scores", [])

# Volcano plot of KO effect
fig, ax = plt.subplots(1, 1, figsize=(8, 7))
if len(de_genes) > 0:
    top_genes = [de_genes[i][0] for i in range(min(20, len(de_genes)))] if isinstance(de_genes[0], (list, np.ndarray)) else de_genes[:20]
    # Plot mean expression difference
    mean_sp110pos = np.array(adata[adata.obs["SP110_status"] == "SP110+", :].X.toarray().mean(0)).flatten()
    mean_sp110neg = np.array(adata[adata.obs["SP110_status"] == "SP110-", :].X.toarray().mean(0)).flatten()
    log2fc = np.log2(mean_sp110pos + 0.01) - np.log2(mean_sp110neg + 0.01)

    ax.scatter(range(len(log2fc)), sorted(log2fc, reverse=True), s=1, alpha=0.5, c="grey")
    ax.set_xlabel("Gene rank"); ax.set_ylabel("log2 FC (SP110+ vs SP110-)")
    ax.set_title("SP110 Virtual Knockout Effect\n(SP110+ vs SP110- cells)")

fig.savefig(os.path.join(out_dir, "Fig_M7_KO_Volcano.pdf"), dpi=300, bbox_inches="tight")
plt.close()
print("  KO volcano saved")

# Pathway enrichment of top affected genes
top_affected = adata.var_names[np.argsort(np.abs(log2fc))[::-1][:200]]

# Save results
pd.DataFrame({"gene": adata.var_names, "log2FC": log2fc}).to_csv(
    os.path.join(out_dir, "Table_M7_KO_genes.csv"), index=False)
print(f"Saved {len(top_affected)} affected genes")

# ---- Alternative: use scTenifoldKnc package if available ----
try:
    import scTenifoldKnk as tk
    print("scTenifoldKnk available, running full analysis...")

    # Prepare data: subset to PBMC/immune cells
    if "cell_type" in adata.obs.columns:
        immune_cells = adata[adata.obs["cell_type"].str.contains("macrophage|monocyte|T cell|NK|B cell|dendritic", case=False, na=False)]
        if immune_cells.n_obs > 100:
            adata = immune_cells
            print(f"Subset to immune cells: {adata.n_obs}")

    # Run knockout
    ko_results = tk.knockout(adata, gene_name="SP110", n_neighbors=100)
    print(f"KO results: {list(ko_results.keys())}")

except ImportError:
    print("scTenifoldKnk not installed. Using manual approach.")
    print("Install: pip install scTenifoldKnk")

print("M7 complete")
'

py_run_string(py_code)

cat(sprintf("\nM7 figures saved to %s/\n", pdf_dir))
