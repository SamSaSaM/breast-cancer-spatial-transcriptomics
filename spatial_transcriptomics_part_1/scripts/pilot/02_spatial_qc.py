import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import *

adata = load_raw()
adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'] if adata.var['mt'].any() else [], inplace=True, percent_top=None)
if 'pct_counts_mt' not in adata.obs: adata.obs['pct_counts_mt'] = np.nan
cols = ['total_counts','n_genes_by_counts','pct_counts_mt']
fig, axes = plt.subplots(1,3,figsize=(14,4))
for ax, c in zip(axes,cols):
    x = adata.obs[c].dropna()
    ax.hist(x, bins=45, color='#4B7AA3'); ax.set_title(c); ax.set_xlabel(c); ax.set_ylabel('Spots')
fig.suptitle(f'{SAMPLE} spatial spot QC')
fig.tight_layout(); fig.savefig(FIG/f'S2_{SAMPLE}_qc_distributions.png',dpi=180); plt.close(fig)
for c in cols:
    if adata.obs[c].notna().any():
        plot_spatial(adata, adata.obs[c].to_numpy(), title=f'{SAMPLE}: {c}',outfile=FIG/f'S2_{SAMPLE}_{c}_spatial.png')
adata.obs[cols].describe().to_csv(RESULTS/f'{SAMPLE}_qc_summary.csv')
# Keep all in-tissue spots in this exploratory workflow; no arbitrary QC threshold.
adata.layers['counts'] = adata.X.copy()
sc.pp.filter_genes(adata,min_cells=3)
sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata.copy()
adata.write_h5ad(RESULTS/f'{SAMPLE}_normalized_spatial.h5ad',compression='gzip')
print('S2 QC and normalized data saved; no spots filtered.')
