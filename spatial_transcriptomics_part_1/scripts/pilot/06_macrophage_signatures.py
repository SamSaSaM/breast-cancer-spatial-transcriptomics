import numpy as np
import pandas as pd
import scanpy as sc
from common import *

p=RESULTS/f'{SAMPLE}_normalized_spatial.h5ad'
if not p.exists(): raise FileNotFoundError('Run 02 first')
adata=sc.read_h5ad(p)
signatures={
    'macrophage_core':['C1QA','C1QB','C1QC','CD68','CSF1R'],
    'antigen_presentation':['HLA-DRA','HLA-DRB1','CD74'],
    'inflammatory_myeloid':['IL1B','CXCL8','TNF'],
    'lipid_associated_myeloid':['APOE','TREM2','SPP1'],
    'immunoregulatory_myeloid':['CD274','PDCD1LG2','MRC1'],
}
rows=[]
for name, genes in signatures.items():
    present=[g for g in genes if g in adata.var_names]
    rows.append({'signature':name,'available_genes':','.join(present),'missing_genes':','.join(g for g in genes if g not in present)})
    if len(present)<2:
        print(f'Skipping {name}: only {len(present)} marker(s) available');continue
    # Mean log-normalized expression: simple descriptive signature, not cell abundance.
    mat=adata[:,present].X
    score=np.asarray(mat.mean(axis=1)).ravel()
    adata.obs[name]=score
    plot_spatial(adata,score,title=f'{SAMPLE}: {name} signature',outfile=FIG/f'S6_{SAMPLE}_{name}.png')
pd.DataFrame(rows).to_csv(RESULTS/f'{SAMPLE}_signature_gene_coverage.csv',index=False)
cols=[k for k in signatures if k in adata.obs]
if cols: adata.obs[cols].to_csv(RESULTS/f'{SAMPLE}_signature_scores.csv')
print('S6 signature maps saved. Scores are not deconvolved macrophage proportions.')
