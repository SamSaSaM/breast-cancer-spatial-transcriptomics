import numpy as np
import pandas as pd
import scanpy as sc
from common import *

p=RESULTS/f'{SAMPLE}_normalized_spatial.h5ad'
if not p.exists(): raise FileNotFoundError('Run 02 first')
adata=sc.read_h5ad(p)
# Use count layer for HVG detection, without dense conversion.
sc.pp.highly_variable_genes(adata,n_top_genes=min(2000,adata.n_vars),flavor='seurat',subset=False)
if adata.var['highly_variable'].sum()<30: raise ValueError('Too few HVGs')
work=adata[:,adata.var['highly_variable']].copy()
sc.pp.pca(work,n_comps=min(30,work.n_obs-1,work.n_vars-1),zero_center=False)
sc.pp.neighbors(work,n_neighbors=min(15,work.n_obs-1),n_pcs=min(30,work.obsm['X_pca'].shape[1]))
sc.tl.leiden(work,resolution=0.5,key_added='leiden',flavor='igraph',directed=False)
sc.tl.umap(work)
adata.obs['leiden']=work.obs['leiden'].copy()
adata.obsm['X_umap']=work.obsm['X_umap'].copy()
adata.write_h5ad(RESULTS/f'{SAMPLE}_processed_spatial.h5ad',compression='gzip')
labels=adata.obs['leiden'].astype('category')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
im,coords,sf=image_and_coords(adata)
radius=float(sf['spot_diameter_fullres'])*float(sf['tissue_hires_scalef'])/2
fig,ax=plt.subplots(figsize=(9,9)); ax.imshow(im)
cats=labels.cat.categories
palette=plt.get_cmap('tab20',len(cats))
colors=[palette(int(labels.cat.codes.iloc[i])) for i in range(len(labels))]
patches=PatchCollection([Circle((x,y),radius) for x,y in coords],facecolors=colors,edgecolor='none',alpha=.88)
ax.add_collection(patches)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=palette(i),label=str(c)) for i,c in enumerate(cats)],title='Leiden',loc='upper left',bbox_to_anchor=(1,1),fontsize=8)
ax.set_xlim(0,im.shape[1]);ax.set_ylim(im.shape[0],0);ax.axis('off');ax.set_title(f'{SAMPLE}: spatial Leiden clusters')
fig.tight_layout();fig.savefig(FIG/f'S4_{SAMPLE}_spatial_leiden.png',dpi=180,bbox_inches='tight');plt.close(fig)
sc.pl.umap(adata,color='leiden',show=False,return_fig=True).savefig(FIG/f'S4_{SAMPLE}_umap_leiden.png',dpi=180,bbox_inches='tight')
plt.close('all')
adata.obs['leiden'].value_counts().sort_index().to_csv(RESULTS/f'{SAMPLE}_cluster_counts.csv')
print(f'S4 saved: {len(cats)} spatial clusters')
