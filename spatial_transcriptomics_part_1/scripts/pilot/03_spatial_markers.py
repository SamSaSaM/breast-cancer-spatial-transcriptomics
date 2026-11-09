import scanpy as sc
from common import *

p=RESULTS/f'{SAMPLE}_normalized_spatial.h5ad'
if not p.exists(): raise FileNotFoundError('Run 02 first')
adata=sc.read_h5ad(p)
for gene in ['EPCAM','KRT19','COL1A1','DCN','PTPRC','CD3D','LYZ','C1QA']:
    vec=gene_vector(adata,gene)
    if vec is None:
        print(f'Skipping missing gene {gene}'); continue
    plot_spatial(adata,vec,title=f'{SAMPLE}: {gene} (log-normalized)',outfile=FIG/f'S3_{SAMPLE}_{gene}.png')
print('S3 spatial marker plots saved.')
