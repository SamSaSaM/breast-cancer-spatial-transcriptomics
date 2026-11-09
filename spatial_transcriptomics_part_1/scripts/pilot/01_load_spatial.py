import gzip
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from scipy import sparse
from anndata import AnnData
from common import *

def open_by_magic(path):
    """Some Zenodo .gz files contain uncompressed Matrix Market text."""
    with open(path, 'rb') as f:
        magic = f.read(2)
    if magic == b'\x1f\x8b':
        return gzip.open(path, 'rb')
    print(f'Note: {path.name} is plain text despite .gz suffix', flush=True)
    return open(path, 'rb')

def read_table(path):
    with open_by_magic(path) as f:
        return pd.read_csv(f, sep='\t', header=None)

check_inputs()
print(f'Loading {SAMPLE} from {COUNT}', flush=True)
with open_by_magic(COUNT / 'matrix.mtx.gz') as f:
    matrix = mmread(f)
if not sparse.issparse(matrix):
    matrix = sparse.coo_matrix(matrix)
barcodes = read_table(COUNT / 'barcodes.tsv.gz').iloc[:, 0].astype(str)
features = read_table(COUNT / 'features.tsv.gz')
if matrix.shape != (len(features), len(barcodes)):
    raise ValueError(f'Matrix {matrix.shape} does not match features ({len(features)}) and barcodes ({len(barcodes)})')
if features.shape[1] >= 3:
    keep = features.iloc[:, 2].astype(str).eq('Gene Expression').to_numpy()
else:
    keep = np.ones(len(features), dtype=bool)
if not keep.any():
    raise ValueError('No Gene Expression features found')
X = matrix.tocsr()[keep, :].T.tocsr()
selected_features = features.loc[keep].reset_index(drop=True)
adata = AnnData(X=X)
adata.obs_names = pd.Index(barcodes.to_numpy(), dtype=str)
adata.var_names = pd.Index(selected_features.iloc[:, 1].astype(str).to_numpy() if selected_features.shape[1] > 1 else selected_features.iloc[:, 0].astype(str).to_numpy())
adata.var['gene_ids'] = selected_features.iloc[:, 0].astype(str).to_numpy()
adata.var_names_make_unique()

positions = pd.read_csv(SPATIAL/'tissue_positions_list.csv', header=None,
    names=['barcode','in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres']).set_index('barcode')
missing = adata.obs_names.difference(positions.index)
if len(missing):
    raise ValueError(f'{len(missing)} expression barcodes have no spatial coordinates')
positions = positions.loc[adata.obs_names]
for col in ['in_tissue','array_row','array_col']:
    adata.obs[col] = positions[col].astype(int).to_numpy()
adata.obsm['spatial'] = positions[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy(dtype=float)
adata = adata[adata.obs['in_tissue'] == 1].copy()
adata.obs['sample_id'] = SAMPLE
adata.write_h5ad(RESULTS/f'{SAMPLE}_raw_spatial.h5ad', compression='gzip')
plot_spatial(adata, title=f'{SAMPLE}: H&E and Visium spots', outfile=FIG/f'S1_{SAMPLE}_tissue_spots.png')
print(f'S1 saved; {adata.n_obs} tissue spots × {adata.n_vars} genes')
