from pathlib import Path
import os, json
import numpy as np
import pandas as pd
import scanpy as sc
from PIL import Image

ROOT = Path(os.environ.get('WU_SPATIAL_ROOT', str(Path.home() / 'Downloads' / 'wu2021_spatial'))).expanduser()
SAMPLE = os.environ.get('WU_SPATIAL_SAMPLE', 'CID44971')
DATA = ROOT / 'data'
FIG = ROOT / 'figures'
RESULTS = ROOT / 'results'
for d in (FIG, RESULTS): d.mkdir(parents=True, exist_ok=True)
COUNT = DATA / 'filtered_count_matrices' / f'{SAMPLE}_filtered_count_matrix'
SPATIAL = DATA / 'spatial' / f'{SAMPLE}_spatial'
META = DATA / 'metadata' / f'{SAMPLE}_metadata.csv'

def check_inputs():
    required = [COUNT/'matrix.mtx.gz', COUNT/'barcodes.tsv.gz', COUNT/'features.tsv.gz', SPATIAL/'tissue_positions_list.csv', SPATIAL/'tissue_hires_image.png', SPATIAL/'scalefactors_json.json']
    missing = [str(p) for p in required if not p.exists()]
    if missing: raise FileNotFoundError('Missing required files:\n' + '\n'.join(missing))

def image_and_coords(adata):
    with (SPATIAL/'scalefactors_json.json').open() as f: sf = json.load(f)
    im = np.asarray(Image.open(SPATIAL/'tissue_hires_image.png').convert('RGB'))
    scale = float(sf['tissue_hires_scalef'])
    return im, adata.obsm['spatial'] * scale, sf

def plot_spatial(adata, values=None, title='', outfile=None, cmap='viridis', vmin=None, vmax=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Circle
    im, coords, sf = image_and_coords(adata)
    radius = float(sf['spot_diameter_fullres']) * float(sf['tissue_hires_scalef']) / 2
    fig, ax = plt.subplots(figsize=(9,9))
    ax.imshow(im)
    circles = [Circle((float(x),float(y)), radius=radius) for x,y in coords]
    coll = PatchCollection(circles, linewidth=0.25, alpha=0.85)
    if values is None:
        coll.set_facecolor('none'); coll.set_edgecolor('#00FFFF'); coll.set_linewidth(0.55)
    else:
        values = np.asarray(values, dtype=float)
        coll.set_array(values); coll.set_cmap(cmap); coll.set_edgecolor('none')
        if vmin is not None or vmax is not None: coll.set_clim(vmin=vmin, vmax=vmax)
    ax.add_collection(coll)
    if values is not None: fig.colorbar(coll, ax=ax, fraction=.035, pad=.015)
    ax.set_xlim(0, im.shape[1]); ax.set_ylim(im.shape[0], 0)
    ax.set_aspect('equal'); ax.axis('off'); ax.set_title(title)
    fig.tight_layout()
    if outfile: fig.savefig(outfile, dpi=180, bbox_inches='tight')
    plt.close(fig)

def load_raw():
    p = RESULTS/f'{SAMPLE}_raw_spatial.h5ad'
    if not p.exists(): raise FileNotFoundError(f'Run 01 first: {p}')
    return sc.read_h5ad(p)

def load_processed():
    p = RESULTS/f'{SAMPLE}_processed_spatial.h5ad'
    if not p.exists(): raise FileNotFoundError(f'Run 02 first: {p}')
    return sc.read_h5ad(p)

def gene_vector(adata, gene):
    if gene not in adata.var_names: return None
    col = adata[:,gene].X
    return col.toarray().ravel() if hasattr(col,'toarray') else np.asarray(col).ravel()
