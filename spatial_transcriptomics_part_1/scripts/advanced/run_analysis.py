#!/usr/bin/env python3
"""Wu etal 2021 Visium exploratory analysis. See README for scientific limitations."""
import argparse, csv, gzip, json, logging, re, shutil, tempfile, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
import scanpy as sc

SAMPLES = ['CID44971','CID4535','1160920F','CID4290','1142243F','CID4465']
SIGNATURES = {
 'epithelial': ['EPCAM','KRT8','KRT18','KRT19','KRT7'],
 'stromal': ['COL1A1','COL1A2','DCN','LUM','COL3A1'],
 'immune': ['PTPRC','LST1','TYROBP','CD3D','MS4A1'],
 'macrophage': ['C1QA','C1QB','C1QC','CD68','CSF1R','LST1','TYROBP'],
 'antigen_presentation': ['HLA-DRA','HLA-DRB1','HLA-DPA1','CD74'],
}

def setup_logger(logfile):
 logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s',handlers=[logging.StreamHandler(),logging.FileHandler(logfile,encoding='utf-8')],force=True)

def open_maybe_gzip(p):
 with open(p,'rb') as f: magic=f.read(2)
 return gzip.open(p,'rb') if magic == b'\x1f\x8b' else open(p,'rb')

def load_matrix(directory):
 # Some deposited *.gz files contain uncompressed data; inspect magic bytes.
 with open_maybe_gzip(directory/'matrix.mtx.gz') as f: X=mmread(f)
 X=X.T.tocsr().astype(np.float32) if sp.issparse(X) else sp.csr_matrix(X.T,dtype=np.float32)
 with open_maybe_gzip(directory/'barcodes.tsv.gz') as f:
  barcodes=pd.read_csv(f,sep='\t',header=None,dtype=str).iloc[:,0].tolist()
 with open_maybe_gzip(directory/'features.tsv.gz') as f:
  features=pd.read_csv(f,sep='\t',header=None,dtype=str)
 if X.shape != (len(barcodes),len(features)):
  raise ValueError(f'Matrix shape {X.shape} inconsistent with {len(barcodes)} barcodes and {len(features)} features')
 if features.shape[1]>=3:
  keep=features.iloc[:,2].eq('Gene Expression').to_numpy()
  if keep.any(): X=X[:,keep];features=features.loc[keep].reset_index(drop=True)
 import anndata as ad
 a=ad.AnnData(X=X)
 a.obs_names=pd.Index(barcodes)
 # Some GEO deposits have only one feature column (gene symbols).
 names = features.iloc[:, 1] if features.shape[1] >= 2 else features.iloc[:, 0]
 a.var_names=pd.Index(names.fillna(features.iloc[:,0]).astype(str).to_numpy(), name=None)
 a.var_names_make_unique()
 a.var['gene_id']=features.iloc[:,0].astype(str).values
 return a

def load_positions(path):
 p=path/'tissue_positions_list.csv'
 if not p.exists(): p=path/'tissue_positions.csv'
 if not p.exists(): raise FileNotFoundError('No tissue positions file in '+str(path))
 first=pd.read_csv(p,nrows=1,header=None).iloc[0,0]
 if str(first).lower()=='barcode':
  pos=pd.read_csv(p)
  pos=pos.rename(columns={'pxl_row_in_fullres':'pxl_row_in_fullres','pxl_col_in_fullres':'pxl_col_in_fullres'})
 else:
  pos=pd.read_csv(p,header=None,names=['barcode','in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres'])
 return pos.set_index('barcode')

def plot_spatial(a,spatial_dir,values,title,out,categorical=False):
 imgpath=spatial_dir/'tissue_hires_image.png'
 scale=1.0
 if (spatial_dir/'scalefactors_json.json').exists():
  with open(spatial_dir/'scalefactors_json.json') as f: scale=float(json.load(f).get('tissue_hires_scalef',1.0))
 fig,ax=plt.subplots(figsize=(8,8))
 if imgpath.exists():
  with Image.open(imgpath) as im: img=np.asarray(im.convert('RGB'))
  ax.imshow(img)
 xy=a.obsm['spatial']*scale
 if categorical:
  cats=pd.Categorical(values)
  palette=plt.get_cmap('tab20')
  for i,c in enumerate(cats.categories):
   mask=np.asarray(cats==c)
   ax.scatter(xy[mask,0],xy[mask,1],s=9,alpha=.9,label=str(c),color=palette(i%20),linewidths=0)
  ax.legend(loc='upper left',bbox_to_anchor=(1.0,1.0),fontsize=8,markerscale=2)
 else:
  pts=ax.scatter(xy[:,0],xy[:,1],c=np.asarray(values,dtype=float),s=10,cmap='viridis',linewidths=0)
  fig.colorbar(pts,ax=ax,fraction=.035,pad=.02)
 if imgpath.exists(): ax.set_xlim(0,img.shape[1]);ax.set_ylim(img.shape[0],0)
 else: ax.invert_yaxis()
 ax.set_aspect('equal');ax.axis('off');ax.set_title(title)
 fig.savefig(out,dpi=180,bbox_inches='tight');plt.close(fig)

def score_signature(a,genes):
 present=[g for g in genes if g in a.var_names]
 if len(present)<2: return None,present
 sub=a[:,present].X
 v=np.asarray(sub.mean(axis=1)).ravel()
 return v,present

def detect_pathology(a,metadata_path,outdir):
 info={'metadata_file':str(metadata_path),'metadata_found':metadata_path.exists(),'matched_spots':0,'pathology_column':None,'status':'unavailable'}
 if not metadata_path.exists(): return info
 try: df=pd.read_csv(metadata_path,low_memory=False)
 except Exception as e: info['status']=f'read error: {e}';return info
 df.head(15).to_csv(outdir/'metadata_preview.csv',index=False)
 info['columns']=list(map(str,df.columns))
 best=(0,None,None)
 for col in df.columns:
  series=df[col].astype(str)
  for transform in [lambda x:x,lambda x:x.str.replace(r'-1$','',regex=True)]:
   val=transform(series)
   obs=pd.Index(a.obs_names)
   obs2=transform(pd.Series(obs)).values
   overlap=len(set(val)&set(obs2))
   if overlap>best[0]:best=(overlap,col,transform)
 overlap,col,transform=best
 info['matched_spots']=overlap
 if col is None or overlap<max(10,int(.5*a.n_obs)):
  info['status']='No confidently matching barcode column; metadata comparison skipped'
  return info
 tmp=df.copy();tmp['_join_key']=transform(tmp[col].astype(str));tmp=tmp.drop_duplicates('_join_key').set_index('_join_key')
 obskey=transform(pd.Series(a.obs_names)).values
 candidate=[c for c in df.columns if c!=col and re.search(r'patholog|histolog|annotat|region|morpholog|tissue.type',str(c),re.I)]
 info['candidate_annotation_columns']=list(map(str,candidate))
 for c in candidate:
  v=tmp[c].reindex(obskey)
  if v.notna().sum()<max(10,.5*a.n_obs) or v.nunique(dropna=True)<2 or v.nunique(dropna=True)>30: continue
  info['pathology_column']=str(c)
  a.obs['pathology_annotation']=pd.Categorical(v.astype('string').fillna('Unknown').to_numpy(dtype=str))
  pd.crosstab(a.obs['leiden'].to_numpy(),a.obs['pathology_annotation'].to_numpy(),normalize='index').to_csv(outdir/'leiden_vs_pathology_fraction.csv')
  info['status']='Annotation comparison produced (candidate column; verify biological meaning)'
  return info
 info['status']='Barcode matched, but no clearly named pathology column found; skipped'
 return info

def neighbor_analysis(a,outdir,sample,seed=42,nperm=199):
 xy=np.asarray(a.obsm['spatial'],float)
 tree=cKDTree(xy)
 d,idx=tree.query(xy,k=2)
 nearest=d[:,1]
 median=float(np.median(nearest[nearest>0])) if np.any(nearest>0) else 0
 radius=median*1.35
 pairs=np.array(list(tree.query_pairs(radius)),dtype=int).reshape(-1,2)
 if len(pairs)<10: return {'status':'Insufficient nearby spot pairs','pairs':len(pairs)}
 e=np.asarray(a.obs['epithelial_score'],float);m=np.asarray(a.obs['macrophage_score'],float)
 # Two separate readouts: global spot-level correlation and neighboring-pair cross-signal.
 rho,p=spearmanr(e,m)
 observed=float(np.mean((e[pairs[:,0]]*m[pairs[:,1]]+e[pairs[:,1]]*m[pairs[:,0]])/2))
 rng=np.random.default_rng(seed); null=[]
 for _ in range(nperm):
  perm=rng.permutation(m)
  null.append(float(np.mean((e[pairs[:,0]]*perm[pairs[:,1]]+e[pairs[:,1]]*perm[pairs[:,0]])/2)))
 null=np.asarray(null)
 pperm=(1+np.count_nonzero(np.abs(null-null.mean())>=abs(observed-null.mean())))/(len(null)+1)
 # Descriptive top-quartile neighborhood enrichment, not deconvolution.
 etop=e>=np.quantile(e,.75);mtop=m>=np.quantile(m,.75)
 adjacent=(etop[pairs[:,0]]&mtop[pairs[:,1]])|(etop[pairs[:,1]]&mtop[pairs[:,0]])
 results={'status':'exploratory','n_spots':len(e),'neighbor_pairs':len(pairs),'neighbor_radius_fullres_pixels':radius,
 'spot_score_spearman_rho':float(rho),'spot_score_spearman_p':float(p),
 'cross_neighbor_score_mean':observed,'permutation_null_mean':float(null.mean()),'permutation_two_sided_p':float(pperm),
 'fraction_neighbor_pairs_high_epithelial_high_macrophage':float(adjacent.mean()),
 'interpretation_warning':'Scores are not cell fractions; spatial proximity is spot-level; permutation ignores tissue compartments and spatial autocorrelation.'}
 pd.DataFrame([results]).to_csv(outdir/'spatial_neighborhood_statistics.csv',index=False)
 return results

def analyze(sample,data,output,seed,nperm):
 sd=output/sample;sd.mkdir(parents=True,exist_ok=True)
 count=data/'filtered_count_matrices'/f'{sample}_filtered_count_matrix'
 spatial=data/'spatial'/f'{sample}_spatial'
 metadata=data/'metadata'/f'{sample}_metadata.csv'
 if not count.exists() or not spatial.exists(): raise FileNotFoundError(f'Missing count/spatial directory for {sample}')
 logging.info('%s: reading expression matrix',sample)
 a=load_matrix(count);pos=load_positions(spatial)
 matched=a.obs_names.intersection(pos.index)
 if len(matched)<.5*a.n_obs: raise ValueError(f'Only {len(matched)}/{a.n_obs} barcodes match spatial coordinates')
 a=a[matched].copy();pos=pos.loc[a.obs_names]
 a.obs['in_tissue']=pd.to_numeric(pos['in_tissue'],errors='coerce').fillna(0).astype(int).values
 a.obsm['spatial']=pos[['pxl_col_in_fullres','pxl_row_in_fullres']].astype(float).to_numpy()
 a=a[a.obs['in_tissue']==1].copy()
 if a.n_obs<30: raise ValueError(f'Too few tissue spots: {a.n_obs}')
 a.var['mt']=a.var_names.str.upper().str.startswith('MT-')
 sc.pp.calculate_qc_metrics(a,qc_vars=['mt'],percent_top=None,inplace=True)
 for field in ['total_counts','n_genes_by_counts','pct_counts_mt']:
  plot_spatial(a,spatial,a.obs[field],f'{sample} {field}',sd/f'QC_{field}.png')
 qc=a.obs[['total_counts','n_genes_by_counts','pct_counts_mt']].describe().to_csv(sd/'qc_summary.csv')
 # Retain tissue spots for spatial interpretation; remove genes never expressed.
 sc.pp.filter_genes(a,min_cells=3)
 sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
 # Save normalized data for reproducibility; no image embedded.
 a.write_h5ad(sd/'normalized_spatial.h5ad',compression='gzip')
 for name,genes in SIGNATURES.items():
  score,present=score_signature(a,genes)
  if score is None:
   logging.warning('%s: %s signature has <2 genes, skipping',sample,name);continue
  a.obs[f'{name}_score']=score
  plot_spatial(a,spatial,score,f'{sample}: {name} (mean log1p normalized expression)',sd/f'SIGNATURE_{name}.png')
  (sd/f'SIGNATURE_{name}_genes.txt').write_text('Present genes: '+', '.join(present)+'\nMissing genes: '+', '.join(set(genes)-set(present))+'\n')
 if 'macrophage_score' in a.obs and 'epithelial_score' in a.obs:
  neighbor=neighbor_analysis(a,sd,sample,seed,nperm)
 else: neighbor={'status':'required signatures unavailable'}
 # HVG clustering; fall back when variance selection fails.
 try: sc.pp.highly_variable_genes(a,n_top_genes=min(2000,a.n_vars),flavor='seurat')
 except Exception as exc:
  logging.warning('HVG fallback: %s',exc);a.var['highly_variable']=True
 sc.pp.pca(a,n_comps=min(30,a.n_obs-1,int(a.var['highly_variable'].sum())-1),use_highly_variable=True)
 sc.pp.neighbors(a,n_neighbors=min(15,a.n_obs-1))
 sc.tl.leiden(a,resolution=.5,key_added='leiden',flavor='igraph',n_iterations=2,directed=False)
 plot_spatial(a,spatial,a.obs['leiden'],f'{sample}: Leiden clusters',sd/'LEIDEN_spatial.png',categorical=True)
 a.obs['leiden'].value_counts().rename_axis('cluster').to_csv(sd/'leiden_cluster_sizes.csv')
 # DEGs are descriptive across spots; spatial autocorrelation invalidates naïve spot-wise inference.
 try:
  sc.tl.rank_genes_groups(a,'leiden',method='wilcoxon',use_raw=False,pts=False)
  deg=sc.get.rank_genes_groups_df(a,group=None)
  deg.to_csv(sd/'leiden_marker_genes_all.csv',index=False)
  deg.groupby('group',sort=False).head(20).to_csv(sd/'leiden_top20_markers.csv',index=False)
 except Exception as exc: logging.warning('%s: DEG failed: %s',sample,exc)
 pathinfo=detect_pathology(a,metadata,sd)
 if 'pathology_annotation' in a.obs:
  plot_spatial(a,spatial,a.obs['pathology_annotation'],f'{sample}: candidate pathology annotations',sd/'PATHOLOGY_spatial.png',categorical=True)
 a.obs[['leiden']+[c for c in a.obs if c.endswith('_score')]].to_csv(sd/'spot_clusters_and_scores.csv')
 summary={'sample':sample,'spots':a.n_obs,'genes_after_min_cells':a.n_vars,'leiden_clusters':int(a.obs['leiden'].nunique()),'metadata':pathinfo,'neighborhood':neighbor,
 'limitations':['Leiden clusters are transcriptional regions, not cell types.','Macrophage and epithelial scores are not cell abundances or proof of malignancy.','Spatial spot-level tests are exploratory due to autocorrelation.','Single sections do not support patient-level generalization.']}
 (sd/'analysis_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,default=str))
 logging.info('%s completed: %s tissue spots, %s clusters',sample,a.n_obs,summary['leiden_clusters'])
 return summary

def main():
 p=argparse.ArgumentParser()
 p.add_argument('--data',type=Path,default=Path.home()/'Downloads/wu2021_spatial/data')
 p.add_argument('--output',type=Path,default=Path.home()/'Downloads/wu2021_spatial/advanced_results')
 p.add_argument('--sample',choices=SAMPLES+['all'],default='all')
 p.add_argument('--permutations',type=int,default=199)
 p.add_argument('--seed',type=int,default=42)
 args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
 setup_logger(args.output/'run.log')
 selected=SAMPLES if args.sample=='all' else [args.sample]
 summaries=[];failures=[]
 for sample in selected:
  try: summaries.append(analyze(sample,args.data,args.output,args.seed,args.permutations))
  except Exception as e:
   logging.exception('%s failed',sample);failures.append({'sample':sample,'error':str(e)})
 if summaries:
  rows=[]
  for s in summaries:
   n=s['neighborhood'];rows.append({'sample':s['sample'],'spots':s['spots'],'leiden_clusters':s['leiden_clusters'],
    'macrophage_epithelial_spearman_rho':n.get('spot_score_spearman_rho'),
    'neighbor_permutation_p_exploratory':n.get('permutation_two_sided_p'),
    'pathology_status':s['metadata']['status']})
  pd.DataFrame(rows).to_csv(args.output/'all_samples_summary.csv',index=False)
 (args.output/'failures.json').write_text(json.dumps(failures,indent=2))
 with open(args.output/'RESULTS_README.md','w') as f:
  f.write('# Wu et al. (2021) spatial reanalysis: computed outputs\n\n')
  f.write('These are exploratory spot-level analyses, **not** validated cell-type deconvolution or exact figure reproduction.\n\n')
  for s in summaries:
   f.write(f"## {s['sample']}\n\n- Tissue spots: {s['spots']}\n- Leiden clusters: {s['leiden_clusters']}\n- Pathology: {s['metadata']['status']}\n")
   n=s['neighborhood']
   if 'spot_score_spearman_rho' in n:f.write(f"- Epithelial vs macrophage score Spearman rho: {n['spot_score_spearman_rho']:.3f}\n- Exploratory neighborhood permutation p: {n['permutation_two_sided_p']:.3f}\n")
   f.write('\n')
  if failures:f.write('## Failed samples\n'+json.dumps(failures,indent=2)+'\n')
  f.write('## Interpretation warnings\n\n- The macrophage signature is a proxy for expression, not a macrophage cell fraction.\n- Epithelial signature is not a malignant-cell classifier.\n- A positive spatial relationship can reflect tissue compartments or shared QC gradients.\n- Spot-wise Wilcoxon p-values are not independent biological replicates.\n- Cross-patient inference requires patient-level analysis and validation.\n- Pathology comparison is only performed when a plausible annotation column and barcode match exist; manually inspect the chosen column.\n')
 print('\nFinished. Results:',args.output,'\nSuccess:',len(summaries),'Failures:',len(failures))
 if failures: print('Check run.log and failures.json for errors.')

if __name__=='__main__':main()
