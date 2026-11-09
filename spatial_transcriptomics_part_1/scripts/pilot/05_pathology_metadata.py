import pandas as pd
from common import *

adata=load_processed()
report=RESULTS/f'{SAMPLE}_metadata_inspection.txt'
if not META.exists():
    report.write_text(f'Metadata file not found: {META}\nS5 skipped; no pathology comparison made.\n')
    print(report.read_text()); raise SystemExit(0)
df=pd.read_csv(META)
lines=[f'File: {META}',f'Shape: {df.shape}',f'Columns: {list(df.columns)}','',df.head(5).to_string(index=False)]
barcode_candidates=[c for c in df.columns if 'barcode' in c.lower() or c.lower().startswith('unnamed')]
if not barcode_candidates:
    lines.append('\nNo recognizable barcode column. Cannot align pathology labels; S5 skipped.')
else:
    col=barcode_candidates[0]
    df[col]=df[col].astype(str)
    common=set(adata.obs_names)&set(df[col])
    lines.append(f'\nBarcode column: {col}; matching spots: {len(common)} / {adata.n_obs}')
    # Only descriptive inspection: never assume an arbitrary column is pathology.
    lines.append('Review column names and values to identify a validated pathology annotation before S5 comparison.')
report.write_text('\n'.join(lines)+'\n')
print('\n'.join(lines[:5]))
print('S5 metadata inspection saved; pathology comparison intentionally deferred until label schema is verified.')
