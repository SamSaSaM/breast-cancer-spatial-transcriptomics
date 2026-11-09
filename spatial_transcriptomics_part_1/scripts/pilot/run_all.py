from pathlib import Path
import subprocess, sys
base=Path(__file__).resolve().parent
scripts=['01_load_spatial.py','02_spatial_qc.py','03_spatial_markers.py','04_spatial_clustering.py','05_pathology_metadata.py','06_macrophage_signatures.py']
for s in scripts:
    print(f'\n===== Running {s} =====',flush=True)
    result=subprocess.run([sys.executable,str(base/'scripts'/s)],cwd=base)
    if result.returncode:
        print(f'FAILED: {s} (exit code {result.returncode}). Fix this before proceeding.')
        sys.exit(result.returncode)
print('\nDONE! Figures and results are in ~/Downloads/wu2021_spatial/{figures,results}/')
