# Purpose: export scvi embedding.
import os
from pathlib import Path
ANALYSIS_DIR = Path(os.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input"))

import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad(
    str(ANALYSIS_DIR / "atlas/01_scvi/all/all_scvi_hvg.h5ad")
)

print(adata.obsm.keys())

pd.DataFrame(
    adata.obsm["X_scVI"],
    index=adata.obs_names
).to_csv(
    "all_X_scVI.tsv",
    sep="\t"
)

adata.obs.to_csv(
    "all_metadata.tsv",
    sep="\t"
)
