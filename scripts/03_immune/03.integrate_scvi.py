# Purpose: integrate scvi.
import os
from pathlib import Path
ANALYSIS_DIR = Path(os.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input"))

import tempfile
import scanpy as sc
import scvi
import torch
import matplotlib.pyplot as plt
from rich import print
from sklearn.ensemble import RandomForestClassifier
from sklearn import datasets

scvi.settings.seed = 0
print(torch.cuda.is_available())
print(torch.__version__)
print("Last run with scvi-tools version:", scvi.__version__)

adata = sc.read_h5ad(str(ANALYSIS_DIR / "immune/input/immune.h5ad"))

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata)
adata.layers['cpm']=adata.X.copy()
sc.pp.log1p(adata)
adata.layers['data']=adata.X.copy()
adata.raw = adata

sc.pp.highly_variable_genes(adata, flavor="seurat_v3",n_top_genes=2000,layer="counts",batch_key="orig.ident",)
adata = adata[:, adata.var.highly_variable]

adata1 = adata.copy()
scvi.model.SCVI.setup_anndata(adata1, layer="counts", batch_key="sample", categorical_covariate_keys=["tissue","stage"])
model = scvi.model.SCVI(adata1, n_layers=2, n_latent=30, gene_likelihood="nb")
model.train()

model.save("immune_scvi_model_hvg.pkl",overwrite=True, save_anndata=True)

##########cluster
SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=0.3)

adata.write_h5ad('immune_scvi_hvg.h5ad')
