# Purpose: scvi graph clustering.
import os
from pathlib import Path

import scanpy as sc
import scvi
import torch
from rich import print
scvi.settings.seed = 0
print(torch.cuda.is_available())
print(torch.__version__)
print("Last run with scvi-tools version:", scvi.__version__)


tissues = ["all"]
# Loop over each tissue
for tissue in tissues:
    # Load the model and data for the specific tissue
    model = scvi.model.SCVI.load(f"{tissue}/01_scVI/{tissue}_scvi_model_hvg.pkl")
    adata = sc.read_h5ad(f'{tissue}/01_scVI/{tissue}_scvi_model_hvg.pkl/adata.h5ad')
    
    # Preprocessing steps
    sc.pp.scale(adata)
    adata.layers['scaled'] = adata.X.copy()
    sc.pp.pca(adata)
    
    # SCVI latent representation
    SCVI_LATENT_KEY = "X_scVI"
    adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()
    
    # Neighbors and clustering
    sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
    sc.tl.leiden(adata)
    
    # MDE and UMAP
    SCVI_MDE_KEY = "X_scVI_MDE"
    adata.obsm[SCVI_MDE_KEY] = scvi.model.utils.mde(adata.obsm[SCVI_LATENT_KEY])
    sc.tl.umap(adata, init_pos=adata.obsm[SCVI_MDE_KEY])
    
    # Save the processed data
    adata.write_h5ad(f'{tissue}_scvi_hvg.h5ad')

# cluster


tissues = ["all"]
resolutions = [0.2, 0.4, 0.6, 0.8, 1.0]

# Perform clustering for each tissue
for tissue in tissues:
    adata = sc.read_h5ad(f"{tissue}_scvi_hvg.h5ad")
    
    # Perform Leiden clustering with different resolutions
    for res in resolutions:
        cluster_key = f"leiden_res{str(res).replace('.', '_')}"  # Create a dynamic key
        sc.tl.leiden(adata, key_added=cluster_key, resolution=res)
    
    # Save the clustering results
    adata.write_h5ad(f'{tissue}_scvi_cluster.h5ad')
