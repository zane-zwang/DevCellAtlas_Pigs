# Purpose: integrate scvi all.
import os
from pathlib import Path
ANALYSIS_DIR = Path(os.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input"))

import scanpy as sc
import scvi
import torch
import matplotlib.pyplot as plt
from rich import print


# settings

scvi.settings.seed = 0

print(torch.cuda.is_available())
print(torch.__version__)
print("scvi-tools version:", scvi.__version__)


# tissue list
tissues = [
    #"adipose",
    #"liver",
    #"muscle",
    #"heart",
    #"hypothalamus",
    #"cerebrum",
    #"duodenum",
    "all"
]


input_dir = str(ANALYSIS_DIR / "atlas")
output_dir = str(ANALYSIS_DIR / "atlas/01_scvi")

os.makedirs(output_dir, exist_ok=True)

# loop

for tissue in tissues:

    print("\n")
    print("="*50)
    print("Running:", tissue)
    print("="*50)


    adata = sc.read_h5ad(
        os.path.join(
            input_dir,
            f"{tissue}.h5ad"
        )
    )

    # counts
    adata.layers["counts"] = adata.X.copy()

    # normalization

    sc.pp.normalize_total(adata)

    adata.layers["cpm"] = adata.X.copy()

    sc.pp.log1p(adata)

    adata.layers["data"] = adata.X.copy()

    adata.raw = adata

    # HVG
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=3000,
        layer="counts",
        batch_key="sample",
    )

    adata = adata[
        :,
        adata.var.highly_variable
    ].copy()

    # scVI
    adata1 = adata.copy()

    scvi.model.SCVI.setup_anndata(
        adata1,
        layer="counts",
        batch_key="sample",
        categorical_covariate_keys=["tissue","stage"]
    )

    model = scvi.model.SCVI(
        adata1,
        n_layers=2,
        n_latent=30,
        gene_likelihood="nb"
    )

    model.train()

    # save model
    tissue_out = os.path.join(
        output_dir,
        tissue
    )

    os.makedirs(
        tissue_out,
        exist_ok=True
    )

    model.save(
        os.path.join(
            tissue_out,
            f"{tissue}_scvi_model"
        ),
        overwrite=True,
        save_anndata=True
    )

    # clustering
    adata.obsm["X_scVI"] = (
        model.get_latent_representation()
    )

    sc.pp.neighbors(
        adata,
        use_rep="X_scVI"
    )

    sc.tl.umap(
        adata,
        min_dist=0.3
    )

    # save h5ad
    adata.write_h5ad(
        os.path.join(
            tissue_out,
            f"{tissue}_scvi_hvg.h5ad"
        )
    )

    print(
        tissue,
        "finished"
    )

print("ALL DONE")
