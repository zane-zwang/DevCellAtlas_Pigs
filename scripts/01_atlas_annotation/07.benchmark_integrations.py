# Purpose: benchmark integrations.
import os
from pathlib import Path
ANALYSIS_DIR = Path(os.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input"))

import tempfile
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
from rich import print
from scib_metrics.benchmark import Benchmarker
import pandas as pd
import numpy as np

# 1. import AnnData

print("Loading integrated AnnData ...")

adata = sc.read_h5ad(
    str(ANALYSIS_DIR / "atlas/scib/TissuesAll_batch_all.h5ad")

)

print(
    "Loading scVI embedding ..."
)

scvi_embedding = pd.read_csv(
    str(ANALYSIS_DIR / "atlas/scib/all_X_scVI.tsv"),

    sep="\t",
    index_col=0
)

# Align cells shared across integration embeddings.

print("\nadata shape:")
print(adata.shape)

print(
    "scVI embedding shape:"
)
print(
    scvi_embedding.shape
)

common_cells = (
    adata.obs_names
    .intersection(
        scvi_embedding.index
    )
)

print(
    "Common cells:",
    len(common_cells)
)

print(
    "adata cells:",
    adata.n_obs
)

if len(common_cells) == 0:
    raise ValueError(
        "No common cells found between AnnData and scVI embedding."
    )

if len(common_cells) < adata.n_obs:

    print(
        "Warning: not all cells are common. "
        "Only common cells retained."
    )

# 3. subset main object

adata = adata[
    common_cells
].copy()

# ensure same order

scvi_embedding = (
    scvi_embedding
    .loc[
        adata.obs_names
    ]
)

assert np.all(
    adata.obs_names ==
    scvi_embedding.index
)

# 4. transfer embeddings from adata

print(
    "\nEmbeddings in adata:"
)

for key in adata.obsm.keys():

    print(key)

rename_dict_adata1 = {
    "scaled|original|X_pca":"X_pca",
    "X_combat":"X_ComBat",
    "X_ComBat":"X_ComBat",
    "X_cca":"X_CCA",
    "X_CCA":"X_CCA",
    "X_harmony":"X_Harmony",
    "X_pca_harmony":"X_Harmony",
    "X_Harmony":"X_Harmony",
    "X_scanorama":"X_scanorama",
}

for old_key, new_key in rename_dict_adata1.items():

    if old_key in adata.obsm.keys():
        adata.obsm[new_key] = (
            np.asarray(
                adata.obsm[old_key]
            )
            .copy()
        )
        print(
            f"Kept: {old_key} -> {new_key}"
        )

# 5. transfer scVI embedding

adata.obsm["X_scVI"] = (
    scvi_embedding
    .values
    .astype(
        np.float32
    )
)

print(
    "Transferred: scVI embedding -> X_scVI"
)

# Validate integration embeddings.

print(
    "\nFinal embeddings in unified AnnData:"
)

for key in adata.obsm.keys():
    print(
        key,
        adata.obsm[key].shape
    )

# 7. save unified AnnData

adata.write_h5ad(
    "adata_all_batch_methods.h5ad",
    compression="gzip"
)

print(
    "\nSaved:"
    "adata_all_batch_methods.h5ad"
)

# 8. run scIB benchmark

embedding_keys = [
    "X_pca",
    "X_scVI",
    "X_ComBat",
    "X_CCA",
    "X_Harmony",
    "X_scanorama"
]

# Use only embeddings present in the combined AnnData object.

embedding_keys = [

    x for x in embedding_keys

    if x in adata.obsm.keys()

]

print(
    "\nBenchmark embeddings:"
)

print(
    embedding_keys
)

bm = Benchmarker(
    adata,
    batch_key="sample",
    label_key="celltype",
    embedding_obsm_keys=
        embedding_keys,

    n_jobs=5
)

print(
    "\nRunning scIB benchmark..."
)

bm.benchmark()

# 9. plot result

bm.plot_results_table(
    min_max_scale=False
)

fig = plt.gcf()
fig.savefig(
    "scib_result.pdf",
    bbox_inches="tight",
    dpi=300
)

# 10. save benchmark object

import pickle

with open(
    "bm_results.pkl",
    "wb"
) as f:
    pickle.dump(
        bm,
        f
    )

print(
    "\nALL DONE"
)
