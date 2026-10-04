# Purpose: Supplementary camex embedding.
import os
from pathlib import Path

import scanpy as sc
import matplotlib.pyplot as plt

adata_CAMEX = sc.read_h5ad('adata_CAMEX.h5ad')
sc.pp.neighbors(adata_CAMEX, use_rep='X_CAMEX_Integration')
sc.tl.umap(adata_CAMEX)

sc.pl.umap(adata_CAMEX, color=['batch'], wspace=0.4)
fig = plt.gcf()
fig.set_size_inches(6, 6)
fig.savefig('3_species_CAMEX_batch.png', bbox_inches='tight')
fig.savefig('3_species_CAMEX_batch.pdf', bbox_inches='tight')
