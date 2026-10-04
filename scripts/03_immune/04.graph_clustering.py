# Purpose: graph clustering.
import os
from pathlib import Path

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


adata = sc.read_h5ad('immune_scvi_hvg.h5ad')

sc.tl.leiden(adata, resolution=0.8)
sc.pl.embedding(adata, color="leiden", frameon=False, basis='X_umap', size=6)
fig = plt.gcf()
fig.set_size_inches(8, 8)
fig.savefig("immune_scvi_leiden_reso1.pdf", bbox_inches='tight')

sc.tl.rank_genes_groups(adata,  groupby="leiden", use_raw=False, )
sc.tl.dendrogram(adata, groupby="leiden", use_rep='X_scVI')

adata.write_h5ad('immune_cluster.h5ad')

df = sc.get.rank_genes_groups_df(adata, group=None,)
df.to_csv('immune_clusters_genes.txt', index=0, sep='\t')

sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden", standard_scale="var", n_genes=5, )
fig = plt.gcf()
fig.savefig("immune_fea_dot.pdf", bbox_inches='tight')


sc.pl.embedding(adata, color="leiden", frameon=False, basis='X_umap', size=6, legend_loc='on data')
fig = plt.gcf()
fig.set_size_inches(8, 8)
fig.savefig("immune_scvi_leiden_reso1_ondata.pdf", bbox_inches='tight')
