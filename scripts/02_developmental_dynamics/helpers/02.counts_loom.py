"""Create a pySCENIC loom from the raw count matrix and cisTarget gene universe."""

import os
from pathlib import Path

import loompy as lp
import numpy as np
import pandas as pd
import scanpy as sc


COUNTS_CSV = Path(os.getenv("SCENIC_COUNTS_CSV", "RNAcounts.csv"))
RANKING_FEATHER = Path(
    os.getenv("SCENIC_RANKING_FEATHER", "path/to/reference/10KbUP_10KbDOWN.regions_vs_motifs.rankings.feather")
)
OUTPUT_LOOM = Path(os.getenv("SCENIC_COUNTS_LOOM", "RNAcounts.loom"))

raw_mat = sc.read_csv(str(COUNTS_CSV))
ranking_feather = pd.read_feather(RANKING_FEATHER)
exp_mat = raw_mat[:, raw_mat.var_names.isin(ranking_feather.columns)]

row_attrs = {"Gene": np.array(exp_mat.var_names)}
col_attrs = {
    "CellID": np.array(exp_mat.obs_names),
    "nGene": np.array(np.sum(exp_mat.X.transpose() > 0, axis=0)).flatten(),
    "nUMI": np.array(np.sum(exp_mat.X.transpose(), axis=0)).flatten(),
}
lp.create(str(OUTPUT_LOOM), exp_mat.X.transpose(), row_attrs, col_attrs)
