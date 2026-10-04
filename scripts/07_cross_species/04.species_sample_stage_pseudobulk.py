# Purpose: species sample stage pseudobulk.
import os
from pathlib import Path
ANALYSIS_DIR = Path(os.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input"))

import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc
import anndata as ad


def pseudobulk_from_anndata(
    adata,
    groupby_cols=("species", "sample", "stage"),
    layer="counts",
    use_raw=False,
    min_count=1,
    min_groups=2,
    min_cells_per_group=1,
    sep="//",
    verbose=True
):
    """
    Construct pseudobulk matrix from AnnData by summing counts within groups.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object.
    groupby_cols : tuple/list
        obs columns used to define pseudobulk groups.
    layer : str or None
        Layer name containing raw counts. If None, use adata.X.
    use_raw : bool
        Whether to use adata.raw.X.
    min_count : int or float
        Expression threshold for filtering genes after pseudobulk construction.
    min_groups : int
        Keep genes expressed (> min_count) in at least this many pseudobulk samples.
    min_cells_per_group : int
        Minimum number of cells required for a group to be kept.
    sep : str
        Separator for group labels.
    verbose : bool
        Whether to print progress.

    Returns
    -------
    expr_df : pd.DataFrame
        Gene x pseudobulk_sample matrix.
    group_info : pd.DataFrame
        Metadata for pseudobulk samples.
    """

    # 1. basic checks
    obs = adata.obs.copy()

    missing_cols = [col for col in groupby_cols if col not in obs.columns]
    if len(missing_cols) > 0:
        raise ValueError(f"Columns not found in adata.obs: {missing_cols}")

    # remove NA groups
    keep_cells = obs[list(groupby_cols)].notna().all(axis=1)
    n_removed_na = int((~keep_cells).sum())
    if n_removed_na > 0 and verbose:
        print(f"[Info] Removing {n_removed_na} cells with NA in grouping columns.")
    adata = adata[keep_cells].copy()
    obs = adata.obs.copy()

    # 2. choose expression matrix
    if use_raw:
        if adata.raw is None:
            raise ValueError("use_raw=True but adata.raw is None.")
        X = adata.raw.X
        var_names = adata.raw.var_names
    elif layer is not None:
        if layer not in adata.layers:
            raise ValueError(f"Layer '{layer}' not found in adata.layers.")
        X = adata.layers[layer]
        var_names = adata.var_names
    else:
        X = adata.X
        var_names = adata.var_names

    if sparse.issparse(X):
        X = X.tocsr()
    else:
        X = np.asarray(X)

    # 3. build group labels
    group_labels = obs[list(groupby_cols)].astype(str).agg(sep.join, axis=1)
    obs["_pseudobulk_group"] = group_labels.values
    adata.obs["_pseudobulk_group"] = group_labels.values

    # count cells per group
    group_sizes = obs["_pseudobulk_group"].value_counts(sort=False)

    # keep groups with enough cells
    keep_groups = group_sizes[group_sizes >= min_cells_per_group].index.tolist()
    if len(keep_groups) == 0:
        raise ValueError("No groups left after applying min_cells_per_group filter.")

    dropped_groups = group_sizes[group_sizes < min_cells_per_group]
    if len(dropped_groups) > 0 and verbose:
        print(f"[Info] Dropping {len(dropped_groups)} groups with < {min_cells_per_group} cells.")

    keep_group_mask = obs["_pseudobulk_group"].isin(keep_groups).values
    adata = adata[keep_group_mask].copy()
    obs = adata.obs.copy()
    if sparse.issparse(X):
        X = X[keep_group_mask]
    else:
        X = X[keep_group_mask, :]

    unique_groups = obs["_pseudobulk_group"].unique().tolist()

    # 4. aggregate counts
    pseudobulk_list = []
    group_meta = []

    group_array = obs["_pseudobulk_group"].values

    for g in unique_groups:
        idx = np.where(group_array == g)[0]

        if sparse.issparse(X):
            vec = np.asarray(X[idx].sum(axis=0)).ravel()
        else:
            vec = X[idx].sum(axis=0)

        pseudobulk_list.append(vec)

        meta_row = {
            "pseudobulk_id": g,
            "n_cells": len(idx),
            "dataset": obs.iloc[idx[0]]["dataset"]
        }
        for col in groupby_cols:
            meta_row[col] = obs.iloc[idx[0]][col]
        group_meta.append(meta_row)

    expr_mat = np.column_stack(pseudobulk_list)
    expr_df = pd.DataFrame(expr_mat, index=var_names, columns=unique_groups)

    # 5. gene filtering
    keep_genes = (expr_df > min_count).sum(axis=1) >= min_groups
    n_before = expr_df.shape[0]
    expr_df = expr_df.loc[keep_genes]
    n_after = expr_df.shape[0]

    if verbose:
        print(f"[Info] Pseudobulk matrix shape before filtering: {n_before} genes x {len(unique_groups)} samples")
        print(f"[Info] Pseudobulk matrix shape after filtering : {n_after} genes x {len(unique_groups)} samples")

    group_info = pd.DataFrame(group_meta)

    return expr_df, group_info


def make_pseudobulk_anndata(expr_df, group_info):
    """
    Convert pseudobulk matrix + metadata to AnnData.
    Output AnnData is pseudobulk_sample x gene.
    """
    pb_adata = ad.AnnData(
        X=expr_df.T.values,
        obs=group_info.set_index("pseudobulk_id"),
        var=pd.DataFrame(index=expr_df.index)
    )
    return pb_adata


if __name__ == "__main__":
    # -------- input --------
    input_file = str(ANALYSIS_DIR / "cross_species/camex/adata_CAMEX_umap.h5ad")
    output_prefix = "pseudobulk_species_sample_stage"

    adata = sc.read_h5ad(input_file)

    expr_df, group_info = pseudobulk_from_anndata(
        adata,
        groupby_cols=("species", "sample", "stage"),
        layer="counts",              # layer = None
        use_raw=False,
        min_count=1,
        min_groups=2,
        min_cells_per_group=10,      # 1/5/10
        sep="//",
        verbose=True
    )

    expr_df.to_csv(f"{output_prefix}.csv")
    group_info.to_csv(f"{output_prefix}_meta.csv", index=False)

    pb_adata = make_pseudobulk_anndata(expr_df, group_info)
    pb_adata.write_h5ad(f"{output_prefix}.h5ad")

    print("[Done] Files written:")
    print(f"  - {output_prefix}.csv")
    print(f"  - {output_prefix}_meta.csv")
    print(f"  - {output_prefix}.h5ad")
