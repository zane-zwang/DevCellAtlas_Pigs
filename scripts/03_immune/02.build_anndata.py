# Purpose: build anndata.
import os
from pathlib import Path

import pandas as pd
import anndata as ad
from scipy.io import mmread
from scipy import sparse

def build_anndata(prefix, outdir="."):
    """
    prefix: input file prefix, for example 'mouse_thymus'
    outdir: output directory
    """

    counts_file = f"{prefix}_counts.mtx"
    genes_file = f"{prefix}_genes.txt"
    barcodes_file = f"{prefix}_barcodes.txt"
    meta_file = f"{prefix}_meta.csv"

    print(f"\n=== Processing {prefix} ===")

    counts = mmread(counts_file)

    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)
    else:
        counts = counts.tocsr()

    print("Raw counts shape (genes x cells):", counts.shape)

    genes = pd.read_csv(genes_file, header=None, sep="\t")
    barcodes = pd.read_csv(barcodes_file, header=None, sep="\t")

    genes = genes.iloc[:, 0].astype(str).tolist()
    barcodes = barcodes.iloc[:, 0].astype(str).tolist()

    print("Number of genes:", len(genes))
    print("Number of barcodes:", len(barcodes))

    meta = pd.read_csv(meta_file, index_col=0)
    meta.index = meta.index.astype(str)

    print("Meta shape before alignment:", meta.shape)

    if counts.shape[0] != len(genes):
        raise ValueError(
            f"{prefix}: counts rows({counts.shape[0]}) != number of genes({len(genes)})"
        )
    if counts.shape[1] != len(barcodes):
        raise ValueError(
            f"{prefix}: counts columns({counts.shape[1]}) != number of barcodes({len(barcodes)})"
        )

    missing_barcodes = [bc for bc in barcodes if bc not in meta.index]
    extra_meta = [idx for idx in meta.index if idx not in barcodes]

    print("Barcodes missing in meta:", len(missing_barcodes))
    print("Extra cells in meta:", len(extra_meta))

    if len(missing_barcodes) > 0:
        print("Example missing barcodes:", missing_barcodes[:10])
        raise ValueError(f"{prefix}: metadata are missing barcodes required for AnnData")

    meta = meta.loc[barcodes].copy()

    if not all(meta.index == pd.Index(barcodes)):
        raise ValueError(f"{prefix}: metadata order does not match barcodes")

    var = pd.DataFrame(index=pd.Index([str(x) for x in genes], dtype=object))

    adata = ad.AnnData(
        X=counts.T.tocsr(),
        obs=meta,
        var=var
    )

    adata.var_names_make_unique()

    print("AnnData shape (cells x genes):", adata.shape)

    outfile = os.path.join(outdir, f"{prefix}.h5ad")
    adata.write_h5ad(outfile, compression="gzip")
    print("Saved to:", outfile)

    return adata


if __name__ == "__main__":
    prefixes = [
        "immune"
    ]

    for prefix in prefixes:
        build_anndata(prefix)
