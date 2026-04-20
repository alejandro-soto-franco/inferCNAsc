"""
Adapter for converting AnnData objects into the array inputs expected by CNAInferrer.

This module handles the AnnData-specific logic so that the Rust core and the
CNAInferrer class remain independent of anndata.
"""
from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata as ad


_DENSE_WARN_CELLS_X_GENES = 5e8


def extract_from_anndata(adata: ad.AnnData) -> tuple[np.ndarray, pd.DataFrame]:
    """Extract expression matrix and gene annotation DataFrame from an AnnData object.

    If `adata.var` already contains non-null `chrom`, `start`, and `end`
    columns, the Ensembl lookup is skipped. Otherwise, `annotate_genes` from
    `infercnasc.io` is called to fetch the missing annotations.

    For sparse inputs, densification is performed only on the columns that
    survive coordinate-annotation filtering, which avoids the typical scRNA
    memory blow-up (a 95%-sparse 10k x 30k matrix is ~100 MB sparse but ~2.4
    GB dense).

    Parameters
    ----------
    adata:
        AnnData object with cells as rows and genes as columns.

    Returns
    -------
    expression:
        Dense float64 array of shape (n_cells, n_valid_genes).
    gene_info:
        DataFrame with columns: gene, chrom, start, end, aligned to the
        columns of `expression`.
    """
    import scipy.sparse as sp

    var = adata.var.copy()
    var.index.name = "gene"
    var = var.reset_index()

    annotation_cols = {"chrom", "start", "end"}
    has_annotations = annotation_cols <= set(var.columns) and (
        var[list(annotation_cols)].notna().all(axis=None)
    )

    if not has_annotations:
        from infercnasc.io import annotate_genes

        gene_ids = var["gene"].tolist()
        fetched = annotate_genes(gene_ids)
        var = var.merge(fetched, on="gene", how="left")

    gene_info = var[["gene", "chrom", "start", "end"]]
    valid_mask = gene_info[["chrom", "start", "end"]].notna().all(axis=1)
    gene_info = gene_info[valid_mask].reset_index(drop=True)
    valid_idx = np.flatnonzero(valid_mask.to_numpy())

    X = adata.X
    n_cells = X.shape[0]
    n_valid = len(valid_idx)

    if sp.issparse(X):
        # Slice sparse columns first, then densify only what we need.
        X_valid = X[:, valid_idx]
        if n_cells * n_valid > _DENSE_WARN_CELLS_X_GENES:
            warnings.warn(
                f"Densifying {n_cells} x {n_valid} matrix (~"
                f"{n_cells * n_valid * 8 / 1e9:.1f} GB float64). "
                "Consider subsetting genes or cells before calling "
                "CNAInferrer.from_anndata().",
                stacklevel=2,
            )
        expression = np.asarray(X_valid.toarray(), dtype=np.float64)
    else:
        expression = np.asarray(X[:, valid_idx], dtype=np.float64)

    return expression, gene_info
