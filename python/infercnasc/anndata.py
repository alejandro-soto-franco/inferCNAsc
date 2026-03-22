"""
Adapter for converting AnnData objects into the array inputs expected by CNAInferrer.

This module handles the AnnData-specific logic so that the Rust core and the
CNAInferrer class remain independent of anndata.
"""
from __future__ import annotations

from typing import TYPE_CHECKING
import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata as ad


def extract_from_anndata(adata: ad.AnnData) -> tuple[np.ndarray, pd.DataFrame]:
    """Extract expression matrix and gene annotation DataFrame from an AnnData object.

    If `adata.var` already contains non-null `chrom`, `start`, and `end`
    columns, the Ensembl lookup is skipped. Otherwise, `annotate_genes` from
    `infercnasc.io` is called to fetch the missing annotations.

    Parameters
    ----------
    adata:
        AnnData object with cells as rows and genes as columns.

    Returns
    -------
    expression:
        Dense float64 array of shape (n_cells, n_genes).
    gene_info:
        DataFrame with columns: gene, chrom, start, end.
    """
    import scipy.sparse

    X = adata.X
    if scipy.sparse.issparse(X):
        expression = np.asarray(X.todense(), dtype=np.float64)
    else:
        expression = np.asarray(X, dtype=np.float64)

    var = adata.var.copy()
    var.index.name = "gene"
    var = var.reset_index()

    annotation_cols = {"chrom", "start", "end"}
    has_annotations = (
        annotation_cols <= set(var.columns)
        and var[list(annotation_cols)].notna().all(axis=None)
    )

    if not has_annotations:
        from infercnasc.io import annotate_genes
        gene_ids = var["gene"].tolist()
        fetched = annotate_genes(gene_ids)
        var = var.merge(fetched, on="gene", how="left")

    return expression, var[["gene", "chrom", "start", "end"]]
