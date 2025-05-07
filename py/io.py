import anndata


def read_h5ad(path: str) -> anndata.AnnData:
    """Read an AnnData object from an H5AD file."""
    return anndata.read_h5ad(path)


def write_h5ad(adata: anndata.AnnData, path: str) -> None:
    """Write an AnnData object to an H5AD file."""
    adata.write_h5ad(path)