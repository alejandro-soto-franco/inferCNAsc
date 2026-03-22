"""
infercnasc: copy number alteration inference from single-cell RNA-seq data.

Quick start (Python):

    from infercnasc import CNAInferrer
    import infercnasc.plot as icplot

    inferrer = CNAInferrer.from_anndata(adata)
    # or
    inferrer = CNAInferrer(window_size=50).fit(expression_matrix, gene_df)
    cnas = inferrer.cna_df()
    icplot.cna_matrix(inferrer)

Quick start (Rust): add infercnasc = "0.1" to Cargo.toml. No feature flags
are needed for the pure Rust API.
"""
from infercnasc.inferrer import CNAInferrer
from infercnasc.io import annotate_genes, parse_simulated

__all__ = ["CNAInferrer", "annotate_genes", "parse_simulated"]
