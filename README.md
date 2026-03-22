# inferCNAsc

[![crates.io](https://img.shields.io/crates/v/infercnasc.svg)](https://crates.io/crates/infercnasc)
[![docs.rs](https://docs.rs/infercnasc/badge.svg)](https://docs.rs/infercnasc)
[![PyPI](https://img.shields.io/pypi/v/infercnasc.svg)](https://pypi.org/project/infercnasc/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Tests](https://github.com/alejandro-soto-franco/inferCNAsc/actions/workflows/ci.yml/badge.svg)](https://github.com/alejandro-soto-franco/inferCNAsc/actions)
[![MSRV](https://img.shields.io/badge/MSRV-1.70-blue.svg)](Cargo.toml)

Copy number alteration (CNA) inference from single-cell RNA-seq data.
A Rust core with optional Python bindings via PyO3.

## Python

```bash
pip install infercnasc
```

```python
from infercnasc import CNAInferrer
import infercnasc.plot as icplot

# From an AnnData object (Ensembl lookup runs automatically when needed)
inferrer = CNAInferrer.from_anndata(adata)

# From raw arrays (no AnnData dependency required)
inferrer = CNAInferrer(window_size=50).fit(expression_matrix, gene_df)

cnas = inferrer.cna_df()       # DataFrame of detected CNA regions
icplot.cna_matrix(inferrer)    # per-cell CNA heatmap
```

`gene_df` is a DataFrame with columns `gene`, `chrom`, `start`, `end`.
Use `infercnasc.io.annotate_genes(gene_ids)` to fetch these from Ensembl.

## Rust

```toml
[dependencies]
infercnasc = "0.1"
```

No feature flags are needed for the native Rust API.

```rust
use infercnasc::{smooth_expression, find_cnas, assign_cnas_to_cells, InferError};

let smoothed = smooth_expression(&expression, &chroms, window_size)?;
let (gains, losses) = find_cnas(&smoothed, z_score_threshold);
let cnas = assign_cnas_to_cells(
    &gains, &losses, &chroms, &starts, &ends, &gene_names, min_region_size,
);
```

`smooth_expression` returns `Result<Array2<f64>, InferError>`.
`find_cnas` and `assign_cnas_to_cells` are infallible.

## How it works

1. **Gene annotation**: gene IDs are mapped to genomic coordinates via the
   Ensembl REST API (results cached locally with `requests-cache`).
2. **Smoothing**: a sliding-window mean is applied across neighboring genes
   within each chromosome. The window resets at chromosome boundaries.
3. **CNA calling**: z-scores are computed per gene across cells. Genes above
   the threshold are flagged as gains; genes below are flagged as losses.
4. **Region assembly**: consecutive flagged genes on the same chromosome are
   merged into CNA regions using a run-length scan.

## Authors

**Alejandro J. Soto Franco** (primary author)

The algorithm design and original Python prototype (v0.2) were co-developed with
Raeann Kalinowski and Amy Liu as a final project for 580.447 Computational Stem
Cell Biology, Spring 2025, Johns Hopkins University Department of Biomedical
Engineering. This crate is a full independent rewrite and is no longer associated
with that course or developed for academic submission purposes.

## Citation

> Soto Franco A.J., Kalinowski R., Liu A. *inferCNAsc: a Python toolkit for
> copy number inference from single-cell transcriptomes*. In preparation (2025).

## License

MIT
