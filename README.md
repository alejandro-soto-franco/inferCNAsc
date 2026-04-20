# inferCNAsc

[![crates.io](https://img.shields.io/crates/v/infercnasc.svg)](https://crates.io/crates/infercnasc)
[![docs.rs](https://docs.rs/infercnasc/badge.svg)](https://docs.rs/infercnasc)
[![PyPI](https://img.shields.io/pypi/v/infercnasc.svg)](https://pypi.org/project/infercnasc/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![CI](https://github.com/alejandro-soto-franco/inferCNAsc/actions/workflows/ci.yml/badge.svg)](https://github.com/alejandro-soto-franco/inferCNAsc/actions)
[![MSRV](https://img.shields.io/badge/MSRV-1.85-blue.svg)](Cargo.toml)

Copy number alteration (CNA) inference from single-cell RNA-seq data. A Rust
core with optional Python bindings via PyO3.

The pipeline is a chromosome-aware sliding-window smoother, per-gene z-score
thresholding, and a parallel run-length merge that assembles per-cell CNA
regions. The Rust core is parallelized over genes (smoothing, z-scoring) and
over cells (region assembly) with `rayon`; the Python layer handles AnnData
adaptation, Ensembl annotation lookup, evaluation, and plotting.

## Installation

```bash
pip install infercnasc
```

Wheels are built for Linux x86_64/aarch64, macOS universal2, and Windows x86_64
against Python's abi3-py310 stable ABI, so a single wheel serves Python 3.10+.

For the Rust API:

```toml
[dependencies]
infercnasc = "0.2"
```

No feature flags are needed for the native Rust API.

## Python

```python
from infercnasc import CNAInferrer
import infercnasc.plot as icplot

inferrer = CNAInferrer.from_anndata(adata)

inferrer = CNAInferrer(window_size=50).fit(expression_matrix, gene_df)

cnas = inferrer.cna_df()
icplot.cna_matrix(inferrer)
```

`gene_df` is a DataFrame with columns `gene`, `chrom`, `start`, `end`.
`infercnasc.io.annotate_genes(gene_ids)` fetches these from Ensembl with a
local `requests-cache` backing store.

Sparse `AnnData.X` is supported natively. Coordinate-annotation filtering
runs on the sparse matrix first, so the eventual dense materialization is
limited to genes that survive annotation; this avoids the standard scRNA
memory blow-up of an unconditional `.toarray()`.

## Rust

```rust
use infercnasc::{smooth_expression, find_cnas, assign_cnas_to_cells, InferError};

let smoothed = smooth_expression(&expression, &chroms, window_size)?;
let (gains, losses) = find_cnas(&smoothed, z_score_threshold);
let cnas = assign_cnas_to_cells(
    &gains, &losses, &chroms, &starts, &ends, &gene_names, min_region_size,
);
```

`smooth_expression` returns `Result<Array2<f64>, InferError>`. `find_cnas` and
`assign_cnas_to_cells` are infallible.

## Pipeline

1. **Gene annotation.** Gene identifiers are resolved to genomic coordinates
   via the Ensembl REST API, with responses cached locally under the
   platform cache directory.
2. **Smoothing.** A sliding-window mean is applied along each chromosome.
   The window resets at chromosome boundaries. Columns are processed in
   parallel.
3. **CNA calling.** Per-gene z-scores are computed across cells. Entries
   above `+z_threshold` are flagged as gains, entries below `-z_threshold`
   as losses. Zero-variance genes are skipped.
4. **Region assembly.** Consecutive flagged genes on the same chromosome are
   merged into `CnaRecord` regions by a parallel per-cell run-length scan.
   Runs shorter than `min_region_size` are dropped.

## Benchmarks

End-to-end pipeline (smoothing + calling + region assembly) on real
public tumor scRNA-seq data: the Tirosh 2016 oligodendroglioma dataset
shipped with the inferCNV R package (184 cells x ~10,000 annotated
genes). A planted-chr1-loss synthetic matrix is also run to show scaling
at larger sizes. Single local run on a Ryzen laptop; your numbers will
differ.

| implementation                | Tirosh (184 x 10,338) | synth (2000 x 10,000) |
|-------------------------------|----------------------:|----------------------:|
| `_core` direct FFI call       |              0.032 s  |              0.62 s   |
| `CNAInferrer.fit` (wrapper)   |              0.111 s  |              0.84 s   |
| `infercnvpy.tl.infercnv`      |              1.134 s  |              1.95 s   |
| pure-numpy reference          |              0.547 s  |             (skip)    |

`_core` direct is the straight FFI call; `CNAInferrer.fit` adds the
coordinate filter, DataFrame sort, and DataFrame assembly around it.
`infercnvpy.tl.infercnv` uses a different algorithmic core (log
fold-change on sparse windows) and is the nearest published
Python-ecosystem comparator. The `pure-numpy reference` is a faithful
per-chromosome cumulative-sum smoothing + z-score reimplementation used
as an apples-to-apples control for the algorithm itself.

Reproduce:

```bash
python benchmarks/compare.py                  # real Tirosh data
python benchmarks/compare.py --synth --cells 2000 --genes 10000
cargo bench --no-default-features             # native Rust criterion
```

## Evaluation

```python
metrics = inferrer.evaluate(simulated_df)
# {"true_positives": ..., "precision": ..., "recall": ..., "f1": ...}
```

Matching is any-overlap on genomic coordinates within the same chromosome
and label. The implementation is O(n_inferred + n_truth) via a chromosome-
and-label-indexed bucket sweep.

## Acknowledgements

A pre-release Python prototype predating this crate was developed with Raeann
Kalinowski and Amy Liu as a 2025 course project at Johns Hopkins; this
repository is a full independent rewrite and is not affiliated with that
coursework.

## License

MIT
