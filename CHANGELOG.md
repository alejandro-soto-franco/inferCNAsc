# Changelog

## [0.2.0] - 2026-04-20

### Changed
- MSRV bumped to Rust 1.85; edition upgraded to 2024.
- `pyo3` 0.22 → 0.28, `numpy` 0.22 → 0.28, `thiserror` 1 → 2.
- Python FFI now returns `np.bool_` arrays directly from `find_cnas`; the
  former `u8` intermediate and the Python-side `.astype(bool)` cast are
  gone.

### Added
- `find_cnas` parallelized per gene via `ndarray::Zip::par_for_each`.
- `assign_cnas_to_cells` parallelized per cell via rayon; output order is
  deterministic (cells ascending, gains before losses, position ascending).
- Sparse-aware `extract_from_anndata`: coordinate-annotation filtering is
  applied to the sparse matrix first, so dense materialization is limited
  to annotated genes. A warning is emitted if the post-filter dense size
  exceeds ~4 GB.
- `evaluate()` rewritten as an O(n_inferred + n_truth) chrom-and-label bucket
  sweep; the old nested `iterrows` loop is gone.
- `genome_line()` plot: O(n²) `list.index` reindex replaced with
  `np.lexsort`.
- CI gains `cargo fmt --check` and `cargo clippy -D warnings`; the Python
  job now skips tests marked `integration` so it stays hermetic.
- End-to-end integration test against the Tirosh 2016 oligodendroglioma
  dataset shipped with the inferCNV R package (opt-in via
  `pytest -m integration`).
- Criterion benchmarks (`benches/pipeline.rs`) covering each stage and
  the full end-to-end pipeline at 500 x 2000 and 2000 x 10000.
- `benchmarks/compare.py`: apples-to-apples wall-clock comparison of
  `CNAInferrer.fit`, direct `_core` FFI, a pure-numpy reference, and
  `infercnvpy.tl.infercnv`. Default data source is the real Tirosh 2016
  oligodendroglioma matrix (broadinstitute/infercnv extdata); `--synth`
  swaps in a parameterized synthetic matrix for size sweeps.
- `ndarray` bumped 0.16 -> 0.17 to align with the `numpy 0.28` transitive;
  Cargo.lock is gitignored, so CI was picking up 0.17 from numpy while
  Cargo.toml pinned 0.16.

### Fixed
- README authorship section reduced to a single legacy acknowledgement;
  MSRV badge and citation year brought in line with the codebase.

## [0.1.1] - prior

See git history.
