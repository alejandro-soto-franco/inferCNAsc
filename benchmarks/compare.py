"""
Benchmark inferCNAsc against infercnvpy and a pure-numpy baseline.

Default data source is the Tirosh 2016 oligodendroglioma dataset shipped
with the inferCNV R package (184 cells x ~10,000 genes, real tumor
scRNA-seq with well-characterised 1p/19q co-deletion). Three files are
downloaded on first use and cached under the platform cache directory;
total payload is ~3.2 MB.

A `--synth` option generates a synthetic matrix instead, for when the
Tirosh mirror is unreachable or when you want to sweep matrix sizes
independently of real biology.

For native-Rust reference numbers, run `cargo bench --no-default-features`
(criterion benchmarks in `benches/pipeline.rs`).
"""
from __future__ import annotations

import argparse
import gzip
import io
import time
import urllib.request
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import platformdirs


@dataclass
class Result:
    name: str
    seconds: float
    ok: bool
    note: str = ""


@contextmanager
def timed(name: str, results: list[Result]):
    t0 = time.perf_counter()
    try:
        yield
        elapsed = time.perf_counter() - t0
        results.append(Result(name, elapsed, True))
    except Exception as e:
        elapsed = time.perf_counter() - t0
        results.append(Result(name, elapsed, False, note=f"{type(e).__name__}: {e}"))


_TIROSH_BASE = (
    "https://raw.githubusercontent.com/broadinstitute/infercnv/master/inst/extdata/"
)
_TIROSH_FILES = {
    "counts": "oligodendroglioma_expression_downsampled.counts.matrix.gz",
    "genes": "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt",
}


def _cache_dir() -> Path:
    d = Path(platformdirs.user_cache_dir("infercnasc")) / "benchmark_fixtures" / "tirosh"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _fetch(name: str) -> Path:
    cache = _cache_dir() / _TIROSH_FILES[name]
    if not cache.exists():
        with urllib.request.urlopen(_TIROSH_BASE + _TIROSH_FILES[name]) as resp:
            cache.write_bytes(resp.read())
    return cache


def load_tirosh() -> tuple[np.ndarray, pd.DataFrame]:
    counts_path = _fetch("counts")
    gene_path = _fetch("genes")

    with open(counts_path, "rb") as fh:
        expr_raw = gzip.decompress(fh.read())
    expr_df = pd.read_csv(io.BytesIO(expr_raw), sep="\t", index_col=0)
    expression = expr_df.T.to_numpy(dtype=np.float64)
    gene_symbols = expr_df.index.to_list()

    coords = pd.read_csv(
        gene_path, sep="\t", header=None, names=["gene", "chrom", "start", "end"]
    )
    coords["chrom"] = coords["chrom"].str.replace("chr", "", regex=False)
    coords = coords.drop_duplicates(subset="gene")

    gene_info = (
        pd.DataFrame({"gene": gene_symbols})
        .merge(coords, on="gene", how="left")
        .dropna(subset=["chrom", "start", "end"])
        .reset_index(drop=True)
    )
    expression = expression[:, gene_info.index.to_numpy()]
    return expression, gene_info


def synth_dataset(
    n_cells: int, n_genes: int, n_chroms: int, seed: int
) -> tuple[np.ndarray, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    expression = rng.standard_normal((n_cells, n_genes))
    chrom_of_gene = np.asarray([1 + (g * n_chroms) // n_genes for g in range(n_genes)])
    chr1_mask = chrom_of_gene == 1
    tumor_cells = rng.random(n_cells) < 0.5
    expression[np.ix_(tumor_cells, chr1_mask)] -= 2.0

    gene_info = pd.DataFrame(
        {
            "gene": [f"G{i}" for i in range(n_genes)],
            "chrom": chrom_of_gene.astype(str),
            "start": np.arange(n_genes, dtype=np.int64) * 1000,
            "end": np.arange(n_genes, dtype=np.int64) * 1000 + 999,
        }
    )
    return expression.astype(np.float64), gene_info


def run_infercnasc(expression: np.ndarray, gene_info: pd.DataFrame) -> int:
    from infercnasc import CNAInferrer

    inferrer = CNAInferrer(window_size=100, z_score_threshold=1.0, min_region_size=25).fit(
        expression, gene_info
    )
    return len(inferrer.cna_df())


def run_core_only(expression: np.ndarray, gene_info: pd.DataFrame) -> int:
    from infercnasc import _core

    gi = gene_info.sort_values(["chrom", "start"]).reset_index(drop=True)
    sort_idx = gi.index.to_numpy()
    expr = np.ascontiguousarray(expression[:, sort_idx], dtype=np.float64)
    chroms = gi["chrom"].astype(str).tolist()
    starts = gi["start"].to_numpy(dtype=np.int64)
    ends = gi["end"].to_numpy(dtype=np.int64)
    names = gi["gene"].astype(str).tolist()

    smoothed = _core.smooth_expression(expr, chroms, 100)
    gains, losses = _core.find_cnas(smoothed, 1.0)
    records = _core.assign_cnas_to_cells(gains, losses, chroms, starts, ends, names, 25)
    return len(records)


def run_numpy_baseline(expression: np.ndarray, gene_info: pd.DataFrame) -> int:
    gi = gene_info.sort_values(["chrom", "start"]).reset_index(drop=True)
    sort_idx = gi.index.to_numpy()
    chrom_of_gene = gi["chrom"].to_numpy()
    expr_sorted = expression[:, sort_idx]

    n_cells, n_genes = expr_sorted.shape
    window = 100
    half = window // 2
    smoothed = np.empty_like(expr_sorted)

    chrom_changes = np.flatnonzero(chrom_of_gene[1:] != chrom_of_gene[:-1])
    seg_bounds = np.concatenate([[0], chrom_changes + 1, [n_genes]])
    for seg_i in range(len(seg_bounds) - 1):
        s, e = seg_bounds[seg_i], seg_bounds[seg_i + 1]
        block = expr_sorted[:, s:e]
        n = block.shape[1]
        cumsum = np.concatenate([np.zeros((n_cells, 1)), np.cumsum(block, axis=1)], axis=1)
        for g in range(n):
            w_s = max(0, g - half)
            w_e = min(n, g + half + 1)
            smoothed[:, s + g] = (cumsum[:, w_e] - cumsum[:, w_s]) / (w_e - w_s)

    mean = smoothed.mean(axis=0, keepdims=True)
    std = smoothed.std(axis=0, ddof=0, keepdims=True)
    std_nonzero = np.where(std > 0, std, 1.0)
    z = (smoothed - mean) / std_nonzero
    z = np.where(std > 0, z, 0.0)

    gains = z > 1.0
    losses = z < -1.0

    regions = 0
    chrom_arr = chrom_of_gene
    for cell in range(n_cells):
        for matrix in (gains, losses):
            run_start = None
            for g in range(n_genes + 1):
                chrom_break = (
                    run_start is not None
                    and g < n_genes
                    and chrom_arr[g] != chrom_arr[run_start]
                )
                boundary = (
                    g == n_genes
                    or (g < n_genes and not matrix[cell, g])
                    or chrom_break
                )
                if boundary:
                    if run_start is not None and g - run_start >= 25:
                        regions += 1
                    run_start = None
                    if g < n_genes and matrix[cell, g]:
                        run_start = g
                elif run_start is None and g < n_genes and matrix[cell, g]:
                    run_start = g
    return regions


def run_infercnvpy(expression: np.ndarray, gene_info: pd.DataFrame) -> int:
    import anndata as ad
    import infercnvpy as cnv

    var = gene_info.rename(columns={"chrom": "chromosome"}).copy()
    var["chromosome"] = "chr" + var["chromosome"].astype(str)
    var = var.set_index("gene")
    adata = ad.AnnData(X=expression, var=var)
    cnv.tl.infercnv(adata, window_size=100, exclude_chromosomes=[], inplace=True)
    return int(adata.obsm["X_cnv"].nnz)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--synth",
        action="store_true",
        help="use synthetic matrix instead of real Tirosh data",
    )
    parser.add_argument("--cells", type=int, default=500)
    parser.add_argument("--genes", type=int, default=5000)
    parser.add_argument("--chroms", type=int, default=23)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--skip-infercnvpy", action="store_true")
    parser.add_argument("--skip-numpy-baseline", action="store_true")
    args = parser.parse_args()

    if args.synth:
        print(
            f"Synthetic dataset: {args.cells} cells x {args.genes} genes, "
            f"{args.chroms} chromosomes, seed={args.seed}"
        )
        expression, gene_info = synth_dataset(
            args.cells, args.genes, args.chroms, args.seed
        )
    else:
        print(
            "Real dataset: Tirosh 2016 oligodendroglioma "
            "(broadinstitute/infercnv/inst/extdata)"
        )
        expression, gene_info = load_tirosh()
        print(
            f"  {expression.shape[0]} cells x {expression.shape[1]} "
            f"annotated genes, {gene_info['chrom'].nunique()} chromosomes"
        )

    results: list[Result] = []

    with timed("infercnasc (CNAInferrer.fit)", results):
        n1 = run_infercnasc(expression, gene_info)
    print(f"  infercnasc: {n1} CNA regions")

    with timed("infercnasc (_core direct)", results):
        n2 = run_core_only(expression, gene_info)
    print(f"  _core direct: {n2} CNA regions")

    if not args.skip_numpy_baseline:
        with timed("numpy baseline (pure Python)", results):
            n3 = run_numpy_baseline(expression, gene_info)
        print(f"  numpy baseline: {n3} CNA regions")

    if not args.skip_infercnvpy:
        try:
            import infercnvpy  # noqa: F401

            with timed("infercnvpy.tl.infercnv", results):
                n4 = run_infercnvpy(expression, gene_info)
            print(f"  infercnvpy: {n4} non-zero CNV entries")
        except ImportError:
            results.append(Result("infercnvpy.tl.infercnv", 0.0, False, "not installed"))

    print()
    print(f"{'implementation':<36}{'seconds':>10}  notes")
    print("-" * 72)
    baseline = next(
        (r.seconds for r in results if r.name == "infercnasc (_core direct)" and r.ok),
        None,
    )
    for r in results:
        if r.ok:
            ratio = f"{r.seconds / baseline:5.2f}x" if baseline else ""
            print(f"{r.name:<36}{r.seconds:>10.3f}  {ratio}")
        else:
            print(f"{r.name:<36}{'fail':>10}  {r.note}")


if __name__ == "__main__":
    main()
