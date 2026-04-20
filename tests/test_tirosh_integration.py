"""
End-to-end integration test against the Tirosh 2016 oligodendroglioma
dataset shipped with the inferCNV R package.

Downloads three files on demand (~3.2 MB total), runs the full
CNAInferrer pipeline, and asserts that detectable CNA regions are
emitted for the malignant-cell subset. 1p/19q co-deletion is the
diagnostic hallmark of oligodendroglioma, so we expect losses on
chromosomes 1 and 19 in tumor cells.

Gated behind `@pytest.mark.integration`; skipped by default CI.
Run with `pytest -m integration tests/test_tirosh_integration.py -v`.
"""
from __future__ import annotations

import gzip
import io
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import platformdirs
import pytest

from infercnasc import CNAInferrer


_BASE = "https://raw.githubusercontent.com/broadinstitute/infercnv/master/inst/extdata/"
_FILES = {
    "counts": "oligodendroglioma_expression_downsampled.counts.matrix.gz",
    "annotations": "oligodendroglioma_annotations_downsampled.txt",
    "genes": "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt",
}


def _cache_dir() -> Path:
    d = Path(platformdirs.user_cache_dir("infercnasc")) / "test_fixtures" / "tirosh"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _fetch(name: str) -> Path:
    cache = _cache_dir() / _FILES[name]
    if not cache.exists():
        with urllib.request.urlopen(_BASE + _FILES[name]) as response:
            cache.write_bytes(response.read())
    return cache


def _load_tirosh() -> tuple[np.ndarray, pd.DataFrame, pd.Series]:
    counts_path = _fetch("counts")
    ann_path = _fetch("annotations")
    gene_path = _fetch("genes")

    with open(counts_path, "rb") as fh:
        expr_raw = gzip.decompress(fh.read())
    expr_df = pd.read_csv(io.BytesIO(expr_raw), sep="\t", index_col=0)
    expression = expr_df.T.to_numpy(dtype=np.float64)
    cell_ids = expr_df.columns.to_list()
    gene_symbols = expr_df.index.to_list()

    annotations = pd.read_csv(
        ann_path, sep="\t", header=None, names=["cell", "group"]
    ).set_index("cell")
    groups = annotations.reindex(cell_ids)["group"]

    coords = pd.read_csv(
        gene_path, sep="\t", header=None, names=["gene", "chrom", "start", "end"]
    )
    coords["chrom"] = coords["chrom"].str.replace("chr", "", regex=False)
    coords = coords.drop_duplicates(subset="gene")

    gene_info = pd.DataFrame({"gene": gene_symbols}).merge(
        coords, on="gene", how="left"
    )
    return expression, gene_info, groups


@pytest.mark.integration
def test_tirosh_oligodendroglioma_end_to_end():
    expression, gene_info, groups = _load_tirosh()

    assert expression.shape[0] == len(groups), "cell count aligned"
    assert expression.shape[1] == len(gene_info), "gene count aligned"

    valid = gene_info[["chrom", "start", "end"]].notna().all(axis=1)
    assert valid.sum() > 500, "many genes should have coordinate annotations"

    inferrer = CNAInferrer(
        window_size=100, z_score_threshold=1.0, min_region_size=25
    ).fit(expression, gene_info)

    cnas = inferrer.cna_df()
    assert not cnas.empty, "pipeline must emit at least one CNA region"

    tumor_mask = groups.str.startswith("malignant", na=False)
    tumor_cells = np.flatnonzero(tumor_mask.to_numpy())
    assert len(tumor_cells) > 0, "expected malignant cells in dataset"
    tumor_cnas = cnas[cnas["cell"].isin(tumor_cells)]

    assert "loss" in tumor_cnas["label"].unique(), (
        "malignant cells should show at least one loss region "
        "(1p/19q hallmark of oligodendroglioma)"
    )

    tumor_chroms = set(tumor_cnas[tumor_cnas["label"] == "loss"]["chrom"].astype(str))
    assert tumor_chroms & {"1", "19"}, (
        f"expected chr1 or chr19 losses in tumor cells; got chromosomes {tumor_chroms}"
    )

    expected_cols = {
        "cell",
        "chrom",
        "start_gene",
        "end_gene",
        "start_pos",
        "end_pos",
        "label",
    }
    assert expected_cols <= set(cnas.columns)
