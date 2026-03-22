import pytest
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # non-interactive backend, must be set before importing pyplot
import matplotlib.pyplot as plt
from infercnasc.inferrer import CNAInferrer
import infercnasc.plot as icplot


def make_fitted():
    rng = np.random.default_rng(1)
    expr = rng.poisson(5.0, size=(10, 30)).astype(np.float64)
    gene_info = pd.DataFrame({
        "gene": [f"G{i}" for i in range(30)],
        "chrom": ["1"] * 15 + ["2"] * 15,
        "start": list(range(0, 30_000, 1000)),
        "end":   list(range(999, 30_000, 1000)),
    })
    return CNAInferrer(window_size=5).fit(expr, gene_info)


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def test_pre_fit_smooth_heatmap_raises():
    with pytest.raises(RuntimeError):
        icplot.smooth_heatmap(CNAInferrer())


def test_pre_fit_cnas_scatter_raises():
    with pytest.raises(RuntimeError):
        icplot.cnas_scatter(CNAInferrer())


def test_pre_fit_genome_line_raises():
    with pytest.raises(RuntimeError):
        icplot.genome_line(CNAInferrer())


def test_pre_fit_cna_matrix_raises():
    with pytest.raises(RuntimeError):
        icplot.cna_matrix(CNAInferrer())


def test_smooth_heatmap_returns_figure():
    fig = icplot.smooth_heatmap(make_fitted())
    assert isinstance(fig, plt.Figure)


def test_cnas_scatter_returns_figure():
    fig = icplot.cnas_scatter(make_fitted())
    assert isinstance(fig, plt.Figure)


def test_genome_line_returns_figure():
    fig = icplot.genome_line(make_fitted())
    assert isinstance(fig, plt.Figure)


def test_cna_matrix_returns_figure():
    fig = icplot.cna_matrix(make_fitted())
    assert isinstance(fig, plt.Figure)
