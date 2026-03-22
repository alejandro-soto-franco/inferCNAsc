import pytest
import numpy as np
import pandas as pd
from infercnasc.inferrer import CNAInferrer


def make_synthetic(n_cells: int = 20, n_genes: int = 50):
    """Return a small synthetic expression matrix and gene_info DataFrame."""
    rng = np.random.default_rng(42)
    expression = rng.poisson(5.0, size=(n_cells, n_genes)).astype(np.float64)
    chroms = ["1"] * 25 + ["2"] * 25
    gene_info = pd.DataFrame({
        "gene": [f"G{i}" for i in range(n_genes)],
        "chrom": chroms,
        "start": list(range(0, n_genes * 1000, 1000)),
        "end":   list(range(999, n_genes * 1000, 1000)),
    })
    return expression, gene_info


def test_pre_fit_accessors_raise_runtime_error():
    inferrer = CNAInferrer()
    for method in ("cna_df", "gains", "losses", "smoothed_expression"):
        with pytest.raises(RuntimeError):
            getattr(inferrer, method)()


def test_pre_fit_evaluate_raises():
    with pytest.raises(RuntimeError):
        CNAInferrer().evaluate(pd.DataFrame())


def test_pre_fit_assign_to_clusters_raises():
    with pytest.raises(RuntimeError):
        CNAInferrer().assign_to_clusters(pd.Series(dtype=str))


def test_fit_returns_self():
    expr, gene_info = make_synthetic()
    inferrer = CNAInferrer()
    result = inferrer.fit(expr, gene_info)
    assert result is inferrer


def test_cna_df_has_correct_schema():
    expr, gene_info = make_synthetic()
    inferrer = CNAInferrer(window_size=5, z_score_threshold=1.0).fit(expr, gene_info)
    df = inferrer.cna_df()
    assert isinstance(df, pd.DataFrame)
    required = {"cell", "chrom", "start_gene", "end_gene", "start_pos", "end_pos", "label"}
    assert required <= set(df.columns), f"missing columns: {required - set(df.columns)}"


def test_gains_losses_shape():
    expr, gene_info = make_synthetic()
    inferrer = CNAInferrer(window_size=5).fit(expr, gene_info)
    assert inferrer.gains().shape == (20, 50)
    assert inferrer.losses().shape == (20, 50)
    assert inferrer.gains().dtype == bool


def test_assign_to_clusters_returns_dataframe():
    expr, gene_info = make_synthetic()
    inferrer = CNAInferrer(window_size=5).fit(expr, gene_info)
    labels = pd.Series(["A"] * 10 + ["B"] * 10)
    result = inferrer.assign_to_clusters(labels)
    assert isinstance(result, pd.DataFrame)


def test_from_anndata_integration():
    """from_anndata() extracts expression and annotated gene_info from AnnData."""
    import anndata
    rng = np.random.default_rng(0)
    X = rng.poisson(3.0, size=(10, 20)).astype(np.float64)
    var = pd.DataFrame({
        "chrom": (["1"] * 10 + ["2"] * 10),
        "start": list(range(0, 20_000, 1000)),
        "end":   list(range(999, 20_000, 1000)),
    }, index=[f"ENSG{i:011d}" for i in range(20)])
    adata = anndata.AnnData(X=X, var=var)

    inferrer = CNAInferrer.from_anndata(adata, window_size=3)
    df = inferrer.cna_df()
    assert isinstance(df, pd.DataFrame)
