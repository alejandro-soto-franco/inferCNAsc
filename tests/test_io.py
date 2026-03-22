import pytest
import pandas as pd
import responses as resp_lib
import json
from infercnasc.io import annotate_genes, parse_simulated


@resp_lib.activate
def test_annotate_genes_returns_correct_schema(tmp_path):
    """Mocked Ensembl response produces a DataFrame with the right columns."""
    mock_response = {
        "ENSG00000000001": {
            "seq_region_name": "1",
            "start": 100,
            "end": 200,
        }
    }
    resp_lib.add(
        resp_lib.POST,
        "https://rest.ensembl.org/lookup/id",
        json=mock_response,
        status=200,
    )
    result = annotate_genes(["ENSG00000000001"], cache_path=str(tmp_path / "cache.sqlite"))
    assert isinstance(result, pd.DataFrame)
    assert set(result.columns) >= {"gene", "chrom", "start", "end"}
    assert result.loc[0, "chrom"] == "1"
    assert result.loc[0, "start"] == 100


def test_parse_simulated_loss():
    series = pd.Series(["chr1:100-200 CN 0"])
    result = parse_simulated(series)
    assert result.loc[0, "label"] == "loss"
    assert result.loc[0, "chrom"] == "chr1"
    assert result.loc[0, "start"] == 100
    assert result.loc[0, "end"] == 200


def test_parse_simulated_gain():
    series = pd.Series(["chr3:50-150 CN 4"])
    result = parse_simulated(series)
    assert result.loc[0, "label"] == "gain"


def test_parse_simulated_neutral():
    series = pd.Series(["chr2:10-20 CN 2"])
    result = parse_simulated(series)
    assert result.loc[0, "label"] == "neutral"


def test_parse_simulated_skips_na():
    series = pd.Series([None, "chr1:100-200 CN 1"])
    result = parse_simulated(series)
    assert len(result) == 1
