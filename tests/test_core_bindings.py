import numpy as np
from infercnasc import _core


def test_smooth_expression_shape():
    expr = np.ones((5, 10), dtype=np.float64)
    chroms = ["1"] * 10
    result = _core.smooth_expression(expr, chroms, 3)
    assert result.shape == (5, 10)
    assert result.dtype == np.float64


def test_find_cnas_shape():
    smoothed = np.random.default_rng(0).standard_normal((5, 10))
    gains, losses = _core.find_cnas(smoothed, 1.0)
    assert gains.shape == (5, 10)


def test_assign_cnas_to_cells_returns_list_of_dicts():
    # find_cnas returns u8 arrays; build u8 gains/losses directly
    gains  = np.zeros((3, 6), dtype=np.uint8)
    losses = np.zeros((3, 6), dtype=np.uint8)
    gains[0, 0:4] = 1   # 4-gene gain run for cell 0
    chroms = ["1"] * 6
    starts = np.arange(6, dtype=np.int64) * 1000
    ends   = starts + 999
    names  = [f"G{i}" for i in range(6)]
    records = _core.assign_cnas_to_cells(gains, losses, chroms, starts, ends, names, 3)
    assert isinstance(records, list)
    assert len(records) == 1
    assert set(records[0].keys()) == {"cell", "chrom", "start_gene", "end_gene", "start_pos", "end_pos", "label"}
