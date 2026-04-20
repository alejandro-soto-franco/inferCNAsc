from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

def smooth_expression(
    expression: NDArray[np.float64],
    chroms: list[str],
    window_size: int,
) -> NDArray[np.float64]: ...
def find_cnas(
    smoothed: NDArray[np.float64],
    z_score_threshold: float,
) -> tuple[NDArray[np.uint8], NDArray[np.uint8]]:
    """Returns uint8 arrays with 0/1 values. View as np.bool_ for boolean
    semantics; pyo3-numpy maps Rust `bool` to numpy uint8 over the FFI
    boundary."""
    ...

def assign_cnas_to_cells(
    gains: NDArray[np.uint8],
    losses: NDArray[np.uint8],
    chroms: list[str],
    starts: NDArray[np.int64],
    ends: NDArray[np.int64],
    gene_names: list[str],
    min_region_size: int,
) -> list[dict]: ...
