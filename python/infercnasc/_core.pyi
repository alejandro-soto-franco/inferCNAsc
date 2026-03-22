from __future__ import annotations
from numpy.typing import NDArray
import numpy as np

def smooth_expression(
    expression: NDArray[np.float64],
    chroms: list[str],
    window_size: int,
) -> NDArray[np.float64]: ...

def find_cnas(
    smoothed: NDArray[np.float64],
    z_score_threshold: float,
) -> tuple[NDArray[np.bool_], NDArray[np.bool_]]: ...

def assign_cnas_to_cells(
    gains: NDArray[np.bool_],
    losses: NDArray[np.bool_],
    chroms: list[str],
    starts: NDArray[np.int64],
    ends: NDArray[np.int64],
    gene_names: list[str],
    min_region_size: int,
) -> list[dict]: ...
