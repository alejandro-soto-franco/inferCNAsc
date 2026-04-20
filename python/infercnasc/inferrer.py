"""
CNAInferrer: sklearn-style interface for copy number alteration inference.

The heavy computation (smoothing, z-scoring, cell assignment) runs in the
Rust core via `infercnasc._core`. The Python layer handles data preparation,
AnnData adaptation, and evaluation.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Optional

import numpy as np
import pandas as pd

from infercnasc import _core


_NOT_FIT_MSG = "Call fit() before accessing results."


class CNAInferrer:
    """Infer copy number alterations from single-cell RNA-seq expression data.

    Parameters
    ----------
    window_size:
        Number of genes in the sliding-window average. Default 50.
    z_score_threshold:
        Z-score cutoff for calling a gain or loss. Default 1.0.
    min_region_size:
        Minimum number of consecutive genes required to report a CNA region.
        Default 3.
    """

    def __init__(
        self,
        window_size: int = 50,
        z_score_threshold: float = 1.0,
        min_region_size: int = 3,
    ) -> None:
        self.window_size = window_size
        self.z_score_threshold = z_score_threshold
        self.min_region_size = min_region_size
        self._fitted = False
        self._smoothed: Optional[np.ndarray] = None
        self._gains: Optional[np.ndarray] = None
        self._losses: Optional[np.ndarray] = None
        self._cna_df: Optional[pd.DataFrame] = None
        self._gene_info: Optional[pd.DataFrame] = None

    @classmethod
    def from_anndata(cls, adata, **kwargs) -> CNAInferrer:
        """Construct a CNAInferrer from an AnnData object and call fit().

        Gene coordinate annotations are read from `adata.var` if the columns
        `chrom`, `start`, and `end` are present and fully non-null. Otherwise
        the Ensembl REST API is queried automatically.

        Parameters
        ----------
        adata:
            AnnData object with cells as rows and genes as columns.
        **kwargs:
            Passed to the CNAInferrer constructor (window_size, z_score_threshold,
            min_region_size).
        """
        from infercnasc.anndata import extract_from_anndata

        expression, gene_info = extract_from_anndata(adata)
        return cls(**kwargs).fit(expression, gene_info)

    def fit(
        self,
        expression: np.ndarray,
        gene_info: pd.DataFrame,
    ) -> CNAInferrer:
        """Run smoothing, CNA calling, and cell assignment.

        Parameters
        ----------
        expression:
            Shape (n_cells, n_genes), dtype float64. Raw or normalized counts.
        gene_info:
            DataFrame with columns: gene, chrom, start, end. One row per gene,
            parallel to columns of `expression`.

        Returns
        -------
        self, for method chaining.
        """
        valid = gene_info[["chrom", "start", "end"]].notna().all(axis=1)
        gene_info = gene_info[valid].copy()
        expression = expression[:, valid.values]

        sorted_gene_info = gene_info.sort_values(["chrom", "start"])
        sort_idx = sorted_gene_info.index.to_numpy()
        gene_info = sorted_gene_info.reset_index(drop=True)
        expression = expression[:, sort_idx]

        chroms = gene_info["chrom"].astype(str).tolist()
        starts = gene_info["start"].to_numpy(dtype=np.int64)
        ends = gene_info["end"].to_numpy(dtype=np.int64)
        names = gene_info["gene"].astype(str).tolist()

        expression = np.ascontiguousarray(expression, dtype=np.float64)

        self._smoothed = _core.smooth_expression(expression, chroms, self.window_size)
        gains_raw, losses_raw = _core.find_cnas(self._smoothed, self.z_score_threshold)
        # pyo3-numpy maps Rust `bool` to numpy uint8; reinterpret as bool_
        # with a zero-copy view for user-facing attributes.
        self._gains = gains_raw.view(np.bool_)
        self._losses = losses_raw.view(np.bool_)
        raw_records = _core.assign_cnas_to_cells(
            gains_raw,
            losses_raw,
            chroms,
            starts,
            ends,
            names,
            self.min_region_size,
        )
        self._cna_df = pd.DataFrame(raw_records)
        self._gene_info = gene_info
        self._fitted = True
        return self

    def _require_fit(self) -> None:
        if not self._fitted:
            raise RuntimeError(_NOT_FIT_MSG)

    def cna_df(self) -> pd.DataFrame:
        """Return the per-cell CNA region DataFrame.

        Columns: cell, chrom, start_gene, end_gene, start_pos, end_pos, label.
        """
        self._require_fit()
        return self._cna_df  # type: ignore[return-value]

    def smoothed_expression(self) -> np.ndarray:
        """Return the smoothed expression matrix, shape (n_cells, n_genes)."""
        self._require_fit()
        return self._smoothed  # type: ignore[return-value]

    def gains(self) -> np.ndarray:
        """Return the gains boolean array, shape (n_cells, n_genes)."""
        self._require_fit()
        return self._gains  # type: ignore[return-value]

    def losses(self) -> np.ndarray:
        """Return the losses boolean array, shape (n_cells, n_genes)."""
        self._require_fit()
        return self._losses  # type: ignore[return-value]

    def evaluate(self, simulated_df: pd.DataFrame) -> dict:
        """Compare inferred CNAs to a known ground-truth DataFrame.

        Matching is any-overlap on genomic coordinates within the same
        chromosome and label. Runs in O(n_inferred + n_truth) for the common
        case via a chrom+label-indexed bucket sweep, not O(n x m).

        Parameters
        ----------
        simulated_df:
            DataFrame with columns: chrom, start, end, label.

        Returns
        -------
        dict with keys: true_positives, false_positives, false_negatives,
        precision, recall, f1.
        """
        self._require_fit()
        inferred = self._cna_df
        assert inferred is not None

        buckets: dict[tuple[str, str], list[tuple[int, int, int]]] = defaultdict(list)
        for j, s_chrom, s_start, s_end, s_label in zip(
            simulated_df.index,
            simulated_df["chrom"].astype(str),
            simulated_df["start"].astype(int),
            simulated_df["end"].astype(int),
            simulated_df["label"].astype(str),
            strict=False,
        ):
            buckets[(s_chrom, s_label)].append((int(s_start), int(s_end), int(j)))

        matched: set[int] = set()
        tp = 0
        fp = 0

        for inf in inferred.itertuples(index=False):
            key = (str(inf.chrom), str(inf.label))
            hit = False
            for s_start, s_end, j in buckets.get(key, []):
                if not (s_end < inf.start_pos or s_start > inf.end_pos):
                    tp += 1
                    matched.add(j)
                    hit = True
                    break
            if not hit:
                fp += 1

        fn = len(simulated_df) - len(matched)
        eps = 1e-9
        precision = tp / (tp + fp + eps)
        recall = tp / (tp + fn + eps)
        f1 = 2 * precision * recall / (precision + recall + eps)
        return {
            "true_positives": tp,
            "false_positives": fp,
            "false_negatives": fn,
            "precision": precision,
            "recall": recall,
            "f1": f1,
        }

    def assign_to_clusters(self, cluster_labels: pd.Series) -> pd.DataFrame:
        """Group cell-level CNAs by cluster label.

        Parameters
        ----------
        cluster_labels:
            Series of cluster assignments indexed 0..n_cells-1.

        Returns
        -------
        DataFrame grouped by (cluster, chrom, start_pos, end_pos, label)
        with a count column.
        """
        self._require_fit()
        df = self._cna_df.copy()
        df["cluster"] = df["cell"].map(cluster_labels)
        return (
            df.groupby(["cluster", "chrom", "start_pos", "end_pos", "label"])
            .size()
            .reset_index(name="count")
        )
