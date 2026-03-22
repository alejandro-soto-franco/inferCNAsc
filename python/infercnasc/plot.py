"""
Visualization functions for CNAInferrer results.

Each function takes a fitted CNAInferrer instance and returns a matplotlib
Figure. Calling any function before fit() raises RuntimeError.
"""
from __future__ import annotations

from typing import TYPE_CHECKING
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.figure

if TYPE_CHECKING:
    from infercnasc.inferrer import CNAInferrer


def _chrom_sort_key(chrom: str) -> int:
    c = str(chrom).replace("chr", "")
    if c.isdigit():
        return int(c)
    return {"X": 23, "Y": 24}.get(c, 25)


def smooth_heatmap(inferrer: CNAInferrer) -> matplotlib.figure.Figure:
    """Heatmap of smoothed gene expression across all cells.

    Parameters
    ----------
    inferrer:
        A fitted CNAInferrer instance.

    Returns
    -------
    matplotlib Figure.
    """
    inferrer._require_fit()
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(inferrer.smoothed_expression(), cmap="coolwarm", aspect="auto")
    fig.colorbar(im, ax=ax, label="Smoothed Expression")
    ax.set_title("Smoothed Gene Expression Across Cells")
    ax.set_ylabel("Cells")
    ax.set_xlabel("Genes")
    fig.tight_layout()
    return fig


def cnas_scatter(inferrer: CNAInferrer) -> matplotlib.figure.Figure:
    """Scatter plot of detected CNA positions along chromosomes.

    Parameters
    ----------
    inferrer:
        A fitted CNAInferrer instance.

    Returns
    -------
    matplotlib Figure.
    """
    inferrer._require_fit()
    df = inferrer.cna_df()
    fig, ax = plt.subplots(figsize=(12, 6))
    for label, color in [("gain", "red"), ("loss", "blue")]:
        sub = df[df["label"] == label]
        ax.scatter(sub["chrom"], sub["start_pos"], c=color, label=label, alpha=0.6)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Start Position")
    ax.set_title("Detected CNAs: Gains and Losses")
    ax.tick_params(axis="x", rotation=90)
    ax.legend()
    fig.tight_layout()
    return fig


def genome_line(inferrer: CNAInferrer) -> matplotlib.figure.Figure:
    """Line plot of average smoothed expression along the genome.

    Parameters
    ----------
    inferrer:
        A fitted CNAInferrer instance.

    Returns
    -------
    matplotlib Figure.
    """
    inferrer._require_fit()
    gi = inferrer._gene_info
    smoothed = inferrer.smoothed_expression()
    avg_expr = smoothed.mean(axis=0)

    sorted_var = gi.sort_values(
        ["chrom", "start"],
        key=lambda col: col.map(_chrom_sort_key) if col.name == "chrom" else col,
    )
    chroms = sorted(sorted_var["chrom"].unique(), key=_chrom_sort_key)

    genome_pos: list[float] = []
    tick_positions: list[float] = []
    tick_labels: list[str] = []
    offset = 0.0

    for chrom in chroms:
        chrom_genes = sorted_var[sorted_var["chrom"] == chrom]
        positions = chrom_genes["start"].values + offset
        genome_pos.extend(positions.tolist())
        tick_positions.append(offset + (positions[-1] - positions[0]) / 2)
        tick_labels.append(str(chrom))
        offset = float(positions[-1]) + 1_000_000

    sorted_idx = sorted_var.index.tolist()
    orig_idx = gi.index.tolist()
    reindex = [orig_idx.index(i) for i in sorted_idx]
    sorted_expr = avg_expr[reindex]

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(genome_pos, sorted_expr)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=90)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Avg Smoothed Expression")
    ax.set_title("Gene Expression Along Genome")
    fig.tight_layout()
    return fig


def cna_matrix(inferrer: CNAInferrer) -> matplotlib.figure.Figure:
    """Matrix heatmap of CNA status per cell and genomic region.

    Gains are shown as +1 (red) and losses as -1 (blue).

    Parameters
    ----------
    inferrer:
        A fitted CNAInferrer instance.

    Returns
    -------
    matplotlib Figure.
    """
    inferrer._require_fit()
    df = inferrer.cna_df()
    if df.empty:
        fig, ax = plt.subplots()
        ax.set_title("No CNAs detected")
        return fig

    n_cells = inferrer.gains().shape[0]
    df = df.copy()
    df["region_id"] = df.apply(
        lambda r: f"{r['chrom']}:{r['start_pos']}-{r['end_pos']} ({r['label']})",
        axis=1,
    )
    unique_regions = sorted(df["region_id"].unique())
    region_index = {r: i for i, r in enumerate(unique_regions)}
    matrix = np.zeros((n_cells, len(unique_regions)))
    for _, row in df.iterrows():
        matrix[int(row["cell"]), region_index[row["region_id"]]] = (
            1 if row["label"] == "gain" else -1
        )

    fig, ax = plt.subplots(figsize=(16, 8))
    im = ax.imshow(matrix, aspect="auto", cmap="bwr", vmin=-1, vmax=1)
    fig.colorbar(im, ax=ax, label="CNA Status (1=gain, -1=loss)")
    ax.set_xticks(range(len(unique_regions)))
    ax.set_xticklabels(unique_regions, rotation=90, fontsize=6)
    ax.set_yticks([])
    ax.set_xlabel("Genomic Region")
    ax.set_ylabel("Cells")
    ax.set_title("Copy Number Alterations per Cell")
    fig.tight_layout()
    return fig
