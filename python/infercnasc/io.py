"""
Ensembl gene annotation lookup and simulated CNA parsing.

All network requests are cached locally using requests-cache with a SQLite
backend stored in the platform-appropriate cache directory.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Optional

import pandas as pd
import platformdirs
import requests_cache


_ENSEMBL_URL = "https://rest.ensembl.org/lookup/id"
_BATCH_SIZE = 1000
_CACHE_EXPIRY_SECONDS = 43200  # 12 hours


def _default_cache_path() -> str:
    cache_dir = Path(platformdirs.user_cache_dir("infercnasc"))
    cache_dir.mkdir(parents=True, exist_ok=True)
    return str(cache_dir / "ensembl_cache")


def annotate_genes(
    gene_ids: list[str],
    cache_path: Optional[str] = None,
) -> pd.DataFrame:
    """Fetch genomic coordinates for a list of Ensembl gene IDs.

    Results are cached locally to avoid redundant network requests across runs.
    The cache lives in the platform cache directory by default and expires after
    12 hours. Pass `cache_path` to override the cache file location (useful in
    tests).

    Parameters
    ----------
    gene_ids:
        List of Ensembl gene IDs (e.g. "ENSG00000000001").
    cache_path:
        Path to the SQLite cache file (without extension). Defaults to the
        platform cache directory.

    Returns
    -------
    DataFrame with columns: gene, chrom, start, end.
        Genes not found in Ensembl have NaN for chrom/start/end.
    """
    resolved = cache_path or _default_cache_path()
    session = requests_cache.CachedSession(
        resolved,
        backend="sqlite",
        expire_after=_CACHE_EXPIRY_SECONDS,
    )

    results: dict[str, dict] = {}
    n_batches = math.ceil(len(gene_ids) / _BATCH_SIZE)
    for i in range(n_batches):
        batch = gene_ids[i * _BATCH_SIZE : (i + 1) * _BATCH_SIZE]
        response = session.post(
            _ENSEMBL_URL,
            headers={"Content-Type": "application/json"},
            data=json.dumps({"ids": batch}),
        )
        if response.status_code == 200:
            data = response.json()
            for gid in batch:
                info = data.get(gid, {})
                if info:
                    results[gid] = {
                        "chrom": str(info.get("seq_region_name")),
                        "start": int(info["start"]),
                        "end": int(info["end"]),
                    }

    rows = []
    for gid in gene_ids:
        if gid in results:
            rows.append({"gene": gid, **results[gid]})
        else:
            rows.append({"gene": gid, "chrom": None, "start": None, "end": None})

    return pd.DataFrame(rows)


def parse_simulated(sim_cna_col: pd.Series) -> pd.DataFrame:
    """Parse a column of free-text simulated CNA annotations.

    The expected format per entry is: "chrN:start-end CN X"
    where X is the copy number integer. Copy number 2 is diploid (neutral).
    Copy number below 2 is a loss; above 2 is a gain.

    Parameters
    ----------
    sim_cna_col:
        Series of annotation strings. NaN entries are skipped.

    Returns
    -------
    DataFrame with columns: chrom, start, end, label.
    """
    rows = []
    for entry in sim_cna_col:
        if pd.isna(entry):
            continue
        parts = str(entry).split()
        chrom, coords = parts[0].split(":")
        start, end = map(int, coords.split("-"))
        cn = int(parts[-1])
        label = "loss" if cn < 2 else "gain" if cn > 2 else "neutral"
        rows.append({"chrom": chrom, "start": start, "end": end, "label": label})
    return pd.DataFrame(rows)
