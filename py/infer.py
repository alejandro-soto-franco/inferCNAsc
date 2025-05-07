import json
import numpy as np
import pandas as pd
import requests_cache

def chunks(lst, size):
    """Yield successive chunks from a list."""
    for i in range(0, len(lst), size):
        yield lst[i:i + size]


def find_gen_reg(adata, cache_path='ensembl_cache', expire_after=43200):
    """
    Annotate adata.var with genomic regions (chrom, start, end) via Ensembl REST API.
    """
    session = requests_cache.CachedSession(cache_path, backend='sqlite', expire_after=expire_after)
    extracted = []
    for name in adata.var_names:
        if name.startswith('ENSG'):
            extracted.append(name)
        else:
            # extract ENSG ID in parentheses
            start = name.find('(')
            end = name.find(')')
            gene_id = None
            if start != -1 and end != -1 and name[start+1:start+5] == 'ENSG':
                gene_id = name[start+1:end]
            extracted.append(gene_id)
    results = {}
    for batch in chunks(extracted, 1000):
        url = 'https://rest.ensembl.org/lookup/id'
        headers = {'Content-Type': 'application/json'}
        payload = {'ids': batch}
        resp = session.post(url, headers=headers, data=json.dumps(payload))
        if resp.status_code == 200:
            data = resp.json()
            for gid in batch:
                info = data.get(gid, {}) if gid else {}
                results[gid] = {
                    'chrom': info.get('seq_region_name'),
                    'start': info.get('start'),
                    'end': info.get('end')
                }
        else:
            for gid in batch:
                results[gid] = {'chrom': None, 'start': None, 'end': None}
    # map to adata.var
    var = adata.var
    genelist = []
    chroms, starts, ends = [], [], []
    for gid in extracted:
        genelist.append(gid)
        info = results.get(gid, {})
        chroms.append(str(info.get('chrom')) if info.get('chrom') is not None else None)
        starts.append(int(info.get('start')) if info.get('start') is not None else None)
        ends.append(int(info.get('end')) if info.get('end') is not None else None)
    var['gene'] = genelist
    var['chrom'] = chroms
    var['start'] = starts
    var['end'] = ends
    return adata


def smooth_expression(adata, window_size=50):
    """
    Smooth expression by moving average along sorted genes.
    Returns smoothed numpy array of shape (cells, genes).
    """
    if 'counts' in adata.layers:
        raw = adata.layers['counts']
    else:
        raw = adata.X
    # filter and sort
    mask = adata.var['chrom'].notna().values
    idx = np.where(mask)[0]
    adata_sub = adata[:, idx]
    sorted_idx = np.argsort(list(zip(adata_sub.var['chrom'], adata_sub.var['start'])))
    expr = raw[:, idx][:, sorted_idx]
    n_cells, n_genes = expr.shape
    sm = np.zeros((n_cells, n_genes), dtype=float)
    for i in range(n_genes):
        start = max(0, i - window_size//2)
        end = min(n_genes, i + window_size//2 + 1)
        sm[:, i] = expr[:, start:end].mean(axis=1)
    return sm


def find_cnas(smoothed_expression, z_score_threshold=2.0):
    """
    Identify gains and losses based on z-score threshold.
    Returns boolean masks (gains, losses).
    """
    median = np.median(smoothed_expression, axis=0)
    std = np.std(smoothed_expression, axis=0)
    std = np.where(std < 1e-6, 1e-6, std)
    z = (smoothed_expression - median) / std
    gains = z > z_score_threshold
    losses = z < -z_score_threshold
    return gains, losses


def group_cnas(adata, gains, losses):
    """
    Group CNAs by chromosome and contiguous regions.
    Returns gain_groups, loss_groups, gains_df, losses_df.
    """
    def geno_to_df(geno, label):
        cells, genes = np.where(geno)
        df = pd.DataFrame({
            'cell': cells,
            'gene': genes,
            'chrom': adata.var.iloc[genes]['chrom'].values,
            'start': adata.var.iloc[genes]['start'].values,
            'end': adata.var.iloc[genes]['end'].values,
            'label': label
        })
        return df
    gains_df = geno_to_df(gains, 'gain')
    losses_df = geno_to_df(losses, 'loss')
    def group_rows(df):
        blocks = []
        last = None
        for _, row in df.sort_values(['cell','chrom','start']).iterrows():
            if last is None or row['chrom']!=last['chrom'] or row['start']!=last['end'] or row['cell']!=last['cell']:
                blocks.append([row.to_dict()])
            else:
                blocks[-1].append(row.to_dict())
            last = row
        return blocks
    gain_groups = group_rows(gains_df)
    loss_groups = group_rows(losses_df)
    return gain_groups, loss_groups, gains_df, losses_df