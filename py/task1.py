def chunks(lst, size):  # split data into chunks to make code run faster
  for i in range(0, len(lst), size):
    yield lst[i:i + size]
#%%
# find genomic region

def find_gen_reg(adata, gene_start_char = '(', gene_end_char = ')'):
  session = requests_cache.CachedSession('ensembl_cache', backend = 'sqlite', expire_after = 43200)   # make cache to make re-running code more efficient

  extracted_genes = []
  for i in adata.var.index.tolist():
    if i.startswith('ENSG'):
      gene_id = i

    else:
      gene_start = i.find(gene_start_char)
      gene_end = i.find(gene_end_char)

      if gene_start != -1 and gene_end != -1 and i[gene_start + 1: gene_start + 5] == 'ENSG':
        gene_id = i[gene_start + 1: gene_end]

      else:
        gene_id = None

    extracted_genes.append(gene_id)

  results_dict = {}
  for batch in chunks(extracted_genes, 1000):
    url = "https://rest.ensembl.org/lookup/id"  # limited to extent of Ensembl database
    headers = {'Content-Type': 'application/json'}
    payload = {'ids': batch}

    response = session.post(url, headers = headers, data = json.dumps(payload))

    if response.status_code == 200:
      data = response.json()

      for gene_id in batch:
        if gene_id != None:
          gene_info = data.get(gene_id, {})

          if gene_info:
            results_dict[gene_id] = {'chrom': gene_info.get('seq_region_name'), 'start': gene_info.get('start'), 'end': gene_info.get('end')}

          else:
            print(f"Gene ID {gene_id} not found in Ensembl data.")
            results_dict[gene_id] = {'chrom': None, 'start': None, 'end': None}

        else:
          print(f"Gene ID {gene_id} not found in Ensembl data.")
          results_dict[gene_id] = {'chrom': None, 'start': None, 'end': None}

  map_genes = {'gene': [], 'chrom': [], 'start': [], 'end': []}
  for gene_id in extracted_genes:
    map_genes['gene'].append(gene_id)

    if gene_id and gene_id in results_dict:
      map_genes['chrom'].append(str(results_dict[gene_id]['chrom']))

      if results_dict[gene_id]['start'] != None and results_dict[gene_id]['end'] != None:
        map_genes['start'].append(int(results_dict[gene_id]['start']))
        map_genes['end'].append(int(results_dict[gene_id]['end']))

      else:
        map_genes['start'].append(None)
        map_genes['end'].append(None)

    else:
      map_genes['chrom'].append(None)
      map_genes['start'].append(None)
      map_genes['end'].append(None)

  adata.var['gene'] = map_genes['gene']
  adata.var['chrom'] = map_genes['chrom']
  adata.var['start'] = map_genes['start']
  adata.var['end'] = map_genes['end']
#%%
# average/smooth gene expression to allow detection of CNAs

def smooth_expression(adata, window_size = 50):
  if 'counts' in adata.layers:
    adata.raw = adata.copy()
    adata.X = adata.layers['counts'].copy()

  adata = adata[:, np.where(adata.var['chrom'].notna())[0]].copy()  # remove genes with missing chromosome
  adata.var = adata.var.sort_values(by=['chrom', 'start'])  # sort genes by chromosome and position

  smoothed = []
  for i in range(adata.var.shape[0]):
    window_start = max(0, i - window_size // 2)
    window_end = min(adata.var.shape[0], i + window_size // 2 + 1)

    smoothed_expr = np.mean(adata.raw.X[:, window_start: window_end].toarray(), axis=1)
    smoothed.append(smoothed_expr)

  smoothed_array = np.array(smoothed).reshape(adata.var.shape[0], -1)
  smoothed_expression = smoothed_array.T.copy()  # add smoothed expression to adata object
  return smoothed_expression
#%%
# ID CNAs

def find_cnas(smoothed_expression, z_score_threshold = 1):
  normal_expression = np.median(smoothed_expression, axis=0)
  z_scores = (smoothed_expression - normal_expression) / np.std(smoothed_expression, axis = 0)

  gains = z_scores > z_score_threshold
  losses = z_scores < -z_score_threshold

  return gains, losses
#%%
# group CNAs by genomic region

def group_cnas(adata, gains, losses):
  def geno_to_df(genotype, label):
    indices = np.where(genotype)
    cell_idx, gene_idx = indices
    return pd.DataFrame({'cell': cell_idx, 'gene': gene_idx, 'chrom': adata.var.iloc[gene_idx]['chrom'].values, 'start': adata.var.iloc[gene_idx]['start'].values, 'end': adata.var.iloc[gene_idx]['end'].values, 'label': label})

  gains_df = geno_to_df(gains, 'gain')
  losses_df = geno_to_df(losses, 'loss')

  # grouping by chromosome (you can group more finely if needed)
  # merge consecutive gains/losses
  def group_rows(df):
    grouped = []
    last_row = None

    for idx, row in df.iterrows():
      if last_row is None or row['chrom'] != last_row['chrom'] or row['start'] != last_row['end']:
        grouped.append([row])

      else:
        grouped[-1].append(row)

      last_row = row

    return grouped

  gain_groups = group_rows(gains_df)
  loss_groups = group_rows(losses_df)

  return gain_groups, loss_groups, gains_df, losses_df
#%%
# heatmap smoothed expression

def smooth_heatmap(smoothed_expression):
  plt.figure(figsize = (10, 6))
  plt.imshow(smoothed_expression, cmap = 'coolwarm', aspect = 'auto')
  plt.colorbar(label = 'Smoothed Expression')
  plt.title('Smoothed Gene Expression Across Cells')
  plt.ylabel('Genes')
  plt.xlabel('Cells')
  plt.show()
#%%
# sort by chromosome for visualization

def chrom_sort_key(chrom):
    chrom = str(chrom).replace('chr', '')
    return int(chrom) if chrom.isdigit() else {'X': 23, 'Y': 24}.get(chrom, 25)
#%%
# visualize the detected CNAs (gains/losses) along chromosomes in scatter plot

def plot_cnas(gains_df, losses_df):
    plt.figure(figsize=(12, 6))

    combined = pd.concat([gains_df.assign(type = 'Gain'), losses_df.assign(type = 'Loss')])
    combined['chrom_sort'] = combined['chrom'].apply(chrom_sort_key)
    combined = combined.sort_values(by = 'chrom_sort')

    for ctype, color in [('Gain', 'red'), ('Loss', 'blue')]:
        subset = combined[combined['type'] == ctype]
        plt.scatter(subset['chrom'], subset['start'], c = color, label = ctype, alpha = 0.6)

    plt.xlabel('Chromosome')
    plt.ylabel('Start Position')
    plt.title('Detected CNAs: Gains and Losses')
    plt.xticks(rotation = 90)  # Rotate x-axis labels 90 degrees counterclockwise
    plt.legend()
    plt.tight_layout()  # Prevent label overlap
    plt.show()
#%%
# plot how gene expression changes along a chromosome in line plot

def gene_express_line(smoothed_expression, adata):
    plt.figure(figsize = (12, 6))
    avg_expr = smoothed_expression.mean(axis=0)

    # Sort chromosomes using chrom_sort_key
    sorted_var = adata.var.sort_values(['chrom', 'start'], key = lambda col: col.map(chrom_sort_key) if col.name == 'chrom' else col)
    chroms = sorted(sorted_var['chrom'].unique(), key = chrom_sort_key)

    genome_pos = []
    offset = 0
    tick_positions = []
    tick_labels = []

    for chrom in chroms:
        chrom_genes = sorted_var[sorted_var['chrom'] == chrom]
        positions = chrom_genes['start'].values + offset
        genome_pos.extend(positions)

        tick_positions.append(offset + (positions[-1] - positions[0]) / 2)
        tick_labels.append(chrom)

        offset = positions[-1] + 1_000_000  # gap between chromosomes

    sorted_idx = sorted_var.index
    sorted_expr = avg_expr[adata.var.index.get_indexer(sorted_idx)]

    plt.plot(genome_pos, sorted_expr)
    plt.xticks(tick_positions, tick_labels, rotation = 90)
    plt.xlabel('Chromosome')
    plt.ylabel('Avg Smoothed Expression')
    plt.title('Gene Expression Along Genome')
    plt.tight_layout()
    plt.show()
#%%
# assign CNAs to cell types

def assign_cnas_to_cells(gains, losses, adata, min_region_size=3):
    cnas_by_cell = []

    gene_chrom = adata.var['chrom'].values
    gene_start = adata.var['start'].values
    gene_names = adata.var_names

    for cell_idx in range(gains.shape[0]):
        for cna_type, matrix in [('gain', gains), ('loss', losses)]:
            start = None
            current_chr = None

            for gene_idx in range(matrix.shape[1]):
                if matrix[cell_idx, gene_idx]:
                    chrom = gene_chrom[gene_idx]
                    if start is None:
                        start = gene_idx
                        current_chr = chrom
                    elif chrom != current_chr:
                        if gene_idx - start >= min_region_size:
                            cnas_by_cell.append({
                                'cell': cell_idx,
                                'chrom': current_chr,
                                'start_gene': gene_names[start],
                                'end_gene': gene_names[gene_idx-1],
                                'start_pos': gene_start[start],
                                'end_pos': gene_start[gene_idx-1],
                                'label': cna_type
                            })
                        start = gene_idx
                        current_chr = chrom
                else:
                    if start is not None and gene_idx - start >= min_region_size:
                        cnas_by_cell.append({
                            'cell': cell_idx,
                            'chrom': current_chr,
                            'start_gene': gene_names[start],
                            'end_gene': gene_names[gene_idx-1],
                            'start_pos': gene_start[start],
                            'end_pos': gene_start[gene_idx-1],
                            'label': cna_type
                        })
                    start = None
                    current_chr = None

    return pd.DataFrame(cnas_by_cell)

def assign_cnas_to_clusters(cna_df, cluster_labels):
    cna_df = cna_df.copy()
    cna_df['cluster'] = cna_df['cell'].map(cluster_labels)

    grouped = cna_df.groupby(['cluster', 'chrom', 'start_pos', 'end_pos', 'label']).size().reset_index(name='count')
    return grouped
#%%
# visualize CNAs by cell type

def cna_matrix(cna_df, num_cells):
    # Each unique CNA region becomes a column
    cna_df['region_id'] = cna_df.apply(lambda row: f"{row['chrom']}:{row['start_pos']}-{row['end_pos']} ({row['label']})", axis=1)
    unique_regions = sorted(cna_df['region_id'].unique())

    region_index = {r: i for i, r in enumerate(unique_regions)}
    matrix = np.zeros((num_cells, len(unique_regions)))

    for _, row in cna_df.iterrows():
        i = row['cell']
        j = region_index[row['region_id']]
        matrix[i, j] = 1 if row['label'] == 'gain' else -1

    return matrix, unique_regions

def plot_cna_matrix(matrix, region_labels):
    plt.figure(figsize=(16, 8))
    plt.imshow(matrix, aspect = 'auto', cmap = 'bwr', vmin = -1, vmax = 1)
    plt.colorbar(label = 'CNA Status\n(1 = gain, -1 = loss)')
    plt.xticks(ticks = np.arange(len(region_labels)), labels = region_labels, rotation = 90, fontsize = 6)
    plt.yticks([])
    plt.xlabel('Genomic Region')
    plt.ylabel('Cells')
    plt.title('Copy Number Alterations per Cell')
    plt.tight_layout()
    plt.show()
#%%
# evaluate accuracy in comparison to known CNAs

def parse_simulated(sim_cna_col):
  sim_list = []
  for entry in sim_cna_col:
    if pd.isna(entry):
        continue

    parts = entry.split()
    chrom, coords = parts[0].split(':')
    start, end = map(int, coords.split('-'))
    cn = int(parts[-1])  # e.g., CN 0 => 0
    label = 'loss' if cn < 2 else 'gain' if cn > 2 else 'neutral'
    sim_list.append({'chrom': chrom, 'start': start, 'end': end, 'label': label})
  return pd.DataFrame(sim_list)

def region_overlap(start1, end1, start2, end2):
  return not (end1 < start2 or end2 < start1)

def evaluate_accuracy(simulated_df, inferred_df):
  true_positives = 0
  false_positives = 0
  matched = set()

  for i, inferred in inferred_df.iterrows():
    hit = False
    for j, sim in simulated_df.iterrows():
        if sim['chrom'] == str(inferred['chrom']) and region_overlap(sim['start'], sim['end'], inferred['start'], inferred['end']):
            if sim['label'] == inferred['label']:
                true_positives += 1
                matched.add(j)
                hit = True
                break
    if not hit:
        false_positives += 1

  false_negatives = len(simulated_df) - len(matched)

  print(f"True Positives: {true_positives}")
  print(f"False Positives: {false_positives}")
  print(f"False Negatives: {false_negatives}")

  precision = true_positives / (true_positives + false_positives + 1e-6)
  recall = true_positives / (true_positives + false_negatives + 1e-6)
  f1 = 2 * precision * recall / (precision + recall + 1e-6)

  print(f"Precision: {precision:.2f}")
  print(f"Recall: {recall:.2f}")
  print(f"F1 Score: {f1:.2f}")