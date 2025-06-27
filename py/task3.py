find_gen_reg(data1norm)   # runs in ~1 min

d1_smoothed_expression = smooth_expression(data1norm)  # runs in ~30 mins

d1_gains, d1_losses = find_cnas(d1_smoothed_expression)  # runs in ~10 secs

d1_gain_groups, d1_loss_groups, d1_gain_df, d1_loss_df = group_cnas(data1norm, d1_gains, d1_losses)  # runs in ~10 mins

smooth_heatmap(d1_smoothed_expression)  # runs in ~1 min

plot_cnas(d1_gain_df, d1_loss_df)  # runs in ~2 mins

gene_express_line(d1_smoothed_expression, data1norm)  # runs in ~1 sec

cna_df3 = assign_cnas_to_cells(d1_gains, d1_losses, data1norm)

grouped3 = assign_cnas_to_clusters(cna_df3, data1norm.obs['cell_type'])

matrix3, uni_reg3 = cna_matrix(cna_df3, data1norm.n_obs)

plot_cna_matrix(matrix3, uni_reg3)