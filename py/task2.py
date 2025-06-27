# downsampling

cna_benchnorm_downsampled1 = cna_benchnorm.copy()
sc.pp.downsample_counts(cna_benchnorm_downsampled1, total_counts = 10000000)

find_gen_reg(cna_benchnorm_downsampled1)   # runs in 2 mins

bench_down_smoothed_expression1 = smooth_expression(cna_benchnorm_downsampled1)  # runs in ~2 mins

bench_down_gains1, bench_down_losses1 = find_cnas(bench_down_smoothed_expression1)  # runs in ~5 secs

bench_down_gain_groups1, bench_down_loss_groups1, bench_down_gain_df1, bench_down_loss_df1 = group_cnas(cna_benchnorm_downsampled1, bench_down_gains1, bench_down_losses1)   # runs in > 20 hours

smooth_heatmap(bench_down_smoothed_expression1)   # runs in ~25 mins

plot_cnas(bench_down_gain_df1, bench_down_loss_df1)   # runs in ~1 min

gene_express_line(bench_down_smoothed_expression1, cna_benchnorm_downsampled1)   # runs in 1 sec

cna_df1 = assign_cnas_to_cells(bench_down_gains1, bench_down_losses1, cna_benchnorm_downsampled1)

grouped1 = assign_cnas_to_clusters(cna_df1, cna_benchnorm_downsampled1.obs['cell_type'])

matrix1, uni_reg1 = cna_matrix(cna_df1, cna_benchnorm_downsampled1.n_obs)

plot_cna_matrix(matrix1, uni_reg1)

sim_df1 = parse_simulated(cna_benchnorm_downsampled1.obs['simulated_cnvs'])

evaluate_accuracy(sim_df1, cna_df1)


# extreme downsampling

cna_benchnorm_downsampled2 = cna_benchnorm.copy()
sc.pp.downsample_counts(cna_benchnorm_downsampled2, total_counts = 10000)

find_gen_reg(cna_benchnorm_downsampled2)   # runs in 2 mins

bench_down_smoothed_expression2 = smooth_expression(cna_benchnorm_downsampled2)  # runs in ~10 secs

bench_down_gains2, bench_down_losses2 = find_cnas(bench_down_smoothed_expression2)  # runs in ~5 secs

bench_down_gain_groups2, bench_down_loss_groups2, bench_down_gain_df2, bench_down_loss_df2 = group_cnas(cna_benchnorm_downsampled2, bench_down_gains2, bench_down_losses2)   # runs in ~10 secs

smooth_heatmap(bench_down_smoothed_expression2)   # runs in ~1 min

plot_cnas(bench_down_gain_df2, bench_down_loss_df2)   # runs in ~5 secs

gene_express_line(bench_down_smoothed_expression2, cna_benchnorm_downsampled2)   # runs in 1 sec

cna_df2 = assign_cnas_to_cells(bench_down_gains2, bench_down_losses2, cna_benchnorm_downsampled2)   # runs in ~30 secs

grouped2 = assign_cnas_to_clusters(cna_df2, cna_benchnorm_downsampled2.obs['cell_type'])    # runs in ~ mins

matrix2, uni_reg2 = cna_matrix(cna_df2, cna_benchnorm_downsampled2.n_obs)

plot_cna_matrix(matrix2, uni_reg2)

sim_df2 = parse_simulated(cna_benchnorm_downsampled2.obs['simulated_cnvs'])

evaluate_accuracy(sim_df2, cna_df2)