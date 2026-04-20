use crate::error::InferError;
use ndarray::{Array2, Axis, s};
use rayon::prelude::*;

/// Applies a chromosome-aware sliding-window mean to an expression matrix.
///
/// For each gene position `g`, the average is computed over genes in the window
/// `[g - window_size/2, g + window_size/2]` clamped to the boundaries of the
/// chromosome containing gene `g`. Genes on adjacent chromosomes never
/// contribute to the window.
///
/// The input matrix must be pre-sorted by `(chrom, start)` with columns
/// matching the order of entries in `chroms`. This is enforced by the Python
/// `fit()` method before calling into Rust.
///
/// # Errors
/// Returns `InferError::EmptyMatrix` if either dimension is zero.
/// Returns `InferError::ShapeMismatch` if `chroms.len() != n_genes`.
pub fn smooth_expression(
    expression: &Array2<f64>,
    chroms: &[&str],
    window_size: usize,
) -> Result<Array2<f64>, InferError> {
    let (n_cells, n_genes) = expression.dim();
    if n_cells == 0 || n_genes == 0 {
        return Err(InferError::EmptyMatrix);
    }
    if chroms.len() != n_genes {
        return Err(InferError::ShapeMismatch {
            expr: n_genes,
            chroms: chroms.len(),
        });
    }

    // Identify the start and end (exclusive) of each chromosome segment.
    // gene_seg_start[g] and gene_seg_end[g] are the inclusive start and
    // exclusive end of the chromosome run that contains gene g.
    let mut gene_seg_start = vec![0usize; n_genes];
    let mut gene_seg_end = vec![n_genes; n_genes];

    let mut seg_start = 0usize;
    let mut i = 1usize;
    while i <= n_genes {
        if i == n_genes || chroms[i] != chroms[seg_start] {
            for g in seg_start..i {
                gene_seg_start[g] = seg_start;
                gene_seg_end[g] = i;
            }
            seg_start = i;
        }
        i += 1;
    }

    let half = window_size / 2;

    // Compute one smoothed column per gene in parallel.
    let smoothed_cols: Vec<Vec<f64>> = (0..n_genes)
        .into_par_iter()
        .map(|g| {
            let seg_s = gene_seg_start[g];
            let seg_e = gene_seg_end[g];
            let win_s = g.saturating_sub(half).max(seg_s);
            let win_e = (g + half + 1).min(seg_e);
            expression
                .slice(s![.., win_s..win_e])
                .mean_axis(Axis(1))
                .expect("window is always non-empty")
                .to_vec()
        })
        .collect();

    // Assemble: smoothed_cols[gene][cell] -> output[cell][gene]
    let mut output = Array2::<f64>::zeros((n_cells, n_genes));
    for (g, col) in smoothed_cols.iter().enumerate() {
        for (c, &val) in col.iter().enumerate() {
            output[[c, g]] = val;
        }
    }

    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ndarray::array;

    #[test]
    fn basic_smoothing_shape_preserved() {
        let expr = array![[1.0_f64, 2.0, 3.0, 4.0, 5.0], [2.0, 3.0, 4.0, 5.0, 6.0],];
        let chroms = vec!["1", "1", "1", "1", "1"];
        let result = smooth_expression(&expr, &chroms, 3).unwrap();
        assert_eq!(result.dim(), (2, 5));
    }

    #[test]
    fn chromosome_boundary_resets_window() {
        // 2 cells, 4 genes: chr1 has genes 0-1 (values 10, 20),
        // chr2 has genes 2-3 (values 1, 2).
        // With window_size=4, gene 2 should only average over chr2 genes,
        // so result for cell 0, gene 2 should be (1+2)/2 = 1.5, not (20+1+2)/3.
        let expr = array![[10.0_f64, 20.0, 1.0, 2.0], [10.0, 20.0, 1.0, 2.0],];
        let chroms = vec!["1", "1", "2", "2"];
        let result = smooth_expression(&expr, &chroms, 4).unwrap();
        assert_abs_diff_eq!(result[[0, 2]], 1.5, epsilon = 1e-10);
        assert_abs_diff_eq!(result[[0, 3]], 1.5, epsilon = 1e-10);
    }

    #[test]
    fn empty_matrix_returns_error() {
        use ndarray::Array2;
        let expr = Array2::<f64>::zeros((0, 0));
        assert!(smooth_expression(&expr, &[], 3).is_err());
    }

    #[test]
    fn shape_mismatch_returns_error() {
        use ndarray::Array2;
        let expr = Array2::<f64>::zeros((3, 5));
        let chroms = vec!["1", "1"]; // length 2, not 5
        let result = smooth_expression(&expr, &chroms, 3);
        assert!(matches!(
            result,
            Err(crate::error::InferError::ShapeMismatch { .. })
        ));
    }

    #[test]
    fn single_gene_single_chrom_identity() {
        // A single gene has no neighbors; the window contains only itself.
        let expr = array![[7.0_f64], [3.0]];
        let chroms = vec!["1"];
        let result = smooth_expression(&expr, &chroms, 5).unwrap();
        assert_abs_diff_eq!(result[[0, 0]], 7.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result[[1, 0]], 3.0, epsilon = 1e-10);
    }
}
