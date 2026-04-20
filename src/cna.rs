use ndarray::{Array2, Axis, Zip};

/// Identifies copy number gains and losses from a smoothed expression matrix.
///
/// For each gene (column), computes the column mean across all cells as the
/// baseline and the column standard deviation as the spread. A cell is flagged
/// as a gain if its z-score exceeds `z_score_threshold` and a loss if it falls
/// below `-z_score_threshold`. Genes with zero standard deviation are skipped.
///
/// Columns are processed in parallel with rayon.
///
/// Returns two boolean arrays of the same shape as `smoothed`: `gains` and
/// `losses`. A cell may be flagged in at most one of the two arrays.
pub fn find_cnas(smoothed: &Array2<f64>, z_score_threshold: f64) -> (Array2<bool>, Array2<bool>) {
    let shape = smoothed.dim();
    let mut gains = Array2::<bool>::from_elem(shape, false);
    let mut losses = Array2::<bool>::from_elem(shape, false);

    Zip::from(smoothed.axis_iter(Axis(1)))
        .and(gains.axis_iter_mut(Axis(1)))
        .and(losses.axis_iter_mut(Axis(1)))
        .par_for_each(|col, mut gcol, mut lcol| {
            let mean = col.mean().unwrap_or(0.0);
            let std = col.std(0.0);
            if std == 0.0 {
                return;
            }
            for (i, &v) in col.iter().enumerate() {
                let z = (v - mean) / std;
                gcol[i] = z > z_score_threshold;
                lcol[i] = z < -z_score_threshold;
            }
        });

    (gains, losses)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn gains_and_losses_shape_matches_input() {
        let smoothed = array![[1.0_f64, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0],];
        let (gains, losses) = find_cnas(&smoothed, 1.0);
        assert_eq!(gains.dim(), (3, 3));
        assert_eq!(losses.dim(), (3, 3));
    }

    #[test]
    fn high_value_cell_is_gain() {
        let smoothed = array![[1.0_f64, 5.0], [1.0, 5.0], [100.0, 5.0],];
        let (gains, losses) = find_cnas(&smoothed, 1.0);
        assert!(gains[[2, 0]], "cell 2 gene 0 should be a gain");
        assert!(!losses[[2, 0]], "cell 2 gene 0 should not be a loss");
    }

    #[test]
    fn low_value_cell_is_loss() {
        let smoothed = array![[100.0_f64, 5.0], [100.0, 5.0], [1.0, 5.0],];
        let (gains, losses) = find_cnas(&smoothed, 1.0);
        assert!(losses[[2, 0]], "cell 2 gene 0 should be a loss");
        assert!(!gains[[2, 0]], "cell 2 gene 0 should not be a gain");
    }

    #[test]
    fn zero_variance_column_produces_no_cna() {
        let smoothed = array![[5.0_f64, 1.0], [5.0, 2.0], [5.0, 3.0],];
        let (gains, losses) = find_cnas(&smoothed, 0.5);
        for c in 0..3 {
            assert!(!gains[[c, 0]]);
            assert!(!losses[[c, 0]]);
        }
    }
}
