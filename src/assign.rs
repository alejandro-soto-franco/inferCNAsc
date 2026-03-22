use ndarray::Array2;

/// A single detected copy number alteration region for one cell.
#[derive(Debug, Clone)]
pub struct CnaRecord {
    /// Zero-based index of the cell.
    pub cell: usize,
    /// Chromosome name (e.g. "1", "X").
    pub chrom: String,
    /// Gene identifier at the start of the region.
    pub start_gene: String,
    /// Gene identifier at the end of the region (inclusive).
    pub end_gene: String,
    /// Genomic start coordinate of the first gene in the region.
    pub start_pos: i64,
    /// Genomic end coordinate of the last gene in the region.
    pub end_pos: i64,
    /// Either "gain" or "loss".
    pub label: String,
}

/// Assigns detected CNAs to individual cells using a run-length merge strategy.
///
/// For each cell and each CNA type (gain, loss), the function scans the boolean
/// array gene by gene. Consecutive flagged genes on the same chromosome form a
/// run. A run is emitted as a `CnaRecord` when it ends (by an unflagged gene or
/// a chromosome boundary) and its length meets `min_region_size`. This step
/// replaces the separate `group_cnas` grouping pass from the Python prototype.
///
/// The inputs must be pre-sorted by `(chrom, start)` (enforced by `fit()`).
pub fn assign_cnas_to_cells(
    gains: &Array2<bool>,
    losses: &Array2<bool>,
    chroms: &[&str],
    starts: &[i64],
    ends: &[i64],
    gene_names: &[&str],
    min_region_size: usize,
) -> Vec<CnaRecord> {
    let n_cells = gains.nrows();
    let n_genes = gains.ncols();
    let mut records = Vec::new();

    for cell in 0..n_cells {
        for (label, matrix) in [("gain", gains), ("loss", losses)] {
            let mut run_start: Option<usize> = None;

            for g in 0..=n_genes {
                // A boundary occurs when the run must end: past the last gene,
                // current gene is unflagged, or chromosome changes mid-run.
                let chrom_break = run_start.map_or(false, |s| {
                    g < n_genes && chroms[g] != chroms[s]
                });
                let boundary = g == n_genes
                    || (g < n_genes && !matrix[[cell, g]])
                    || chrom_break;

                if boundary {
                    // Flush any open run.
                    if let Some(s) = run_start {
                        if g - s >= min_region_size {
                            records.push(CnaRecord {
                                cell,
                                chrom: chroms[s].to_string(),
                                start_gene: gene_names[s].to_string(),
                                end_gene: gene_names[g - 1].to_string(),
                                start_pos: starts[s],
                                end_pos: ends[g - 1],
                                label: label.to_string(),
                            });
                        }
                        run_start = None;
                    }
                    // Start a new run if current gene is flagged (chrom break case).
                    if g < n_genes && matrix[[cell, g]] {
                        run_start = Some(g);
                    }
                } else if run_start.is_none() && g < n_genes && matrix[[cell, g]] {
                    // Start a run at the first flagged gene.
                    run_start = Some(g);
                }
            }
        }
    }

    records
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    fn make_bool(rows: Vec<Vec<bool>>) -> Array2<bool> {
        let n_rows = rows.len();
        let n_cols = rows[0].len();
        let flat: Vec<bool> = rows.into_iter().flatten().collect();
        Array2::from_shape_vec((n_rows, n_cols), flat).unwrap()
    }

    #[test]
    fn record_contains_all_seven_fields() {
        // 1 cell, 5 genes on chr1, all gains flagged.
        let gains = make_bool(vec![vec![true, true, true, true, true]]);
        let losses = make_bool(vec![vec![false, false, false, false, false]]);
        let chroms = vec!["1", "1", "1", "1", "1"];
        let starts = vec![100i64, 200, 300, 400, 500];
        let ends   = vec![199i64, 299, 399, 499, 599];
        let names  = vec!["G1", "G2", "G3", "G4", "G5"];

        let records = assign_cnas_to_cells(
            &gains, &losses, &chroms, &starts, &ends, &names, 3,
        );

        assert_eq!(records.len(), 1);
        let r = &records[0];
        assert_eq!(r.cell, 0);
        assert_eq!(r.chrom, "1");
        assert_eq!(r.start_gene, "G1");
        assert_eq!(r.end_gene, "G5");
        assert_eq!(r.start_pos, 100);
        assert_eq!(r.end_pos, 599);
        assert_eq!(r.label, "gain");
    }

    #[test]
    fn run_shorter_than_min_region_size_is_dropped() {
        // 1 cell, 5 genes, only 2 flagged (below min_region_size=3).
        let gains = make_bool(vec![vec![true, true, false, false, false]]);
        let losses = make_bool(vec![vec![false; 5]]);
        let chroms = vec!["1"; 5];
        let starts: Vec<i64> = (0..5).map(|i| i * 100).collect();
        let ends: Vec<i64>   = (0..5).map(|i| i * 100 + 99).collect();
        let names: Vec<&str> = (0..5).map(|i| ["G1","G2","G3","G4","G5"][i]).collect();

        let records = assign_cnas_to_cells(
            &gains, &losses, &chroms, &starts, &ends, &names, 3,
        );
        assert!(records.is_empty(), "run of 2 should be dropped");
    }

    #[test]
    fn chromosome_boundary_splits_runs() {
        // 1 cell, 4 genes: chr1 (2) then chr2 (2), all gains flagged.
        // min_region_size=2 so each chromosome run should emit a record.
        let gains = make_bool(vec![vec![true, true, true, true]]);
        let losses = make_bool(vec![vec![false; 4]]);
        let chroms = vec!["1", "1", "2", "2"];
        let starts: Vec<i64> = vec![100, 200, 100, 200];
        let ends: Vec<i64>   = vec![199, 299, 199, 299];
        let names  = vec!["G1", "G2", "G3", "G4"];

        let records = assign_cnas_to_cells(
            &gains, &losses, &chroms, &starts, &ends, &names, 2,
        );
        assert_eq!(records.len(), 2, "one record per chromosome");
        assert_eq!(records[0].chrom, "1");
        assert_eq!(records[1].chrom, "2");
    }

    #[test]
    fn losses_labeled_correctly() {
        let gains  = make_bool(vec![vec![false; 4]]);
        let losses = make_bool(vec![vec![true, true, true, true]]);
        let chroms = vec!["1"; 4];
        let starts: Vec<i64> = vec![0, 100, 200, 300];
        let ends: Vec<i64>   = vec![99, 199, 299, 399];
        let names  = vec!["G1", "G2", "G3", "G4"];

        let records = assign_cnas_to_cells(
            &gains, &losses, &chroms, &starts, &ends, &names, 2,
        );
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].label, "loss");
    }
}
