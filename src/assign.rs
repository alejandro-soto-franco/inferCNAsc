use ndarray::Array2;

#[derive(Debug, Clone)]
pub struct CnaRecord {
    pub cell: usize,
    pub chrom: String,
    pub start_gene: String,
    pub end_gene: String,
    pub start_pos: i64,
    pub end_pos: i64,
    pub label: String,
}

pub fn assign_cnas_to_cells(
    _gains: &Array2<bool>,
    _losses: &Array2<bool>,
    _chroms: &[&str],
    _starts: &[i64],
    _ends: &[i64],
    _gene_names: &[&str],
    _min_region_size: usize,
) -> Vec<CnaRecord> {
    unimplemented!()
}
