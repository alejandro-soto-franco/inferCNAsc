#[derive(Debug, thiserror::Error)]
pub enum InferError {
    #[error("shape mismatch: expression has {expr} genes but chroms has {chroms}")]
    ShapeMismatch { expr: usize, chroms: usize },
    #[error("empty expression matrix")]
    EmptyMatrix,
}
