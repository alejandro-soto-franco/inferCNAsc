/// Errors returned by the infercnasc Rust core.
#[derive(Debug, thiserror::Error)]
pub enum InferError {
    /// The number of genes in the expression matrix does not match the length
    /// of the `chroms` slice passed to `smooth_expression`.
    #[error("shape mismatch: expression has {expr} genes but chroms has {chroms}")]
    ShapeMismatch { expr: usize, chroms: usize },

    /// The expression matrix has zero rows or zero columns.
    #[error("empty expression matrix")]
    EmptyMatrix,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shape_mismatch_message_contains_counts() {
        let e = InferError::ShapeMismatch { expr: 100, chroms: 50 };
        let msg = format!("{e}");
        assert!(msg.contains("100"), "message: {msg}");
        assert!(msg.contains("50"), "message: {msg}");
    }

    #[test]
    fn empty_matrix_message_is_nonempty() {
        let e = InferError::EmptyMatrix;
        assert!(!format!("{e}").is_empty());
    }
}
