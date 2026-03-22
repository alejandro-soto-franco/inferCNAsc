use ndarray::Array2;
use crate::error::InferError;

pub fn smooth_expression(
    _expression: &Array2<f64>,
    _chroms: &[&str],
    _window_size: usize,
) -> Result<Array2<f64>, InferError> {
    unimplemented!()
}
