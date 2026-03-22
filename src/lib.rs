pub mod error;
pub mod smooth;
pub mod cna;
pub mod assign;

pub use error::InferError;
pub use smooth::smooth_expression;
pub use cna::find_cnas;
pub use assign::{assign_cnas_to_cells, CnaRecord};

#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use numpy::{PyArray2, PyReadonlyArray1, PyReadonlyArray2, IntoPyArray};

/// Python extension module. Only compiled when the `python` feature is active.
#[cfg(feature = "python")]
#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_smooth_expression, m)?)?;
    m.add_function(wrap_pyfunction!(py_find_cnas, m)?)?;
    m.add_function(wrap_pyfunction!(py_assign_cnas_to_cells, m)?)?;
    Ok(())
}

/// Smooth expression matrix with a chromosome-aware sliding window.
#[cfg(feature = "python")]
#[pyfunction(name = "smooth_expression")]
fn py_smooth_expression<'py>(
    py: Python<'py>,
    expression: PyReadonlyArray2<'py, f64>,
    chroms: Vec<String>,
    window_size: usize,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let expr = expression.as_array().to_owned();
    let chrom_refs: Vec<&str> = chroms.iter().map(String::as_str).collect();
    let result = smooth_expression(&expr, &chrom_refs, window_size)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
    Ok(result.into_pyarray_bound(py))
}

/// Identify copy number gains and losses via z-score thresholding.
///
/// Returns two u8 arrays (0/1) representing gains and losses. On the Python
/// side these can be cast to bool with `.astype(bool)`. The u8 intermediate
/// is required because numpy 0.22 does not implement `IntoPyArray` for
/// `Array2<bool>` directly.
#[cfg(feature = "python")]
#[pyfunction(name = "find_cnas")]
fn py_find_cnas<'py>(
    py: Python<'py>,
    smoothed: PyReadonlyArray2<'py, f64>,
    z_score_threshold: f64,
) -> PyResult<(Bound<'py, PyArray2<u8>>, Bound<'py, PyArray2<u8>>)> {
    let arr = smoothed.as_array().to_owned();
    let (gains, losses) = find_cnas(&arr, z_score_threshold);
    let gains_u8 = gains.mapv(|x| x as u8);
    let losses_u8 = losses.mapv(|x| x as u8);
    Ok((gains_u8.into_pyarray_bound(py), losses_u8.into_pyarray_bound(py)))
}

/// Assign CNA regions to individual cells.
///
/// Accepts u8 arrays (0/1) for gains and losses, converting them to bool
/// internally before passing to the Rust core.
#[cfg(feature = "python")]
#[pyfunction(name = "assign_cnas_to_cells")]
fn py_assign_cnas_to_cells<'py>(
    py: Python<'py>,
    gains: PyReadonlyArray2<'py, u8>,
    losses: PyReadonlyArray2<'py, u8>,
    chroms: Vec<String>,
    starts: PyReadonlyArray1<'py, i64>,
    ends: PyReadonlyArray1<'py, i64>,
    gene_names: Vec<String>,
    min_region_size: usize,
) -> PyResult<Bound<'py, pyo3::types::PyList>> {
    use pyo3::types::{PyDict, PyList};

    let gains_bool = gains.as_array().mapv(|x| x != 0);
    let losses_bool = losses.as_array().mapv(|x| x != 0);
    let chrom_refs: Vec<&str> = chroms.iter().map(String::as_str).collect();
    let name_refs: Vec<&str> = gene_names.iter().map(String::as_str).collect();
    let starts_slice = starts.as_slice()?;
    let ends_slice = ends.as_slice()?;

    let records = assign_cnas_to_cells(
        &gains_bool, &losses_bool,
        &chrom_refs, starts_slice, ends_slice,
        &name_refs, min_region_size,
    );

    let list = PyList::empty_bound(py);
    for r in &records {
        let dict = PyDict::new_bound(py);
        dict.set_item("cell", r.cell)?;
        dict.set_item("chrom", &r.chrom)?;
        dict.set_item("start_gene", &r.start_gene)?;
        dict.set_item("end_gene", &r.end_gene)?;
        dict.set_item("start_pos", r.start_pos)?;
        dict.set_item("end_pos", r.end_pos)?;
        dict.set_item("label", &r.label)?;
        list.append(dict)?;
    }
    Ok(list)
}
