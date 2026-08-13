//! Map [`samrust_core::SamRustError`] to Python exceptions.

use pyo3::exceptions::{PyIOError, PyRuntimeError, PyValueError};
use pyo3::PyErr;
use samrust_core::SamRustError;

pub fn to_pyerr(err: SamRustError) -> PyErr {
    match err {
        SamRustError::InvalidArgument(msg) => PyValueError::new_err(msg),
        SamRustError::Io(e) => PyIOError::new_err(e.to_string()),
        SamRustError::MissingIndex(path) => {
            PyValueError::new_err(format!("could not find index for {path}"))
        }
        SamRustError::NotImplemented(msg) => PyRuntimeError::new_err(msg),
    }
}
