//! PyO3 extension module for SAMRust (M3–M9).

mod alignment;
mod error;
mod segment;
mod variant;

use pyo3::prelude::*;
use samrust_core::VERSION;

use alignment::{PyAlignmentFile, PyBatchIterator, PyFetchIterator};
use segment::PyAlignedSegment;
use variant::{PyVariantFile, PyVariantHeader, PyVariantRecord};

/// Return the SAMRust package version string.
#[pyfunction]
fn version() -> &'static str {
    VERSION
}

/// Native extension imported as `samrust._samrust`.
#[pymodule]
fn _samrust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", VERSION)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_class::<PyAlignmentFile>()?;
    m.add_class::<PyAlignedSegment>()?;
    m.add_class::<PyFetchIterator>()?;
    m.add_class::<PyBatchIterator>()?;
    m.add_class::<PyVariantFile>()?;
    m.add_class::<PyVariantRecord>()?;
    m.add_class::<PyVariantHeader>()?;
    Ok(())
}
