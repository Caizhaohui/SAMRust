//! PyO3 `AlignedSegment` wrapping [`samrust_core::Record`].

use pyo3::prelude::*;
use pyo3::types::PyList;
use pyo3::IntoPyObjectExt;
use samrust_core::{Cigar, Record, TagValue};

use crate::error::to_pyerr;

/// pysam-compatible alignment record view.
#[pyclass(name = "AlignedSegment", module = "samrust._samrust")]
#[derive(Clone)]
pub struct PyAlignedSegment {
    inner: Record,
}

impl PyAlignedSegment {
    pub fn new(record: Record) -> Self {
        Self { inner: record }
    }
}

#[pymethods]
impl PyAlignedSegment {
    #[getter]
    fn query_name(&self) -> &str {
        self.inner.query_name()
    }

    #[getter]
    fn flag(&self) -> u16 {
        self.inner.flag()
    }

    #[getter]
    fn reference_id(&self) -> i32 {
        self.inner.reference_id()
    }

    #[getter]
    fn reference_name(&self) -> Option<&str> {
        self.inner.reference_name()
    }

    #[getter]
    fn reference_start(&self) -> i64 {
        self.inner.reference_start()
    }

    #[getter]
    fn mapping_quality(&self) -> u8 {
        self.inner.mapping_quality()
    }

    #[getter]
    fn cigarstring(&self) -> Option<String> {
        self.inner.cigarstring()
    }

    #[getter]
    fn cigar(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        cigar_to_py(py, self.inner.cigar())
    }

    #[getter]
    fn cigartuples(&self) -> Vec<(u8, u32)> {
        self.inner.cigar().cigartuples()
    }

    #[getter]
    fn query_sequence(&self) -> &str {
        self.inner.query_sequence()
    }

    #[getter]
    fn query_length(&self) -> usize {
        self.inner.query_length()
    }

    #[getter]
    fn query_qualities(&self) -> Vec<u8> {
        self.inner.query_qualities().to_vec()
    }

    #[getter]
    fn next_reference_id(&self) -> i32 {
        self.inner.mate_reference_id()
    }

    #[getter]
    fn next_reference_start(&self) -> i64 {
        self.inner.mate_reference_start()
    }

    #[getter]
    fn template_length(&self) -> i32 {
        self.inner.template_length()
    }

    #[getter]
    fn is_paired(&self) -> bool {
        self.inner.is_paired()
    }

    #[getter]
    fn is_proper_pair(&self) -> bool {
        self.inner.is_proper_pair()
    }

    #[getter]
    fn is_unmapped(&self) -> bool {
        self.inner.is_unmapped()
    }

    #[getter]
    fn mate_is_unmapped(&self) -> bool {
        self.inner.mate_is_unmapped()
    }

    #[getter]
    fn is_reverse(&self) -> bool {
        self.inner.is_reverse()
    }

    #[getter]
    fn mate_is_reverse(&self) -> bool {
        self.inner.mate_is_reverse()
    }

    #[getter]
    fn is_read1(&self) -> bool {
        self.inner.is_read1()
    }

    #[getter]
    fn is_read2(&self) -> bool {
        self.inner.is_read2()
    }

    #[getter]
    fn is_secondary(&self) -> bool {
        self.inner.is_secondary()
    }

    #[getter]
    fn is_qcfail(&self) -> bool {
        self.inner.is_qcfail()
    }

    #[getter]
    fn is_duplicate(&self) -> bool {
        self.inner.is_duplicate()
    }

    #[getter]
    fn is_supplementary(&self) -> bool {
        self.inner.is_supplementary()
    }

    fn has_tag(&self, tag: &str) -> bool {
        self.inner.tags().get(tag).is_some()
    }

    fn get_tag(&self, py: Python<'_>, tag: &str) -> PyResult<Py<PyAny>> {
        match self.inner.tags().get(tag) {
            Some(val) => tag_value_to_py(py, val),
            None => Err(pyo3::exceptions::PyKeyError::new_err(format!(
                "tag {tag} not present"
            ))),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "AlignedSegment(query_name={:?}, flag={})",
            self.inner.query_name(),
            self.inner.flag()
        )
    }
}

fn cigar_to_py(py: Python<'_>, cigar: &Cigar) -> PyResult<Py<PyAny>> {
    let list = PyList::empty(py);
    for (op, len) in cigar.cigartuples() {
        list.append((op, len))?;
    }
    Ok(list.into())
}

fn tag_value_to_py(py: Python<'_>, val: &TagValue) -> PyResult<Py<PyAny>> {
    match val {
        TagValue::Char(c) => c.to_string().into_py_any(py),
        TagValue::Int(n) => (*n).into_py_any(py),
        TagValue::Float(f) => (*f).into_py_any(py),
        TagValue::Str(s) | TagValue::Other(s) => s.clone().into_py_any(py),
    }
}

pub fn records_to_py(py: Python<'_>, records: Vec<Record>) -> PyResult<Vec<Py<PyAlignedSegment>>> {
    records
        .into_iter()
        .map(|r| Py::new(py, PyAlignedSegment::new(r)))
        .collect::<PyResult<Vec<_>>>()
}

pub fn map_core_err<T>(result: samrust_core::Result<T>) -> PyResult<T> {
    result.map_err(to_pyerr)
}
