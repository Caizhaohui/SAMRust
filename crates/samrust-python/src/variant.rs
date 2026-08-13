//! PyO3 `VariantFile` / `VariantRecord` (M9).

use std::path::PathBuf;
use std::sync::Arc;

use pyo3::exceptions::{PyKeyError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use pyo3::IntoPyObjectExt;
use samrust_core::{Interval, VariantInfoValue, VariantReader, VariantRecord, VariantSample};

use crate::alignment::checked_coord;
use crate::error::to_pyerr;

/// pysam-compatible VCF/BCF reader.
///
/// `unsendable`: noodles index trait objects are not `Send`.
#[pyclass(name = "VariantFile", module = "samrust._samrust", unsendable)]
pub struct PyVariantFile {
    path: PathBuf,
    reader: Option<VariantReader>,
    closed: bool,
    iter_records: std::vec::IntoIter<VariantRecord>,
}

#[pymethods]
impl PyVariantFile {
    #[new]
    #[pyo3(signature = (filename, mode = "r"))]
    fn new(filename: &str, mode: &str) -> PyResult<Self> {
        if mode != "r" && mode != "rb" {
            return Err(PyValueError::new_err(format!(
                "only modes 'r'/'rb' are supported, got {mode:?}"
            )));
        }
        let path = PathBuf::from(filename);
        let reader = VariantReader::open(&path).map_err(to_pyerr)?;
        Ok(Self {
            path,
            reader: Some(reader),
            closed: false,
            iter_records: Vec::new().into_iter(),
        })
    }

    #[getter]
    fn filename(&self) -> String {
        self.path.display().to_string()
    }

    #[getter]
    fn header(&self, py: Python<'_>) -> PyResult<Py<PyVariantHeader>> {
        let reader = self.ensure_open()?;
        let h = reader.header();
        Py::new(
            py,
            PyVariantHeader {
                samples: h.samples.clone(),
                contigs: h.contigs.clone(),
            },
        )
    }

    #[getter]
    fn samples(&self) -> PyResult<Vec<String>> {
        Ok(self.ensure_open()?.header().samples.clone())
    }

    fn close(&mut self) {
        self.reader = None;
        self.closed = true;
        self.iter_records = Vec::new().into_iter();
    }

    fn __enter__(slf: PyRefMut<'_, Self>) -> PyRefMut<'_, Self> {
        slf
    }

    fn __exit__(
        &mut self,
        _exc_type: Bound<'_, PyAny>,
        _exc_val: Bound<'_, PyAny>,
        _exc_tb: Bound<'_, PyAny>,
    ) -> PyResult<bool> {
        self.close();
        Ok(false)
    }

    fn __iter__<'a>(mut slf: PyRefMut<'a, Self>, py: Python<'_>) -> PyResult<PyRefMut<'a, Self>> {
        let reader = slf.ensure_open()?;
        let recs = py.allow_threads(|| reader.records().map_err(to_pyerr))?;
        slf.iter_records = recs.into_iter();
        Ok(slf)
    }

    fn __next__(
        mut slf: PyRefMut<'_, Self>,
        py: Python<'_>,
    ) -> PyResult<Option<Py<PyVariantRecord>>> {
        let Some(rec) = slf.iter_records.next() else {
            return Ok(None);
        };
        let names = slf
            .reader
            .as_ref()
            .map(|r| Arc::new(r.header().samples.clone()))
            .unwrap_or_default();
        Ok(Some(Py::new(
            py,
            PyVariantRecord {
                inner: rec,
                sample_names: names,
            },
        )?))
    }

    /// Region fetch, 0-based half-open.
    ///
    /// pysam semantics: negative coordinates raise `ValueError`; `stop` is
    /// clamped to the contig length when the header provides one. When the
    /// header has no contig length and `stop` is omitted, an unbounded fetch
    /// falls back to a sequential scan (tabix/CSI cannot express it).
    #[pyo3(signature = (contig = None, start = None, stop = None))]
    fn fetch(
        &mut self,
        py: Python<'_>,
        contig: Option<&str>,
        start: Option<i64>,
        stop: Option<i64>,
    ) -> PyResult<PyVariantFetchIterator> {
        let reader = self.ensure_open()?;
        let records = if let Some(contig) = contig {
            // 0 = header carries no length for this contig.
            let len = reader.header().contig_length(contig).map_err(to_pyerr)?;
            let start = start
                .map(|v| checked_coord(v, "start"))
                .transpose()?
                .unwrap_or(0);
            let stop = stop.map(|v| checked_coord(v, "stop")).transpose()?;
            match (stop, len) {
                (None, 0) => {
                    py.allow_threads(|| reader.fetch_from(contig, start).map_err(to_pyerr))?
                }
                (s, l) => {
                    let (start, stop) = if l > 0 {
                        (start.min(l), s.unwrap_or(l).min(l))
                    } else {
                        (start, s.expect("stop is Some when len == 0"))
                    };
                    let interval = Interval::new(start, stop).map_err(to_pyerr)?;
                    py.allow_threads(|| reader.fetch(contig, interval).map_err(to_pyerr))?
                }
            }
        } else {
            py.allow_threads(|| reader.records().map_err(to_pyerr))?
        };
        let names = Arc::new(reader.header().samples.clone());
        Ok(PyVariantFetchIterator {
            records: records.into_iter(),
            names,
        })
    }
}

impl PyVariantFile {
    fn ensure_open(&self) -> PyResult<&VariantReader> {
        self.reader
            .as_ref()
            .ok_or_else(|| PyValueError::new_err("I/O operation on closed file"))
    }
}

#[pyclass(name = "VariantHeader", module = "samrust._samrust")]
#[derive(Clone)]
pub struct PyVariantHeader {
    samples: Vec<String>,
    contigs: Vec<String>,
}

#[pymethods]
impl PyVariantHeader {
    #[getter]
    fn samples(&self) -> Vec<String> {
        self.samples.clone()
    }

    #[getter]
    fn contigs(&self) -> Vec<String> {
        self.contigs.clone()
    }
}

#[pyclass(name = "VariantRecord", module = "samrust._samrust")]
#[derive(Clone)]
pub struct PyVariantRecord {
    inner: VariantRecord,
    sample_names: Arc<Vec<String>>,
}

#[pymethods]
impl PyVariantRecord {
    #[getter]
    fn chrom(&self) -> &str {
        &self.inner.chrom
    }

    #[getter]
    fn contig(&self) -> &str {
        &self.inner.chrom
    }

    /// 1-based POS (pysam `VariantRecord.pos`).
    #[getter]
    fn pos(&self) -> u64 {
        self.inner.pos()
    }

    /// 0-based start.
    #[getter]
    fn start(&self) -> u64 {
        self.inner.start
    }

    /// 0-based exclusive stop.
    #[getter]
    fn stop(&self) -> u64 {
        self.inner.stop
    }

    #[getter]
    fn id(&self) -> Option<&str> {
        self.inner.id.as_deref()
    }

    #[getter]
    fn r#ref(&self) -> &str {
        &self.inner.ref_allele
    }

    #[getter]
    fn alts(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        PyTuple::new(py, &self.inner.alts)?.into_py_any(py)
    }

    #[getter]
    fn alleles(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let mut all = Vec::with_capacity(1 + self.inner.alts.len());
        all.push(self.inner.ref_allele.clone());
        all.extend(self.inner.alts.iter().cloned());
        PyTuple::new(py, &all)?.into_py_any(py)
    }

    #[getter]
    fn qual(&self) -> Option<f32> {
        self.inner.qual
    }

    #[getter]
    fn filter(&self) -> Vec<String> {
        self.inner.filter.clone()
    }

    #[getter]
    fn format(&self) -> Vec<String> {
        self.inner.format.clone()
    }

    #[getter]
    fn info(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let dict = pyo3::types::PyDict::new(py);
        for (k, v) in &self.inner.info {
            dict.set_item(k, info_to_py(py, v)?)?;
        }
        Ok(dict.into())
    }

    #[getter]
    fn samples(&self) -> PyVariantRecordSamples {
        PyVariantRecordSamples {
            names: self.sample_names.clone(),
            samples: self.inner.samples.clone(),
        }
    }
}

#[pyclass(name = "_VariantRecordSamples", module = "samrust._samrust")]
#[derive(Clone)]
pub struct PyVariantRecordSamples {
    names: Arc<Vec<String>>,
    samples: Vec<VariantSample>,
}

#[pymethods]
impl PyVariantRecordSamples {
    fn __len__(&self) -> usize {
        self.samples.len()
    }

    fn __getitem__(&self, py: Python<'_>, key: Bound<'_, PyAny>) -> PyResult<Py<PyVariantSample>> {
        let idx = if let Ok(i) = key.extract::<isize>() {
            if i < 0 {
                return Err(PyKeyError::new_err(i));
            }
            i as usize
        } else if let Ok(name) = key.extract::<String>() {
            self.names
                .iter()
                .position(|n| n == &name)
                .ok_or_else(|| PyKeyError::new_err(name))?
        } else {
            return Err(PyKeyError::new_err("sample key must be int or str"));
        };
        let sample = self
            .samples
            .get(idx)
            .cloned()
            .ok_or_else(|| PyKeyError::new_err(idx))?;
        Py::new(py, PyVariantSample { inner: sample })
    }
}

#[pyclass(name = "_VariantSample", module = "samrust._samrust")]
#[derive(Clone)]
pub struct PyVariantSample {
    inner: VariantSample,
}

#[pymethods]
impl PyVariantSample {
    fn __getitem__(&self, py: Python<'_>, key: &str) -> PyResult<Py<PyAny>> {
        match key {
            "GT" => match &self.inner.gt {
                Some(gt) => alleles_to_tuple(py, gt),
                None => Ok(py.None()),
            },
            "DP" => match self.inner.dp {
                Some(n) => n.into_py_any(py),
                None => Ok(py.None()),
            },
            "AD" => match &self.inner.ad {
                Some(ad) => optional_ints_to_tuple(py, ad),
                None => Ok(py.None()),
            },
            _ => Err(PyKeyError::new_err(key.to_string())),
        }
    }
}

#[pyclass(name = "_VariantFetchIterator", module = "samrust._samrust")]
pub struct PyVariantFetchIterator {
    records: std::vec::IntoIter<VariantRecord>,
    names: Arc<Vec<String>>,
}

#[pymethods]
impl PyVariantFetchIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<Py<PyVariantRecord>>> {
        match self.records.next() {
            Some(rec) => Ok(Some(Py::new(
                py,
                PyVariantRecord {
                    inner: rec,
                    sample_names: self.names.clone(),
                },
            )?)),
            None => Ok(None),
        }
    }
}

fn info_to_py(py: Python<'_>, value: &VariantInfoValue) -> PyResult<Py<PyAny>> {
    match value {
        VariantInfoValue::Integer(n) => n.into_py_any(py),
        VariantInfoValue::Float(n) => (*n as f64).into_py_any(py),
        VariantInfoValue::Flag => true.into_py_any(py),
        VariantInfoValue::Character(c) => c.to_string().into_py_any(py),
        VariantInfoValue::String(s) => s.into_py_any(py),
        VariantInfoValue::IntegerArray(v) => optional_ints_to_tuple(py, v),
        VariantInfoValue::FloatArray(v) => {
            let mut items = Vec::with_capacity(v.len());
            for x in v {
                items.push(match x {
                    Some(n) => (*n as f64).into_py_any(py)?,
                    None => py.None(),
                });
            }
            PyTuple::new(py, items)?.into_py_any(py)
        }
        VariantInfoValue::StringArray(v) => {
            let mut items = Vec::with_capacity(v.len());
            for x in v {
                items.push(match x {
                    Some(s) => s.into_py_any(py)?,
                    None => py.None(),
                });
            }
            PyTuple::new(py, items)?.into_py_any(py)
        }
    }
}

fn alleles_to_tuple(py: Python<'_>, gt: &[Option<i32>]) -> PyResult<Py<PyAny>> {
    let mut items = Vec::with_capacity(gt.len());
    for a in gt {
        items.push(match a {
            Some(n) => n.into_py_any(py)?,
            None => py.None(),
        });
    }
    PyTuple::new(py, items)?.into_py_any(py)
}

fn optional_ints_to_tuple(py: Python<'_>, values: &[Option<i32>]) -> PyResult<Py<PyAny>> {
    let mut items = Vec::with_capacity(values.len());
    for a in values {
        items.push(match a {
            Some(n) => n.into_py_any(py)?,
            None => py.None(),
        });
    }
    PyTuple::new(py, items)?.into_py_any(py)
}
