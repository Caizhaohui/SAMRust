//! PyO3 `AlignmentFile` and stats APIs (M3–M6; M11 CRAM read path).

use std::path::PathBuf;

use numpy::{Element, PyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyList, PyTuple};
use pyo3::IntoPyObjectExt;
use samrust_core::{
    coverage_profile_with_filter, depth_blocks, depth_profile, is_cram_path, parallel_count,
    parallel_coverage_profile_with_filter, parallel_depth_profile, parallel_fetch_records,
    parallel_pileup_counts, pileup_counts, AlignmentReader, CoverageProfile, CramAlignmentReader,
    DepthProfile, IndexedAlignmentReader, Interval, PileupCounts, PileupFilter, ReadFilter, Record,
    Scheduler,
};

use crate::error::to_pyerr;
use crate::segment::{map_core_err, records_to_py, PyAlignedSegment};

const DEFAULT_BATCH_SIZE: usize = 256;
const CRAM_STATS_MSG: &str =
    "CRAM count/depth/coverage/pileup is not in M11 scope (BAM stats only; see DEVELOPMENT_PLAN M11)";

/// pysam-compatible BAM / CRAM reader.
///
/// `unsendable`: noodles CSI/CRAM index trait objects are not `Send`/`Sync`.
#[pyclass(name = "AlignmentFile", module = "samrust._samrust", unsendable)]
pub struct PyAlignmentFile {
    path: PathBuf,
    mode: String,
    reader: Option<AlignmentReader>,
    cram: Option<CramAlignmentReader>,
    closed: bool,
    iter_batch: Vec<Record>,
    iter_pos: usize,
}

#[pymethods]
impl PyAlignmentFile {
    #[new]
    #[pyo3(signature = (filename, mode = "rb", reference_filename = None))]
    fn new(filename: &str, mode: &str, reference_filename: Option<&str>) -> PyResult<Self> {
        let path = PathBuf::from(filename);
        let is_cram = is_cram_path(&path);
        if mode != "rb" && !(is_cram && mode == "rc") {
            return Err(PyValueError::new_err(format!(
                "only mode 'rb' is supported (CRAM also accepts 'rc'), got {mode:?}"
            )));
        }
        if is_cram {
            let fasta = reference_filename
                .map(PathBuf::from)
                .or_else(|| sibling_fasta(&path))
                .ok_or_else(|| {
                    PyValueError::new_err(
                        "CRAM requires reference_filename= or a sibling .fa/.fasta/.fna",
                    )
                })?;
            let cram = map_core_err(CramAlignmentReader::open(&path, &fasta))?;
            return Ok(Self {
                path,
                mode: mode.to_string(),
                reader: None,
                cram: Some(cram),
                closed: false,
                iter_batch: Vec::new(),
                iter_pos: 0,
            });
        }
        let reader = map_core_err(AlignmentReader::open(&path))?;
        Ok(Self {
            path,
            mode: mode.to_string(),
            reader: Some(reader),
            cram: None,
            closed: false,
            iter_batch: Vec::new(),
            iter_pos: 0,
        })
    }

    #[getter]
    fn filename(&self) -> String {
        self.path.display().to_string()
    }

    #[getter]
    fn mode(&self) -> &str {
        &self.mode
    }

    #[getter]
    fn references(&self) -> Vec<String> {
        self.header_view().references().to_vec()
    }

    #[getter]
    fn lengths(&self) -> Vec<u64> {
        self.header_view().lengths().to_vec()
    }

    #[getter]
    fn nreferences(&self) -> usize {
        self.header_view().nreferences()
    }

    #[getter]
    fn header(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let header = self.header_view();
        let dict = pyo3::types::PyDict::new(py);
        let refs = PyList::new(py, header.references())?;
        let lens = PyList::new(py, header.lengths())?;
        dict.set_item("nreferences", header.nreferences())?;
        dict.set_item("references", refs)?;
        dict.set_item("lengths", lens)?;
        Ok(dict.into())
    }

    fn close(&mut self) {
        self.reader = None;
        self.cram = None;
        self.closed = true;
        self.iter_batch.clear();
        self.iter_pos = 0;
    }

    fn reset(&mut self) -> PyResult<()> {
        if self.closed {
            return Err(PyValueError::new_err("I/O operation on closed file"));
        }
        if self.cram.is_some() {
            let fasta = self
                .cram
                .as_ref()
                .map(|c| c.fasta_path().to_path_buf())
                .expect("cram");
            let cram = map_core_err(CramAlignmentReader::open(&self.path, &fasta))?;
            self.cram = Some(cram);
        } else {
            let reader = map_core_err(AlignmentReader::open(&self.path))?;
            self.reader = Some(reader);
        }
        self.iter_batch.clear();
        self.iter_pos = 0;
        Ok(())
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

    fn __iter__(mut slf: PyRefMut<'_, Self>) -> PyRefMut<'_, Self> {
        slf.iter_batch.clear();
        slf.iter_pos = 0;
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<Py<PyAlignedSegment>>> {
        if self.closed {
            return Err(PyValueError::new_err("I/O operation on closed file"));
        }
        if self.iter_pos >= self.iter_batch.len() {
            let path = self.path.clone();
            // Prefer existing sequential reader when available; otherwise reopen.
            let batch = if let Some(cram) = self.cram.as_mut() {
                map_core_err(cram.read_batch(DEFAULT_BATCH_SIZE))?
            } else if let Some(reader) = self.reader.as_mut() {
                map_core_err(reader.read_batch(DEFAULT_BATCH_SIZE))?
            } else {
                py.allow_threads(|| {
                    let mut reader = AlignmentReader::open(&path).map_err(to_pyerr)?;
                    reader.read_batch(DEFAULT_BATCH_SIZE).map_err(to_pyerr)
                })?
            };
            if batch.is_empty() {
                return Ok(None);
            }
            self.iter_batch = batch;
            self.iter_pos = 0;
        }
        let record = self.iter_batch[self.iter_pos].clone();
        self.iter_pos += 1;
        Ok(Some(Py::new(py, PyAlignedSegment::new(record))?))
    }

    #[pyo3(signature = (contig, start = None, stop = None))]
    fn fetch(
        &mut self,
        py: Python<'_>,
        contig: &str,
        start: Option<u64>,
        stop: Option<u64>,
    ) -> PyResult<PyFetchIterator> {
        self.ensure_open()?;
        let (interval, _) = self.resolve_region(contig, start, stop)?;
        let records = if let Some(cram) = self.cram.as_mut() {
            map_core_err(cram.fetch_records(contig, interval))?
        } else {
            let path = self.path.clone();
            let contig = contig.to_string();
            py.allow_threads(|| {
                let mut indexed = IndexedAlignmentReader::open(&path).map_err(to_pyerr)?;
                indexed.fetch_records(&contig, interval).map_err(to_pyerr)
            })?
        };
        Ok(PyFetchIterator { records, pos: 0 })
    }

    #[pyo3(signature = (contig, start = None, stop = None, read_callback = "nofilter", threads = 1))]
    fn count(
        &mut self,
        py: Python<'_>,
        contig: &str,
        start: Option<u64>,
        stop: Option<u64>,
        read_callback: &str,
        threads: usize,
    ) -> PyResult<u64> {
        self.ensure_open()?;
        self.reject_cram_stats()?;
        let (interval, _) = self.resolve_region(contig, start, stop)?;
        let filter = parse_read_filter(read_callback)?;
        let path = self.path.clone();
        let contig = contig.to_string();
        py.allow_threads(|| {
            parallel_count(&path, &contig, interval, filter, threads).map_err(to_pyerr)
        })
    }

    #[pyo3(signature = (contig, start = None, stop = None, quality_threshold = 15, read_callback = "all", threads = 1))]
    #[allow(clippy::too_many_arguments)]
    fn count_coverage(
        &mut self,
        py: Python<'_>,
        contig: &str,
        start: Option<u64>,
        stop: Option<u64>,
        quality_threshold: u8,
        read_callback: &str,
        threads: usize,
    ) -> PyResult<Py<PyAny>> {
        self.ensure_open()?;
        self.reject_cram_stats()?;
        let filter = parse_read_filter(read_callback)?;
        let (interval, _) = self.resolve_region(contig, start, stop)?;
        let path = self.path.clone();
        let contig = contig.to_string();
        let profile = if threads <= 1 {
            py.allow_threads(|| {
                let mut indexed = IndexedAlignmentReader::open(&path).map_err(to_pyerr)?;
                coverage_profile_with_filter(
                    &mut indexed,
                    &contig,
                    interval,
                    quality_threshold,
                    filter,
                )
                .map_err(to_pyerr)
            })?
        } else {
            py.allow_threads(|| {
                parallel_coverage_profile_with_filter(
                    &path,
                    &contig,
                    interval,
                    quality_threshold,
                    filter,
                    threads,
                )
                .map_err(to_pyerr)
            })?
        };
        coverage_to_py(py, profile)
    }

    #[pyo3(signature = (contig, start = None, stop = None, threads = 1))]
    fn depth_blocks(
        &mut self,
        py: Python<'_>,
        contig: &str,
        start: Option<u64>,
        stop: Option<u64>,
        threads: usize,
    ) -> PyResult<Py<PyAny>> {
        let profile = self.depth_profile_impl(py, contig, start, stop, threads)?;
        let blocks = depth_blocks(&profile);
        let list = PyList::empty(py);
        for (start, length, depth) in blocks {
            list.append((start, length, depth))?;
        }
        Ok(list.into())
    }

    #[pyo3(signature = (contig, start = None, stop = None, threads = 1))]
    fn depth_numpy(
        &mut self,
        py: Python<'_>,
        contig: &str,
        start: Option<u64>,
        stop: Option<u64>,
        threads: usize,
    ) -> PyResult<Py<PyAny>> {
        let profile = self.depth_profile_impl(py, contig, start, stop, threads)?;
        depth_to_numpy(py, profile)
    }

    #[pyo3(signature = (contig, start = None, stop = None, min_base_quality = 0, min_mapping_quality = 0, threads = 1))]
    #[allow(clippy::too_many_arguments)]
    fn pileup_counts(
        &mut self,
        py: Python<'_>,
        contig: &str,
        start: Option<u64>,
        stop: Option<u64>,
        min_base_quality: u8,
        min_mapping_quality: u8,
        threads: usize,
    ) -> PyResult<Py<PyAny>> {
        self.ensure_open()?;
        self.reject_cram_stats()?;
        let (interval, _) = self.resolve_region(contig, start, stop)?;
        let filter = PileupFilter {
            min_base_quality,
            min_mapping_quality,
            ..PileupFilter::default()
        };
        let path = self.path.clone();
        let contig = contig.to_string();
        let counts = if threads <= 1 {
            py.allow_threads(|| {
                let mut indexed = IndexedAlignmentReader::open(&path).map_err(to_pyerr)?;
                pileup_counts(&mut indexed, &contig, interval, filter).map_err(to_pyerr)
            })?
        } else {
            py.allow_threads(|| {
                parallel_pileup_counts(&path, &contig, interval, filter, threads).map_err(to_pyerr)
            })?
        };
        pileup_to_py(py, counts)
    }

    #[pyo3(signature = (batch_size = DEFAULT_BATCH_SIZE, threads = 1, ordered = true))]
    fn iter_batches(
        &self,
        py: Python<'_>,
        batch_size: usize,
        threads: usize,
        ordered: bool,
    ) -> PyResult<Py<PyAny>> {
        self.ensure_open()?;
        self.reject_cram_stats()?;
        let path = self.path.clone();
        let batch_size = batch_size.max(1);
        if threads <= 1 {
            let batch = py.allow_threads(|| {
                let mut reader = AlignmentReader::open(&path).map_err(to_pyerr)?;
                reader.read_batch(batch_size).map_err(to_pyerr)
            })?;
            let py_records = records_to_py(py, batch)?;
            let list = PyList::empty(py);
            list.append(py_records)?;
            return Ok(list.into());
        }
        let header = self.header_view();
        let mut region_chunks = Vec::new();
        let sched = Scheduler::default();
        for (name, len) in header.references().iter().zip(header.lengths()) {
            region_chunks.extend(map_core_err(sched.chunk_interval(
                name.clone(),
                0,
                *len,
                threads,
            ))?);
        }
        let batches = py.allow_threads(|| {
            samrust_core::parallel_map_regions(
                &path,
                region_chunks,
                threads,
                ordered,
                |_h, indexed, chunk| indexed.fetch_records(&chunk.contig, chunk.interval),
            )
            .map_err(to_pyerr)
        })?;
        let list = PyList::empty(py);
        for batch in batches {
            let py_records = records_to_py(py, batch)?;
            list.append(py_records)?;
        }
        Ok(list.into())
    }

    #[pyo3(signature = (regions, threads = 1, ordered = true))]
    fn parallel_fetch(
        &self,
        py: Python<'_>,
        regions: &Bound<'_, PyAny>,
        threads: usize,
        ordered: bool,
    ) -> PyResult<Py<PyAny>> {
        self.ensure_open()?;
        self.reject_cram_stats()?;
        let path = self.path.clone();
        let header = self.header_view();
        let mut chunks = Vec::new();
        let sched = Scheduler::default();
        for item in regions.try_iter()? {
            let item = item?;
            let tuple = item.downcast::<PyTuple>()?;
            if tuple.len() != 3 {
                return Err(PyValueError::new_err(
                    "each region must be (contig, start, stop)",
                ));
            }
            let contig: String = tuple.get_item(0)?.extract()?;
            let start: u64 = tuple.get_item(1)?.extract()?;
            let stop: u64 = tuple.get_item(2)?.extract()?;
            header.reference_id(&contig).map_err(to_pyerr)?;
            chunks.extend(map_core_err(
                sched.chunk_interval(contig, start, stop, threads),
            )?);
        }
        let records = py.allow_threads(|| {
            parallel_fetch_records(&path, chunks, threads, ordered).map_err(to_pyerr)
        })?;
        let py_records = records_to_py(py, records)?;
        py_records.into_py_any(py)
    }
}

impl PyAlignmentFile {
    fn ensure_open(&self) -> PyResult<()> {
        if self.closed || (self.reader.is_none() && self.cram.is_none()) {
            Err(PyValueError::new_err("I/O operation on closed file"))
        } else {
            Ok(())
        }
    }

    fn reject_cram_stats(&self) -> PyResult<()> {
        if self.cram.is_some() {
            Err(to_pyerr(samrust_core::SamRustError::NotImplemented(
                CRAM_STATS_MSG,
            )))
        } else {
            Ok(())
        }
    }

    fn header_view(&self) -> &samrust_core::Header {
        if let Some(cram) = self.cram.as_ref() {
            cram.header()
        } else {
            self.reader
                .as_ref()
                .map(AlignmentReader::header)
                .expect("closed file without header")
        }
    }

    fn resolve_region(
        &self,
        contig: &str,
        start: Option<u64>,
        stop: Option<u64>,
    ) -> PyResult<(Interval, u64)> {
        let header = self.header_view();
        let ref_id = header.reference_id(contig).map_err(to_pyerr)?;
        let ref_len = header.lengths()[ref_id as usize];
        let start = start.unwrap_or(0);
        let stop = stop.unwrap_or(ref_len);
        let interval = map_core_err(Interval::new(start, stop))?;
        Ok((interval, ref_len))
    }

    fn depth_profile_impl(
        &mut self,
        py: Python<'_>,
        contig: &str,
        start: Option<u64>,
        stop: Option<u64>,
        threads: usize,
    ) -> PyResult<DepthProfile> {
        self.ensure_open()?;
        self.reject_cram_stats()?;
        let (interval, _) = self.resolve_region(contig, start, stop)?;
        let path = self.path.clone();
        let contig = contig.to_string();
        if threads <= 1 {
            py.allow_threads(|| {
                let mut indexed = IndexedAlignmentReader::open(&path).map_err(to_pyerr)?;
                depth_profile(&mut indexed, &contig, interval).map_err(to_pyerr)
            })
        } else {
            py.allow_threads(|| {
                parallel_depth_profile(&path, &contig, interval, threads).map_err(to_pyerr)
            })
        }
    }
}

/// Iterator returned by `fetch`.
#[pyclass(name = "FetchIterator", module = "samrust._samrust", unsendable)]
pub struct PyFetchIterator {
    records: Vec<Record>,
    pos: usize,
}

#[pymethods]
impl PyFetchIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<Py<PyAlignedSegment>>> {
        if self.pos >= self.records.len() {
            return Ok(None);
        }
        let record = self.records[self.pos].clone();
        self.pos += 1;
        Ok(Some(Py::new(py, PyAlignedSegment::new(record))?))
    }
}

fn parse_read_filter(name: &str) -> PyResult<ReadFilter> {
    match name {
        "nofilter" => Ok(ReadFilter::NoFilter),
        "all" => Ok(ReadFilter::All),
        other => Err(PyValueError::new_err(format!(
            "unsupported read_callback: {other:?}"
        ))),
    }
}

fn sibling_fasta(cram: &std::path::Path) -> Option<std::path::PathBuf> {
    let stem = cram.file_stem()?.to_string_lossy();
    let dir = cram.parent().unwrap_or_else(|| std::path::Path::new("."));
    for ext in ["fa", "fasta", "fna"] {
        let candidate = dir.join(format!("{stem}.{ext}"));
        if candidate.is_file() {
            return Some(candidate);
        }
    }
    None
}

fn coverage_to_py(py: Python<'_>, profile: CoverageProfile) -> PyResult<Py<PyAny>> {
    if numpy_available(py) {
        let a = vec_to_numpy(py, profile.a)?;
        let c = vec_to_numpy(py, profile.c)?;
        let g = vec_to_numpy(py, profile.g)?;
        let t = vec_to_numpy(py, profile.t)?;
        Ok(PyTuple::new(py, [a, c, g, t])?.into())
    } else {
        Ok(PyTuple::new(py, [profile.a, profile.c, profile.g, profile.t])?.into())
    }
}

fn depth_to_numpy(py: Python<'_>, profile: DepthProfile) -> PyResult<Py<PyAny>> {
    if numpy_available(py) {
        Ok(vec_to_numpy(py, profile.depth)?.into())
    } else {
        profile.depth.into_py_any(py)
    }
}

fn pileup_to_py(py: Python<'_>, counts: PileupCounts) -> PyResult<Py<PyAny>> {
    let dict = pyo3::types::PyDict::new(py);
    if numpy_available(py) {
        dict.set_item("A", vec_to_numpy(py, counts.a)?)?;
        dict.set_item("C", vec_to_numpy(py, counts.c)?)?;
        dict.set_item("G", vec_to_numpy(py, counts.g)?)?;
        dict.set_item("T", vec_to_numpy(py, counts.t)?)?;
        dict.set_item("N", vec_to_numpy(py, counts.n)?)?;
        dict.set_item("depth", vec_to_numpy(py, counts.depth)?)?;
    } else {
        dict.set_item("A", counts.a)?;
        dict.set_item("C", counts.c)?;
        dict.set_item("G", counts.g)?;
        dict.set_item("T", counts.t)?;
        dict.set_item("N", counts.n)?;
        dict.set_item("depth", counts.depth)?;
    }
    Ok(dict.into())
}

fn numpy_available(py: Python<'_>) -> bool {
    py.import("numpy").is_ok()
}

fn vec_to_numpy<T>(py: Python<'_>, data: Vec<T>) -> PyResult<Py<PyArray1<T>>>
where
    T: Element,
{
    Ok(PyArray1::from_vec(py, data).unbind())
}
