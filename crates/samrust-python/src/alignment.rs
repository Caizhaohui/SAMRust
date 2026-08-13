//! PyO3 `AlignmentFile` and stats APIs (M3–M6; M11 CRAM read path).

use std::collections::VecDeque;
use std::path::PathBuf;

use numpy::{Element, PyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use pyo3::IntoPyObjectExt;
use samrust_core::{
    coverage_profile_with_filter, depth_blocks, depth_profile, is_cram_path, parallel_count,
    parallel_coverage_profile_with_filter, parallel_depth_profile, parallel_fetch_regions,
    parallel_fetch_wave, parallel_pileup_counts, pileup_counts, whole_file_windows,
    AlignmentReader, CoverageProfile, CramAlignmentReader, DepthProfile, FetchWindow, HeaderDict,
    IndexedAlignmentReader, Interval, PileupCounts, PileupFilter, ReadFilter, Record,
};

use crate::error::to_pyerr;
use crate::segment::{map_core_err, records_to_py, PyAlignedSegment};

const DEFAULT_BATCH_SIZE: usize = 256;
const CRAM_STATS_MSG: &str =
    "CRAM count/depth/coverage/pileup is not in M11 scope (BAM stats only; see DEVELOPMENT_PLAN M11)";

/// Validate a user-supplied 0-based coordinate (pysam raises ValueError on negatives).
pub(crate) fn checked_coord(value: i64, name: &str) -> PyResult<u64> {
    u64::try_from(value)
        .map_err(|_| PyValueError::new_err(format!("{name} must be >= 0, got {value}")))
}

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
    iter_batch: std::vec::IntoIter<Record>,
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
                iter_batch: Vec::new().into_iter(),
            });
        }
        let reader = map_core_err(AlignmentReader::open(&path))?;
        Ok(Self {
            path,
            mode: mode.to_string(),
            reader: Some(reader),
            cram: None,
            closed: false,
            iter_batch: Vec::new().into_iter(),
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
    fn references(&self) -> PyResult<Vec<String>> {
        Ok(self.header_view()?.references().to_vec())
    }

    #[getter]
    fn lengths(&self) -> PyResult<Vec<u64>> {
        Ok(self.header_view()?.lengths().to_vec())
    }

    #[getter]
    fn nreferences(&self) -> PyResult<usize> {
        Ok(self.header_view()?.nreferences())
    }

    /// pysam-style header dict: `{'HD': {...}, 'SQ': [{'SN', 'LN', ...}], 'RG': [...], 'PG': [...]}`.
    #[getter]
    fn header(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let dict = self.header_dict()?;
        header_dict_to_py(py, &dict)
    }

    fn close(&mut self) {
        self.reader = None;
        self.cram = None;
        self.closed = true;
        self.iter_batch = Vec::new().into_iter();
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
        self.iter_batch = Vec::new().into_iter();
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

    /// pysam semantics: iteration continues from the current position;
    /// records already prefetched into the batch buffer are kept.
    fn __iter__(slf: PyRefMut<'_, Self>) -> PyRefMut<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<Py<PyAlignedSegment>>> {
        if self.closed {
            return Err(PyValueError::new_err("I/O operation on closed file"));
        }
        if let Some(record) = self.iter_batch.next() {
            return Ok(Some(Py::new(py, PyAlignedSegment::new(record))?));
        }
        let batch = if let Some(cram) = self.cram.as_mut() {
            py.allow_threads(|| cram.read_batch(DEFAULT_BATCH_SIZE))
                .map_err(to_pyerr)?
        } else if let Some(reader) = self.reader.as_mut() {
            py.allow_threads(|| reader.read_batch(DEFAULT_BATCH_SIZE))
                .map_err(to_pyerr)?
        } else {
            return Err(PyValueError::new_err("I/O operation on closed file"));
        };
        if batch.is_empty() {
            return Ok(None);
        }
        self.iter_batch = batch.into_iter();
        let record = self.iter_batch.next().expect("nonempty batch");
        Ok(Some(Py::new(py, PyAlignedSegment::new(record))?))
    }

    #[pyo3(signature = (contig, start = None, stop = None))]
    fn fetch(
        &mut self,
        py: Python<'_>,
        contig: &str,
        start: Option<i64>,
        stop: Option<i64>,
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
        Ok(PyFetchIterator {
            records: records.into_iter(),
        })
    }

    #[pyo3(signature = (contig, start = None, stop = None, read_callback = "nofilter", threads = 1))]
    fn count(
        &mut self,
        py: Python<'_>,
        contig: &str,
        start: Option<i64>,
        stop: Option<i64>,
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
        start: Option<i64>,
        stop: Option<i64>,
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
        start: Option<i64>,
        stop: Option<i64>,
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
        start: Option<i64>,
        stop: Option<i64>,
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
        start: Option<i64>,
        stop: Option<i64>,
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

    /// Iterate records in batches of `batch_size`.
    ///
    /// `threads=1` streams from a single reader (memory O(batch_size)).
    /// `threads>1` fetches ~1 Mb windows in parallel waves of `threads`
    /// chunks (memory O(threads × window records)) and finishes with the
    /// unmapped tail. Both paths yield every record exactly once, in file
    /// order. `ordered` is accepted for compatibility; output is always
    /// ordered.
    #[pyo3(signature = (batch_size = DEFAULT_BATCH_SIZE, threads = 1, ordered = true))]
    fn iter_batches(
        &self,
        batch_size: usize,
        threads: usize,
        ordered: bool,
    ) -> PyResult<PyBatchIterator> {
        self.ensure_open()?;
        self.reject_cram_stats()?;
        let _ = ordered;
        let batch_size = batch_size.max(1);
        let state = if threads <= 1 {
            let reader = map_core_err(AlignmentReader::open(&self.path))?;
            BatchState::Sequential {
                reader: Box::new(reader),
            }
        } else {
            let windows = whole_file_windows(self.header_view()?);
            BatchState::Parallel {
                path: self.path.clone(),
                windows: windows.into_iter(),
                threads,
                buffer: VecDeque::new(),
                tail_done: false,
            }
        };
        Ok(PyBatchIterator { state, batch_size })
    }

    /// Parallel fetch over multiple regions.
    ///
    /// Overlapping regions are merged; each record is returned exactly once,
    /// in genomic order (contigs in header order). `ordered` is accepted for
    /// compatibility and no longer changes the output order.
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
        let _ = ordered;
        let path = self.path.clone();
        let mut parsed: Vec<(String, Interval)> = Vec::new();
        for item in regions.try_iter()? {
            let item = item?;
            let tuple = item.downcast::<PyTuple>()?;
            if tuple.len() != 3 {
                return Err(PyValueError::new_err(
                    "each region must be (contig, start, stop)",
                ));
            }
            let contig: String = tuple.get_item(0)?.extract()?;
            let start = checked_coord(tuple.get_item(1)?.extract()?, "start")?;
            let stop = checked_coord(tuple.get_item(2)?.extract()?, "stop")?;
            parsed.push((contig, map_core_err(Interval::new(start, stop))?));
        }
        let records =
            py.allow_threads(|| parallel_fetch_regions(&path, &parsed, threads).map_err(to_pyerr))?;
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

    fn header_view(&self) -> PyResult<&samrust_core::Header> {
        if let Some(cram) = self.cram.as_ref() {
            Ok(cram.header())
        } else if let Some(reader) = self.reader.as_ref() {
            Ok(reader.header())
        } else {
            Err(PyValueError::new_err("I/O operation on closed file"))
        }
    }

    fn header_dict(&self) -> PyResult<HeaderDict> {
        if let Some(cram) = self.cram.as_ref() {
            Ok(HeaderDict::from_noodles(cram.raw_header()))
        } else if let Some(reader) = self.reader.as_ref() {
            Ok(HeaderDict::from_noodles(reader.raw_header()))
        } else {
            Err(PyValueError::new_err("I/O operation on closed file"))
        }
    }

    /// Resolve user coordinates to a clamped 0-based half-open interval.
    ///
    /// pysam semantics: negative coordinates raise `ValueError`; `stop`
    /// beyond the contig length is clamped; `start` beyond the contig length
    /// yields an empty interval (zero records).
    fn resolve_region(
        &self,
        contig: &str,
        start: Option<i64>,
        stop: Option<i64>,
    ) -> PyResult<(Interval, u64)> {
        let header = self.header_view()?;
        let ref_id = header.reference_id(contig).map_err(to_pyerr)?;
        let ref_len = header.lengths()[ref_id as usize];
        let start = start
            .map(|v| checked_coord(v, "start"))
            .transpose()?
            .unwrap_or(0)
            .min(ref_len);
        let stop = stop
            .map(|v| checked_coord(v, "stop"))
            .transpose()?
            .unwrap_or(ref_len)
            .min(ref_len);
        let interval = map_core_err(Interval::new(start, stop))?;
        Ok((interval, ref_len))
    }

    fn depth_profile_impl(
        &mut self,
        py: Python<'_>,
        contig: &str,
        start: Option<i64>,
        stop: Option<i64>,
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
    records: std::vec::IntoIter<Record>,
}

#[pymethods]
impl PyFetchIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<Py<PyAlignedSegment>>> {
        match self.records.next() {
            Some(record) => Ok(Some(Py::new(py, PyAlignedSegment::new(record))?)),
            None => Ok(None),
        }
    }
}

/// Batch source for [`PyBatchIterator`].
enum BatchState {
    /// Single-threaded streaming decode (memory O(batch_size)).
    Sequential { reader: Box<AlignmentReader> },
    /// Parallel waves over ~1 Mb windows + the unmapped tail.
    Parallel {
        path: PathBuf,
        windows: std::vec::IntoIter<FetchWindow>,
        threads: usize,
        buffer: VecDeque<Record>,
        tail_done: bool,
    },
}

/// Iterator of record batches returned by `AlignmentFile.iter_batches`.
#[pyclass(name = "BatchIterator", module = "samrust._samrust", unsendable)]
pub struct PyBatchIterator {
    state: BatchState,
    batch_size: usize,
}

impl PyBatchIterator {
    fn next_batch(&mut self, py: Python<'_>) -> PyResult<Option<Vec<Record>>> {
        match &mut self.state {
            BatchState::Sequential { reader } => {
                let batch = py
                    .allow_threads(|| reader.read_batch(self.batch_size))
                    .map_err(to_pyerr)?;
                if batch.is_empty() {
                    Ok(None)
                } else {
                    Ok(Some(batch))
                }
            }
            BatchState::Parallel {
                path,
                windows,
                threads,
                buffer,
                tail_done,
            } => {
                let path: &PathBuf = &*path;
                let wave_size = (*threads).max(1);
                while buffer.len() < self.batch_size {
                    let wave: Vec<FetchWindow> = windows.by_ref().take(wave_size).collect();
                    if wave.is_empty() {
                        if *tail_done {
                            break;
                        }
                        *tail_done = true;
                        let tail = py.allow_threads(|| {
                            let mut indexed =
                                IndexedAlignmentReader::open(path).map_err(to_pyerr)?;
                            indexed.unmapped_tail_records().map_err(to_pyerr)
                        })?;
                        buffer.extend(tail);
                        if buffer.is_empty() {
                            break;
                        }
                        continue;
                    }
                    let results = py
                        .allow_threads(|| parallel_fetch_wave(path, &wave, *threads))
                        .map_err(to_pyerr)?;
                    for records in results {
                        buffer.extend(records);
                    }
                }
                if buffer.is_empty() {
                    return Ok(None);
                }
                let n = self.batch_size.min(buffer.len());
                Ok(Some(buffer.drain(..n).collect()))
            }
        }
    }
}

#[pymethods]
impl PyBatchIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<Py<PyAny>>> {
        match self.next_batch(py)? {
            Some(batch) => {
                let list = PyList::empty(py);
                for seg in records_to_py(py, batch)? {
                    list.append(seg)?;
                }
                Ok(Some(list.into()))
            }
            None => Ok(None),
        }
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

fn header_dict_to_py(py: Python<'_>, dict: &HeaderDict) -> PyResult<Py<PyAny>> {
    let out = PyDict::new(py);
    if !dict.hd.is_empty() {
        let hd = PyDict::new(py);
        for (k, v) in &dict.hd {
            hd.set_item(k, v)?;
        }
        out.set_item("HD", hd)?;
    }
    if !dict.sq.is_empty() {
        let sq = PyList::empty(py);
        for entry in &dict.sq {
            let rec = PyDict::new(py);
            rec.set_item("SN", &entry.sn)?;
            rec.set_item("LN", entry.ln)?;
            for (k, v) in &entry.extra {
                rec.set_item(k, v)?;
            }
            sq.append(rec)?;
        }
        out.set_item("SQ", sq)?;
    }
    if !dict.rg.is_empty() {
        let rg = PyList::empty(py);
        for (id, fields) in &dict.rg {
            let rec = PyDict::new(py);
            rec.set_item("ID", id)?;
            for (k, v) in fields {
                rec.set_item(k, v)?;
            }
            rg.append(rec)?;
        }
        out.set_item("RG", rg)?;
    }
    if !dict.pg.is_empty() {
        let pg = PyList::empty(py);
        for (id, fields) in &dict.pg {
            let rec = PyDict::new(py);
            rec.set_item("ID", id)?;
            for (k, v) in fields {
                rec.set_item(k, v)?;
            }
            pg.append(rec)?;
        }
        out.set_item("PG", pg)?;
    }
    Ok(out.into())
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
