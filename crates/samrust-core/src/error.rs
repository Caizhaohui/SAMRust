//! Error types for SAMRust core.

use thiserror::Error;

/// Result alias used throughout `samrust-core`.
pub type Result<T> = std::result::Result<T, SamRustError>;

/// Top-level error type for SAMRust.
#[derive(Debug, Error)]
pub enum SamRustError {
    /// Invalid user input or configuration.
    #[error("invalid argument: {0}")]
    InvalidArgument(String),

    /// I/O failure (files, indexes, etc.).
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// BAM index missing or unreadable.
    #[error("missing index for {0}")]
    MissingIndex(String),

    /// Feature not yet implemented (milestones after M0).
    #[error("not implemented: {0}")]
    NotImplemented(&'static str),
}

impl SamRustError {
    /// Map a low-level I/O error from index loading into [`MissingIndex`] when appropriate.
    pub fn from_index_io(path: &str, err: std::io::Error) -> Self {
        if err.kind() == std::io::ErrorKind::NotFound {
            Self::MissingIndex(path.to_string())
        } else {
            Self::Io(err)
        }
    }
}
