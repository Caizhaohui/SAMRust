//! Auxiliary tag (BAM data field) helpers.

use std::collections::BTreeMap;
use std::fmt;

use noodles::sam::alignment::record::data::field::Value;

use crate::error::{Result, SamRustError};

/// Owned auxiliary tags keyed by two-letter tag name.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Tags {
    values: BTreeMap<String, TagValue>,
}

/// Simplified tag value for parity dumps / later Python mapping.
#[derive(Debug, Clone, PartialEq)]
pub enum TagValue {
    /// Character (`A`).
    Char(char),
    /// Integer (any BAM integer width).
    Int(i64),
    /// Float (`f`).
    Float(f32),
    /// String / hex (`Z` / `H`).
    Str(String),
    /// Other / array — stringified for M2 dumps.
    Other(String),
}

impl Tags {
    /// Empty tag map.
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of tags.
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Whether no tags are present.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Get a tag by name (`"NM"`, `"RG"`, ...).
    pub fn get(&self, tag: &str) -> Option<&TagValue> {
        self.values.get(tag)
    }

    /// Iterate tags in sorted name order.
    pub fn iter(&self) -> impl Iterator<Item = (&str, &TagValue)> + '_ {
        self.values.iter().map(|(k, v)| (k.as_str(), v))
    }

    /// Insert a tag.
    pub fn insert(&mut self, tag: impl Into<String>, value: TagValue) {
        self.values.insert(tag.into(), value);
    }

    /// Build from noodles BAM data fields.
    pub fn from_noodles_data<'a, I>(iter: I) -> Result<Self>
    where
        I: IntoIterator<
            Item = std::io::Result<(noodles::sam::alignment::record::data::field::Tag, Value<'a>)>,
        >,
    {
        let mut tags = Self::new();
        for item in iter {
            let (tag, value) = item.map_err(SamRustError::from)?;
            let name = std::str::from_utf8(tag.as_ref())
                .unwrap_or("??")
                .to_string();
            tags.insert(name, TagValue::from_noodles(value));
        }
        Ok(tags)
    }
}

impl TagValue {
    fn from_noodles(value: Value<'_>) -> Self {
        if let Some(n) = value.as_int() {
            return Self::Int(n);
        }
        match value {
            Value::Character(b) => Self::Char(char::from(b)),
            Value::Float(f) => Self::Float(f),
            Value::String(s) | Value::Hex(s) => Self::Str(s.to_string()),
            other => Self::Other(format!("{other:?}")),
        }
    }

    /// JSON-ish literal for dump/parity scripts.
    pub fn to_parity_string(&self) -> String {
        match self {
            Self::Char(c) => format!("\"{c}\""),
            Self::Int(n) => n.to_string(),
            Self::Float(f) => {
                // Match Python-ish compact floats when possible.
                if f.fract() == 0.0 {
                    format!("{f:.1}")
                } else {
                    f.to_string()
                }
            }
            Self::Str(s) => format!("\"{}\"", escape_json(s)),
            Self::Other(s) => format!("\"{}\"", escape_json(s)),
        }
    }
}

impl fmt::Display for TagValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Char(c) => write!(f, "{c}"),
            Self::Int(n) => write!(f, "{n}"),
            Self::Float(x) => write!(f, "{x}"),
            Self::Str(s) | Self::Other(s) => write!(f, "{s}"),
        }
    }
}

fn escape_json(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}
