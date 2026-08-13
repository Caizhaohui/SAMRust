//! Integration tests for linear BAM reading (M2).

use std::path::PathBuf;

use samrust_core::AlignmentReader;

fn fixture_bam() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/small.bam")
}

#[test]
fn opens_fixture_and_reads_header() {
    let reader = AlignmentReader::open(fixture_bam()).expect("open fixture");
    let header = reader.header();
    assert_eq!(header.nreferences(), 2);
    assert_eq!(
        header.references(),
        &["chr1".to_string(), "chr2".to_string()]
    );
    assert_eq!(header.lengths(), &[1000, 500]);
}

#[test]
fn iterates_fixture_records_with_expected_flags() {
    let mut reader = AlignmentReader::open(fixture_bam()).expect("open fixture");
    let records: Vec<_> = reader.records().map(|r| r.expect("record")).collect();
    assert_eq!(records.len(), 14);

    let by_name = |name: &str| {
        records
            .iter()
            .find(|r| r.query_name() == name)
            .unwrap_or_else(|| panic!("missing {name}"))
    };

    let pair = by_name("pair1");
    assert_eq!(pair.flag(), 99);
    assert_eq!(pair.reference_id(), 0);
    assert_eq!(pair.reference_start(), 10);
    assert_eq!(pair.mapping_quality(), 60);
    assert_eq!(pair.cigarstring().as_deref(), Some("50M"));
    assert_eq!(pair.query_length(), 50);
    assert_eq!(pair.tags().get("NM"), Some(&samrust_core::TagValue::Int(0)));
    assert_eq!(
        pair.tags().get("RG"),
        Some(&samrust_core::TagValue::Str("synth".into()))
    );

    assert_eq!(by_name("dup1").flag(), 1024);
    assert!(by_name("dup1").is_duplicate());
    assert_eq!(by_name("sec1").flag(), 256);
    assert!(by_name("sec1").is_secondary());
    assert_eq!(by_name("sup1").flag(), 2048);
    assert!(by_name("sup1").is_supplementary());
    assert_eq!(by_name("qcfail1").flag(), 512);
    assert!(by_name("qcfail1").is_qcfail());

    let unmapped = by_name("unmap1");
    assert_eq!(unmapped.flag(), 4);
    assert_eq!(unmapped.reference_id(), -1);
    assert_eq!(unmapped.reference_start(), -1);
    assert!(unmapped.cigarstring().is_none());
    assert_eq!(unmapped.query_length(), 10);

    let placed = by_name("unmap_placed");
    assert_eq!(placed.flag(), 133);
    assert!(placed.is_unmapped());
    assert_eq!(placed.reference_name(), Some("chr1"));
    assert_eq!(placed.reference_start(), 50);
    assert!(placed.cigarstring().is_none());

    assert_eq!(by_name("eqx1").cigarstring().as_deref(), Some("20=2X18="));
    assert_eq!(by_name("clip1").cigarstring().as_deref(), Some("5S20M5H"));
    assert_eq!(by_name("mapq255").mapping_quality(), 255);
}
