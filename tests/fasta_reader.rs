//! Tests for the FASTA reader: soft-masked (lowercase) bases come out uppercase, headers as they
//! are, and multi-line sequences are joined.

use std::io::Cursor;
use std::sync::{Arc, Mutex};

use bioreader::fasta_byte_reader::FastaByteReader;
use bioreader::fasta_reader::FastaReader;
use bioreader::sequence::fasta_record::OwnedFastaRecord;

fn records(fasta: &str) -> Vec<(String, String)> {
    let mut bytes = Arc::new(Mutex::new(FastaByteReader::new(Cursor::new(fasta.as_bytes().to_vec()), 1 << 16).unwrap()));
    let mut reader = FastaReader::with_capacity(1 << 16);
    let mut record = OwnedFastaRecord::new();
    let mut out = Vec::new();
    while let Some(()) = reader.load_batch_par(&mut bytes).unwrap() {
        while let Some(()) = reader.next(&mut record) {
            out.push((String::from_utf8(record.head().to_vec()).unwrap(), String::from_utf8(record.seq().to_vec()).unwrap()));
        }
    }
    out
}

#[test]
fn soft_masked_bases_are_uppercased() {
    let got = records(">chr1A Triticum aestivum\nctaaacccTAAACCnnNNacgt\nAACCggtt\n>c2 lower header\nACGT\n");
    assert_eq!(got, vec![
        (">chr1A Triticum aestivum".to_string(), "CTAAACCCTAAACCNNNNACGTAACCGGTT".to_string()),
        (">c2 lower header".to_string(), "ACGT".to_string()),
    ]);
}

#[test]
fn uppercased_records_pass_valid_extended() {
    let mut bytes = Arc::new(Mutex::new(FastaByteReader::new(Cursor::new(b">r\nacgtrykn\n>s\nACGT\n".to_vec()), 1 << 16).unwrap()));
    let mut reader = FastaReader::with_capacity(1 << 16);
    let mut record = OwnedFastaRecord::new();
    let mut seen = 0;
    while let Some(()) = reader.load_batch_par(&mut bytes).unwrap() {
        while let Some(()) = reader.next(&mut record) {
            assert!(record.valid_extended(), "{:?}", record.to_string());
            seen += 1;
        }
    }
    assert_eq!(seen, 2);
}
