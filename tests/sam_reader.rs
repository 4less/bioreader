//! Tests for the multithreaded SAM reader.
//!
//! Two load-bearing properties: every worker together sees each alignment line exactly once
//! regardless of thread count and `buffer_size`, and `SamReadId` is a pure function of the
//! input byte stream (same ordinal and same absolute byte offset no matter how the file was
//! chunked or which thread happened to parse it).

use std::io::Write;
use std::sync::{Arc, Mutex};

use bioreader::parallel::fastq::Merge;
use bioreader::parallel::sam::{
    read_sam_par, read_sam_state_par, read_sam_state_par_ids, read_sam_state_par_shared,
    read_sam_state_par_with_header,
};
use bioreader::sam_byte_reader::{SamByteReader, SamReadId};
use bioreader::sequence::sam_record::RefSamRecord;

const REFERENCES: usize = 3;

/// One well-formed alignment line with a deliberately varying read length and a varying number
/// of optional fields, so record byte offsets are non-uniform and a naive
/// "offset = i * fixed_len" would be wrong.
fn push_record(buf: &mut Vec<u8>, i: usize) {
    let len = 30 + (i * 7) % 50; // 30..=79, varies per record
    let seq: Vec<u8> = (0..len).map(|j| b"ACGT"[(i + j) % 4]).collect();
    let qual = vec![b'I'; len];
    let flag = if i % 2 == 0 { 99 } else { 147 };
    let reference = format!("chr{}", i % REFERENCES + 1);

    write!(
        buf,
        "read{}\t{}\t{}\t{}\t60\t{}M\t=\t{}\t{}\t",
        i,
        flag,
        reference,
        i + 1,
        len,
        i + 200,
        len * 2
    )
    .unwrap();
    buf.extend_from_slice(&seq);
    buf.push(b'\t');
    buf.extend_from_slice(&qual);
    // Every third record carries optional fields, the rest stop after QUAL.
    if i % 3 == 0 {
        write!(buf, "\tNM:i:{}\tMD:Z:{}\tAS:i:-{}", i % 5, len, i % 7).unwrap();
    }
    buf.push(b'\n');
}

/// Build a SAM file with a header and `n` alignment records.
fn make_sam(n: usize) -> Vec<u8> {
    let mut sam = Vec::new();
    sam.extend_from_slice(b"@HD\tVN:1.6\tSO:unsorted\n");
    for r in 0..REFERENCES {
        write!(sam, "@SQ\tSN:chr{}\tLN:{}\n", r + 1, 100_000 * (r + 1)).unwrap();
    }
    sam.extend_from_slice(b"@RG\tID:sample1\tSM:sample1\n");
    sam.extend_from_slice(b"@PG\tID:bioreader\tPN:bioreader\n");
    for i in 0..n {
        push_record(&mut sam, i);
    }
    sam
}

/// The QNAME of the record whose first byte is at `off` in `data`. Panics if `off` does not
/// point at the start of a line.
fn qname_at(data: &[u8], off: u64) -> &[u8] {
    let off = off as usize;
    assert!(
        off == 0 || data[off - 1] == b'\n',
        "byte offset {} is not a record start",
        off
    );
    let tab = data[off..]
        .iter()
        .position(|&c| c == b'\t')
        .expect("record has no tab")
        + off;
    &data[off..tab]
}

type Row = (u64, u64, Vec<u8>, u64); // (ordinal, byte_offset, qname, pos)

#[derive(Default, Clone)]
struct Collector {
    rows: Vec<Row>,
}

impl Merge for Collector {
    fn merge_from(&mut self, other: &mut Self) {
        self.rows.append(&mut other.rows);
    }
}

#[derive(Default, Clone)]
struct Names {
    names: Vec<Vec<u8>>,
}

impl Merge for Names {
    fn merge_from(&mut self, other: &mut Self) {
        self.names.append(&mut other.names);
    }
}

/// Run the id variant and return rows sorted by ordinal (workers merge out of order).
fn collect_ids(sam: &[u8], buffer_size: usize, threads: u32) -> Vec<Row> {
    let state = read_sam_state_par_ids(
        sam,
        buffer_size,
        threads,
        Collector::default(),
        |record: &RefSamRecord, id: SamReadId, state: &mut Collector| {
            state
                .rows
                .push((id.ordinal, id.byte_offset, record.qname().to_vec(), record.pos()));
        },
    )
    .unwrap();

    let mut rows = state.rows;
    rows.sort_by_key(|row| row.0);
    rows
}

/// Run the plain variant and return the QNAMEs in sorted order.
fn collect_names(sam: &[u8], buffer_size: usize, threads: u32) -> Vec<Vec<u8>> {
    let state = read_sam_state_par(
        sam,
        buffer_size,
        threads,
        Names::default(),
        |record: &RefSamRecord, state: &mut Names| state.names.push(record.qname().to_vec()),
    )
    .unwrap();

    let mut names = state.names;
    names.sort();
    names
}

#[test]
fn every_record_is_read_exactly_once() {
    let sam = make_sam(5000);
    let expected: Vec<Vec<u8>> = {
        let mut names: Vec<Vec<u8>> = (0..5000).map(|i| format!("read{}", i).into_bytes()).collect();
        names.sort();
        names
    };

    for (buffer_size, threads) in [(256, 1), (256, 8), (1 << 16, 8), (1 << 16, 1)] {
        assert_eq!(
            collect_names(&sam, buffer_size, threads),
            expected,
            "buffer_size {} / {} threads lost or duplicated records",
            buffer_size,
            threads
        );
    }
}

#[test]
fn ids_are_identical_across_thread_counts() {
    let sam = make_sam(5000);
    let single = collect_ids(&sam, 512, 1);
    let many = collect_ids(&sam, 512, 16);

    assert_eq!(single.len(), 5000);
    assert_eq!(single, many, "id->record map differs between 1 and 16 threads");
}

#[test]
fn ids_are_identical_across_buffer_sizes() {
    // Byte offset and ordinal are intrinsic to the file, so changing how it is chunked must not
    // move them.
    let sam = make_sam(5000);
    let small = collect_ids(&sam, 512, 8);
    let large = collect_ids(&sam, 65536, 8);

    assert_eq!(small, large, "id->record map differs between buffer sizes");
}

#[test]
fn ordinals_are_dense_and_in_file_order() {
    let sam = make_sam(3000);
    let rows = collect_ids(&sam, 1024, 8);

    for (i, row) in rows.iter().enumerate() {
        assert_eq!(row.0, i as u64, "ordinal is not the dense file-order index");
        assert_eq!(
            row.2,
            format!("read{}", i).into_bytes(),
            "ordinal {} maps to the wrong record",
            i
        );
        assert_eq!(row.3, i as u64 + 1, "ordinal {} maps to the wrong POS", i);
    }
}

#[test]
fn byte_offsets_point_at_their_records() {
    let sam = make_sam(4000);
    let rows = collect_ids(&sam, 777, 8); // odd buffer size on purpose

    for (ordinal, byte_offset, qname, _pos) in &rows {
        assert_eq!(
            qname_at(&sam, *byte_offset),
            &qname[..],
            "byte offset for ordinal {} points at the wrong record",
            ordinal
        );
    }
}

#[test]
fn id_variant_visits_same_records_as_plain_variant() {
    let sam = make_sam(5000);

    let plain = collect_names(&sam, 1024, 8);
    let mut with_ids: Vec<Vec<u8>> = collect_ids(&sam, 1024, 8).into_iter().map(|r| r.2).collect();
    with_ids.sort();

    assert_eq!(plain.len(), 5000);
    assert_eq!(plain, with_ids, "id variant saw a different set of records");
}

#[test]
fn header_is_available_alongside_the_records() {
    let sam = make_sam(1000);

    let (header, state) = read_sam_state_par_with_header(
        sam.as_slice(),
        4096,
        4,
        Names::default(),
        |record: &RefSamRecord, state: &mut Names| state.names.push(record.qname().to_vec()),
    )
    .unwrap();

    assert_eq!(state.names.len(), 1000);
    assert_eq!(header.lines().count(), 3 + REFERENCES); // @HD, @RG, @PG and one @SQ per reference
    assert_eq!(
        header.reference_sequences(),
        vec![
            (b"chr1".as_slice(), 100_000),
            (b"chr2".as_slice(), 200_000),
            (b"chr3".as_slice(), 300_000),
        ]
    );
}

#[test]
fn shared_reader_exposes_the_header_before_processing() {
    let sam = make_sam(1000);
    let byte_reader = Arc::new(Mutex::new(SamByteReader::new(sam.as_slice(), 4096).unwrap()));

    let references = byte_reader.lock().unwrap().header().reference_sequences().len();
    assert_eq!(references, REFERENCES);

    let state = read_sam_state_par_shared(
        &byte_reader,
        4096,
        4,
        Names::default(),
        |record: &RefSamRecord, state: &mut Names| state.names.push(record.qname().to_vec()),
    );

    assert_eq!(state.names.len(), 1000);
}

#[test]
fn records_and_optional_fields_parse_under_every_chunking() {
    let sam = make_sam(2000);

    for buffer_size in [128, 512, 4096, 1 << 16] {
        let tagged = read_sam_par(sam.as_slice(), buffer_size, 4, |record: &RefSamRecord| {
            assert!(record.valid_extended(), "invalid record {}", record);
            assert!(record.is_paired() && record.is_proper_pair());
            assert_eq!(record.rnext(), b"=");
            match record.tag(b"NM") {
                Some(tag) => {
                    assert!(tag.as_int().is_some());
                    assert_eq!(record.tags().count(), 3);
                    1
                }
                None => {
                    assert!(record.optional().is_empty());
                    0
                }
            }
        })
        .unwrap();

        // Every third record carries the optional fields.
        assert_eq!(tagged, 2000usize.div_ceil(3), "buffer_size {}", buffer_size);
    }
}

#[test]
fn reads_gzipped_input() {
    use flate2::write::GzEncoder;
    use flate2::Compression;

    let sam = make_sam(3000);
    let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
    encoder.write_all(&sam).unwrap();
    let gz = encoder.finish().unwrap();

    let rows = read_sam_state_par_ids(
        flate2::read::GzDecoder::new(gz.as_slice()),
        4096,
        4,
        Collector::default(),
        |record: &RefSamRecord, id: SamReadId, state: &mut Collector| {
            state
                .rows
                .push((id.ordinal, id.byte_offset, record.qname().to_vec(), record.pos()));
        },
    )
    .unwrap();

    let mut rows = rows.rows;
    rows.sort_by_key(|row| row.0);

    // Offsets are into the decompressed stream, so they must match the plain-text run exactly.
    assert_eq!(rows, collect_ids(&sam, 4096, 4));
}

/// A file whose records are far larger than the buffer must still be cut on record boundaries.
#[test]
fn records_larger_than_the_buffer() {
    let sam = make_sam(200);
    let rows = collect_ids(&sam, 16, 4);

    assert_eq!(rows.len(), 200);
    assert_eq!(rows, collect_ids(&sam, 1 << 16, 4));
}

#[test]
fn empty_body_yields_no_records() {
    let mut sam = Vec::new();
    sam.extend_from_slice(b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n");

    let (header, state) = read_sam_state_par_with_header(
        sam.as_slice(),
        1024,
        4,
        Names::default(),
        |record: &RefSamRecord, state: &mut Names| state.names.push(record.qname().to_vec()),
    )
    .unwrap();

    assert!(state.names.is_empty());
    assert_eq!(header.reference_sequences(), vec![(b"chr1".as_slice(), 1000)]);
}

#[test]
fn headerless_input_is_read_as_records() {
    let mut sam = Vec::new();
    for i in 0..100 {
        push_record(&mut sam, i);
    }

    let rows = collect_ids(&sam, 256, 4);
    assert_eq!(rows.len(), 100);
    assert_eq!(rows[0].1, 0, "first record of a headerless file starts at byte 0");
    for (ordinal, byte_offset, qname, _pos) in &rows {
        assert_eq!(qname_at(&sam, *byte_offset), &qname[..], "ordinal {}", ordinal);
    }
}
