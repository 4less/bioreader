//! Tests for the deterministic read-id plumbing (`*_ids` paired reader).
//!
//! The load-bearing property is that `read_id` is a pure function of the input byte stream:
//! identical regardless of thread count and of `buffer_size`. These tests pin that down, plus
//! that the byte offsets actually point at the records they claim to, and that the id variant
//! visits exactly the same records as the plain variant.

use std::io::Write;

use bioreader::fastq_byte_reader::ReadId;
use bioreader::parallel::fastq::{
    read_fastq_paired_end_state_par, read_fastq_paired_end_state_par_ids, Merge,
};
use bioreader::sequence::fastq_record::RefFastqRecord;

/// One well-formed FASTQ record with a deliberately varying read length, so record byte
/// offsets are non-uniform and a naive "offset = i * fixed_len" would be wrong.
fn push_record(buf: &mut Vec<u8>, i: usize, mate: u8) {
    let len = 30 + (i * 7) % 50; // 30..=79, varies per record
    let seq: Vec<u8> = (0..len).map(|j| b"ACGT"[(i + j) % 4]).collect();
    let qual = vec![b'I'; len];
    write!(buf, "@read{}/{}\n", i, mate).unwrap();
    buf.extend_from_slice(&seq);
    buf.push(b'\n');
    buf.extend_from_slice(b"+\n");
    buf.extend_from_slice(&qual);
    buf.push(b'\n');
}

/// Build a paired FASTQ (R1, R2) with `n` records each.
fn make_pair(n: usize) -> (Vec<u8>, Vec<u8>) {
    let mut d1 = Vec::new();
    let mut d2 = Vec::new();
    for i in 0..n {
        push_record(&mut d1, i, 1);
        push_record(&mut d2, i, 2);
    }
    (d1, d2)
}

/// The header (without the leading `@` and trailing newline) of the record whose first byte is
/// at `off` in `data`. Panics if `off` does not point at a record start.
fn header_at(data: &[u8], off: u64) -> &[u8] {
    let off = off as usize;
    assert_eq!(data[off], b'@', "byte offset {} is not a record start", off);
    let nl = data[off + 1..]
        .iter()
        .position(|&c| c == b'\n')
        .expect("record header has no newline")
        + off
        + 1;
    &data[off + 1..nl]
}

type Row = (u64, u64, u64, Vec<u8>, Vec<u8>); // (ordinal, byte1, byte2, head1, head2)

#[derive(Default, Clone)]
struct Collector {
    rows: Vec<Row>,
}

impl Merge for Collector {
    fn merge_from(&mut self, other: &mut Self) {
        self.rows.append(&mut other.rows);
    }
}

/// Run the id variant and return rows sorted by ordinal (workers merge out of order).
fn collect_ids(d1: &[u8], d2: &[u8], buffer_size: usize, threads: u32) -> Vec<Row> {
    let st = read_fastq_paired_end_state_par_ids(
        d1,
        d2,
        buffer_size,
        threads,
        Collector::default(),
        |r1: &RefFastqRecord, r2: &RefFastqRecord, id: ReadId, st: &mut Collector| {
            st.rows.push((
                id.ordinal,
                id.byte_offset1,
                id.byte_offset2,
                r1.head().to_vec(),
                r2.head().to_vec(),
            ));
        },
    )
    .unwrap();
    let mut rows = st.rows;
    rows.sort_by_key(|r| r.0);
    rows
}

#[test]
fn ids_are_identical_across_thread_counts() {
    let (d1, d2) = make_pair(5000);
    let single = collect_ids(&d1, &d2, 512, 1);
    let many = collect_ids(&d1, &d2, 512, 16);

    assert_eq!(single.len(), 5000);
    assert_eq!(single, many, "id->record map differs between 1 and 16 threads");
}

#[test]
fn ids_are_identical_across_buffer_sizes() {
    // Byte offset and ordinal are intrinsic to the file, so changing how it is chunked must
    // not move them. This is the property (batch_id, idx) could NOT provide.
    let (d1, d2) = make_pair(5000);
    let small = collect_ids(&d1, &d2, 512, 8);
    let large = collect_ids(&d1, &d2, 65536, 8);

    assert_eq!(small, large, "id->record map differs between buffer sizes");
}

#[test]
fn ordinals_are_dense_and_in_file_order() {
    let (d1, d2) = make_pair(3000);
    let rows = collect_ids(&d1, &d2, 1024, 8);

    for (i, row) in rows.iter().enumerate() {
        assert_eq!(row.0, i as u64, "ordinal is not the dense file-order index");
        let expected1 = format!("read{}/1", i);
        let expected2 = format!("read{}/2", i);
        assert_eq!(row.3, expected1.as_bytes(), "ordinal {} maps to wrong R1", i);
        assert_eq!(row.4, expected2.as_bytes(), "ordinal {} maps to wrong R2", i);
    }
}

#[test]
fn byte_offsets_point_at_their_records() {
    let (d1, d2) = make_pair(4000);
    let rows = collect_ids(&d1, &d2, 777, 8); // odd buffer size on purpose

    for (ordinal, b1, b2, head1, head2) in &rows {
        assert_eq!(
            header_at(&d1, *b1),
            &head1[..],
            "R1 byte offset for ordinal {} points at the wrong record",
            ordinal
        );
        assert_eq!(
            header_at(&d2, *b2),
            &head2[..],
            "R2 byte offset for ordinal {} points at the wrong record",
            ordinal
        );
    }
}

#[derive(Default, Clone)]
struct HeadSet {
    heads: Vec<Vec<u8>>,
}

impl Merge for HeadSet {
    fn merge_from(&mut self, other: &mut Self) {
        self.heads.append(&mut other.heads);
    }
}

#[test]
fn id_variant_visits_same_records_as_plain_variant() {
    let (d1, d2) = make_pair(5000);

    let plain = read_fastq_paired_end_state_par(
        d1.as_slice(),
        d2.as_slice(),
        1024,
        8,
        HeadSet::default(),
        |r1: &RefFastqRecord, _r2: &RefFastqRecord, st: &mut HeadSet| {
            st.heads.push(r1.head().to_vec());
        },
    )
    .unwrap();

    let mut plain_heads = plain.heads;
    plain_heads.sort();

    let mut id_heads: Vec<Vec<u8>> = collect_ids(&d1, &d2, 1024, 8)
        .into_iter()
        .map(|r| r.3)
        .collect();
    id_heads.sort();

    assert_eq!(plain_heads.len(), 5000);
    assert_eq!(plain_heads, id_heads, "id variant saw a different set of records");
}
