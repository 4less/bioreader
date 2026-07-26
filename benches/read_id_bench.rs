//! Head-to-head timing of the plain paired reader vs. the read-id variant.
//!
//! Both benches do the *same* per-record work (sum sequence lengths) over the same in-memory
//! FASTQ with the same thread count and buffer size. The only difference is that `with_ids`
//! goes through `*_ids` (batch-base counters under the lock + per-record id derivation +
//! passing the id to the closure). The delta between the two is the cost of the id plumbing.
//!
//! Run: `cargo +nightly bench --bench read_id_bench`

#![feature(test)]
extern crate test;

use std::io::Write;

use bioreader::fastq_byte_reader::ReadId;
use bioreader::parallel::fastq::{
    read_fastq_paired_end_state_par, read_fastq_paired_end_state_par_ids, Merge,
};
use bioreader::sequence::fastq_record::RefFastqRecord;
use test::{black_box, Bencher};

const N_RECORDS: usize = 500_000;
const BUFFER_SIZE: usize = 1 << 16;
const THREADS: u32 = 4;

fn push_record(buf: &mut Vec<u8>, i: usize, mate: u8) {
    let len = 90 + (i * 7) % 60; // ~90..150, realistic short-read lengths
    let seq: Vec<u8> = (0..len).map(|j| b"ACGT"[(i + j) % 4]).collect();
    let qual = vec![b'I'; len];
    write!(buf, "@read{}/{}\n", i, mate).unwrap();
    buf.extend_from_slice(&seq);
    buf.push(b'\n');
    buf.extend_from_slice(b"+\n");
    buf.extend_from_slice(&qual);
    buf.push(b'\n');
}

fn make_pair(n: usize) -> (Vec<u8>, Vec<u8>) {
    let mut d1 = Vec::new();
    let mut d2 = Vec::new();
    for i in 0..n {
        push_record(&mut d1, i, 1);
        push_record(&mut d2, i, 2);
    }
    (d1, d2)
}

#[derive(Default, Clone)]
struct Counter(u64);

impl Merge for Counter {
    fn merge_from(&mut self, other: &mut Self) {
        self.0 = self.0.wrapping_add(other.0);
    }
}

#[bench]
fn no_ids(b: &mut Bencher) {
    let (d1, d2) = make_pair(N_RECORDS);
    b.iter(|| {
        let st = read_fastq_paired_end_state_par(
            d1.as_slice(),
            d2.as_slice(),
            BUFFER_SIZE,
            THREADS,
            Counter::default(),
            |r1: &RefFastqRecord, _r2: &RefFastqRecord, st: &mut Counter| {
                st.0 = st.0.wrapping_add(r1.seq().len() as u64);
            },
        )
        .unwrap();
        black_box(st.0)
    });
}

#[bench]
fn with_ids(b: &mut Bencher) {
    let (d1, d2) = make_pair(N_RECORDS);
    b.iter(|| {
        let st = read_fastq_paired_end_state_par_ids(
            d1.as_slice(),
            d2.as_slice(),
            BUFFER_SIZE,
            THREADS,
            Counter::default(),
            |r1: &RefFastqRecord, _r2: &RefFastqRecord, id: ReadId, st: &mut Counter| {
                st.0 = st.0.wrapping_add(r1.seq().len() as u64);
                black_box(id); // ensure the id is actually materialized, not optimized away
            },
        )
        .unwrap();
        black_box(st.0)
    });
}
