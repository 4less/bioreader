// How fast can the prefetcher alone deliver batches, with consumers that do nothing?
// If this matches rg-msa's throughput, the workers are starved and the producer is the
// wall; if it is far faster, the workers are the wall.
use bioreader::parallel::prefetch::{read_sam_state_par_prefetch, BatchPrefetcher};
use bioreader::parallel::fastq::Merge;
use bioreader::sam_byte_reader::SamByteReader;
use bioreader::sequence::sam_record::RefSamRecord;
use std::fs::File;
use std::io::Read;
use std::time::Instant;

#[derive(Default, Clone)]
struct Count(u64);
impl Merge for Count { fn merge_from(&mut self, o: &mut Self) { self.0 += o.0; } }

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let threads: u32 = a[1].parse().unwrap();
    let bufsz: usize = a[2].parse().unwrap();
    let mut total = 0u64;
    let t0 = Instant::now();
    for p in &a[3..] {
        // gz goes through the same parallel decoder rg-msa uses, so this measures the
        // producer *plus* the decoder with no handler behind it.
        let src: Box<dyn Read + Send> = if p.ends_with(".gz") {
            let d = rapidgzip_core::Decoder::builder().decoder_threads(16).build().unwrap();
            Box::new(d.open(p).unwrap())
        } else {
            Box::new(File::open(p).unwrap())
        };
        let br = SamByteReader::new(src, bufsz).unwrap();
        let pf = BatchPrefetcher::new(br, bufsz, 4);
        let c = read_sam_state_par_prefetch::<_, Box<dyn Read + Send>, _>(
            &pf, bufsz, threads, Count::default(),
            // Touch one byte so the record is not optimised away, and nothing more.
            |r: &RefSamRecord, st: &mut Count| {
                // Simulate a filter's per-record cost: touch the fields it would and do
                // a little arithmetic, so the consumer is heavy rather than free.
                let mut h = 0u64;
                for b in r.qname() { h = h.wrapping_mul(31).wrapping_add(*b as u64); }
                for b in r.cigar() { h = h.wrapping_mul(31).wrapping_add(*b as u64); }
                for b in r.seq()   { h = h.wrapping_mul(31).wrapping_add(*b as u64); }
                for b in r.qual()  { h = h.wrapping_mul(31).wrapping_add(*b as u64); }
                st.0 = st.0.wrapping_add(h);
            },
        );
        total += c.0;
    }
    let secs = t0.elapsed().as_secs_f64();
    eprintln!("{threads} workers: {secs:.2} s  (qname bytes {total})");
}
