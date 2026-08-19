// How fast can the prefetcher alone deliver batches, with consumers that do nothing?
// If this matches rg-msa's throughput, the workers are starved and the producer is the
// wall; if it is far faster, the workers are the wall.
use bioreader::parallel::prefetch::{read_sam_state_par_prefetch, BatchPrefetcher};
use bioreader::parallel::fastq::Merge;
use bioreader::sam_byte_reader::SamByteReader;
use bioreader::sequence::sam_record::RefSamRecord;
use std::fs::File;
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
        let f = File::open(p).unwrap();
        let br = SamByteReader::new(f, bufsz).unwrap();
        let pf = BatchPrefetcher::new(br, bufsz, 4);
        let c = read_sam_state_par_prefetch::<_, File, _>(
            &pf, bufsz, threads, Count::default(),
            // Touch one byte so the record is not optimised away, and nothing more.
            |r: &RefSamRecord, st: &mut Count| st.0 += r.qname().len() as u64,
        );
        total += c.0;
    }
    let secs = t0.elapsed().as_secs_f64();
    eprintln!("{threads} workers: {secs:.2} s  (qname bytes {total})");
}
