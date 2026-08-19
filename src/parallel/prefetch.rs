//! Batch prefetching for the parallel SAM reader.
//!
//! [`read_sam_state_par_shared`](super::sam::read_sam_state_par_shared) has every worker
//! call `load_batch_par`, which holds the shared reader's lock across the whole of
//! `cut_batch`: the read from the source, the gzip inflation when the source is
//! compressed, the scan back to the last newline, and the copy into the worker's buffer.
//! Byte *delivery* is therefore serial and only parsing runs in parallel, which on a SAM
//! whose records are cheap to reject is a small fraction of the work. Measured on rg-msa
//! over 6.3 GB of plain SAM: 10.39 s at one worker and 6.90 s at sixteen, never above
//! 1.9 cores, and unchanged by `--read-buffer` because the cost under the lock is linear
//! in bytes however it is chunked.
//!
//! Here one thread does the reading and cutting, and hands finished batches to workers
//! through a channel. A worker's only synchronisation is taking a batch off the queue --
//! a pointer swap, not a read.
//!
//! Memory is bounded by construction rather than by trust: `depth` buffers are allocated
//! up front and recycled, and the producer must take an empty one before it can fill it.
//! A slow consumer therefore blocks the producer instead of growing a queue.

use std::io::{Error, Read};
use std::sync::mpsc::{channel, Receiver, Sender};
use std::sync::Mutex;
use std::thread::{self, JoinHandle};

use crate::fastq_byte_reader::FillBuffer;
use crate::sam_byte_reader::SamByteReader;
use crate::sam_byte_reader::SamHeader;

/// A filled batch travelling from the producer to a worker, or the error that ended it.
type Delivery = Result<(Vec<u8>, usize), Error>;

/// Reads and cuts SAM batches on its own thread so workers never block on IO.
pub struct BatchPrefetcher {
    /// Behind a `Mutex` because `Receiver` is not `Sync`. The lock is held for a `recv`
    /// and nothing else, so unlike the reader lock it does not span the read.
    filled: Mutex<Receiver<Delivery>>,
    empty: Sender<Vec<u8>>,
    header: SamHeader,
    producer: Mutex<Option<JoinHandle<()>>>,
}

impl BatchPrefetcher {
    /// Start prefetching from `reader`, keeping at most `depth` batches in flight.
    ///
    /// `depth` must be at least one; two is the smallest that lets the producer work
    /// while a consumer parses, and more only helps when batch times are uneven.
    pub fn new<T: Read + Send + 'static>(
        mut reader: SamByteReader<T>,
        buffer_size: usize,
        depth: usize,
    ) -> Self {
        let header = reader.header().clone();
        let (filled_tx, filled_rx) = channel::<Delivery>();
        let (empty_tx, empty_rx) = channel::<Vec<u8>>();

        for _ in 0..depth.max(1) {
            let _ = empty_tx.send(vec![0u8; buffer_size]);
        }

        let producer = thread::spawn(move || {
            // Ends when the empty queue is dropped (every worker gone) or the source is
            // exhausted; dropping `filled_tx` on the way out is what tells workers so.
            while let Ok(mut buf) = empty_rx.recv() {
                match reader.fill_buf(&mut buf) {
                    Ok(Some(bytes)) => {
                        if filled_tx.send(Ok((buf, bytes))).is_err() {
                            return;
                        }
                    }
                    Ok(None) => return,
                    Err(e) => {
                        let _ = filled_tx.send(Err(e));
                        return;
                    }
                }
            }
        });

        Self {
            filled: Mutex::new(filled_rx),
            empty: empty_tx,
            header,
            producer: Mutex::new(Some(producer)),
        }
    }

    pub fn header(&self) -> &SamHeader {
        &self.header
    }

    /// Take the next batch into `buf`, returning its length in bytes.
    ///
    /// `buf` is swapped with the delivered buffer rather than copied into, and the old one
    /// goes back to the producer to be refilled. `Ok(None)` means the source is done.
    pub fn next_batch(&self, buf: &mut Vec<u8>) -> Result<Option<usize>, Error> {
        let msg = {
            let rx = self
                .filled
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner());
            rx.recv()
        };
        match msg {
            // The producer dropped its sender: the source ran out or it returned early.
            Err(_) => Ok(None),
            Ok(Err(e)) => Err(e),
            Ok(Ok((batch, bytes))) => {
                let spent = std::mem::replace(buf, batch);
                // Failing to recycle is not an error: it only means the producer has
                // already stopped, and the next `recv` will report that.
                let _ = self.empty.send(spent);
                Ok(Some(bytes))
            }
        }
    }

    /// Stop the producer and wait for it, so a caller that abandons a file early does not
    /// leave a thread reading it.
    pub fn finish(&self) {
        let handle = self
            .producer
            .lock()
            .unwrap_or_else(|p| p.into_inner())
            .take();
        if let Some(h) = handle {
            let _ = h.join();
        }
    }
}

impl Drop for BatchPrefetcher {
    fn drop(&mut self) {
        self.finish();
    }
}

/// Like [`read_sam_state_par_shared`](super::sam::read_sam_state_par_shared), but the
/// batches come off a [`BatchPrefetcher`] instead of from the reader lock.
///
/// Same contract: one `State` per worker, merged at the end, so a caller sees exactly the
/// records it would have seen either way. The difference is only where the reading happens.
pub fn read_sam_state_par_prefetch<G, T, State>(
    prefetcher: &BatchPrefetcher,
    buffer_size: usize,
    num_threads: u32,
    mut global_state: State,
    f: G,
) -> State
where
    G: FnMut(&crate::sequence::sam_record::RefSamRecord, &mut State) + Clone + Send,
    T: Read + Send,
    State: Default + Clone + Send + crate::parallel::fastq::Merge,
{
    let _ = std::marker::PhantomData::<T>;
    std::thread::scope(|scope| {
        let mut threads = Vec::new();
        for _ in 0..num_threads.max(1) {
            let mut sam_reader = crate::sam_reader::SamReader::with_capacity(buffer_size);
            let mut state = global_state.clone();
            let mut f_local = f.clone();
            threads.push(scope.spawn(move || {
                loop {
                    match prefetcher.next_batch(&mut sam_reader.buffer) {
                        Ok(Some(bytes)) => {
                            sam_reader.set_batch(bytes);
                            while let Some(record) = sam_reader.next() {
                                f_local(&record, &mut state);
                            }
                        }
                        Ok(None) => break,
                        Err(_) => break,
                    }
                }
                state
            }));
        }
        for t in threads {
            let mut state = t.join().expect("worker panicked");
            global_state.merge_from(&mut state);
        }
    });
    global_state
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parallel::fastq::Merge;
    use crate::parallel::sam::read_sam_state_par_shared;
    use crate::sequence::sam_record::RefSamRecord;
    use std::sync::Arc;

    #[derive(Default, Clone)]
    struct Seen {
        names: Vec<String>,
        bases: usize,
    }
    impl Merge for Seen {
        fn merge_from(&mut self, other: &mut Self) {
            self.names.append(&mut other.names);
            self.bases += other.bases;
        }
    }

    fn sam(n: usize) -> Vec<u8> {
        let mut s = String::from("@HD\tVN:1.6\n@SQ\tSN:ref1\tLN:1000\n");
        for i in 0..n {
            // A record wide enough that a small batch buffer cuts mid-stream.
            s.push_str(&format!(
                "read{i}\t0\tref1\t{}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n",
                i % 900 + 1
            ));
        }
        s.into_bytes()
    }

    fn via_lock(data: &[u8], buffer_size: usize, threads: u32) -> Seen {
        let br = Arc::new(Mutex::new(
            SamByteReader::new(std::io::Cursor::new(data.to_vec()), buffer_size).unwrap(),
        ));
        read_sam_state_par_shared(
            &br,
            buffer_size,
            threads,
            Seen::default(),
            |r: &RefSamRecord, st: &mut Seen| {
                st.names.push(String::from_utf8_lossy(r.qname()).into_owned());
                st.bases += r.seq().len();
            },
        )
    }

    fn via_prefetch(data: &[u8], buffer_size: usize, threads: u32) -> Seen {
        let reader =
            SamByteReader::new(std::io::Cursor::new(data.to_vec()), buffer_size).unwrap();
        let pf = BatchPrefetcher::new(reader, buffer_size, 4);
        read_sam_state_par_prefetch::<_, std::io::Cursor<Vec<u8>>, _>(
            &pf,
            buffer_size,
            threads,
            Seen::default(),
            |r: &RefSamRecord, st: &mut Seen| {
                st.names.push(String::from_utf8_lossy(r.qname()).into_owned());
                st.bases += r.seq().len();
            },
        )
    }

    /// The property the whole change rests on: moving the read off the lock must not
    /// change which records a caller sees, only when they arrive.
    #[test]
    fn prefetching_yields_exactly_what_the_locked_reader_yields() {
        for &(n, buf, threads) in &[
            (1_000usize, 1 << 12, 1u32),
            (1_000, 1 << 12, 4),
            (10_000, 1 << 14, 8),
            (10_000, 1 << 20, 8), // one batch, more workers than work
        ] {
            let data = sam(n);
            let mut a = via_lock(&data, buf, threads);
            let mut b = via_prefetch(&data, buf, threads);
            assert_eq!(a.bases, b.bases, "{n} records, {buf} buffer, {threads} threads");
            a.names.sort();
            b.names.sort();
            assert_eq!(a.names, b.names, "{n} records, {buf} buffer, {threads} threads");
            assert_eq!(a.names.len(), n);
        }
    }

    /// A record longer than the batch buffer forces `cut_batch` to grow it. The producer
    /// runs that path on its own thread, so the growth must survive the handoff.
    #[test]
    fn a_record_longer_than_the_buffer_still_arrives_whole() {
        let mut s = String::from("@HD\tVN:1.6\n@SQ\tSN:ref1\tLN:100000\n");
        let seq = "A".repeat(50_000);
        s.push_str(&format!("long\t0\tref1\t1\t60\t50000M\t*\t0\t0\t{seq}\t{}\n", "I".repeat(50_000)));
        s.push_str("short\t0\tref1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
        let data = s.into_bytes();
        let got = via_prefetch(&data, 1 << 10, 4);
        assert_eq!(got.names.len(), 2);
        assert_eq!(got.bases, 50_004);
    }

    /// No trailing newline: the last record must not be dropped or duplicated.
    #[test]
    fn a_file_without_a_final_newline_loses_nothing() {
        let mut data = sam(50);
        assert_eq!(data.pop(), Some(b'\n'));
        let a = via_lock(&data, 1 << 12, 4);
        let b = via_prefetch(&data, 1 << 12, 4);
        assert_eq!(a.names.len(), 50);
        assert_eq!(b.names.len(), 50);
        assert_eq!(a.bases, b.bases);
    }

    /// Header but no alignments: workers must finish rather than block on an empty queue.
    #[test]
    fn an_alignment_free_file_terminates() {
        let data = b"@HD\tVN:1.6\n@SQ\tSN:ref1\tLN:10\n".to_vec();
        assert_eq!(via_prefetch(&data, 1 << 12, 4).names.len(), 0);
    }

    /// Buffers are recycled, so the number in flight is what was allocated -- a long file
    /// through a shallow queue must not grow the process.
    #[test]
    fn a_shallow_queue_still_reads_a_long_file() {
        let data = sam(20_000);
        let reader =
            SamByteReader::new(std::io::Cursor::new(data.clone()), 1 << 12).unwrap();
        let pf = BatchPrefetcher::new(reader, 1 << 12, 1);
        let got = read_sam_state_par_prefetch::<_, std::io::Cursor<Vec<u8>>, _>(
            &pf, 1 << 12, 8, Seen::default(),
            |r: &RefSamRecord, st: &mut Seen| {
                st.names.push(String::from_utf8_lossy(r.qname()).into_owned());
                st.bases += r.seq().len();
            },
        );
        assert_eq!(got.names.len(), 20_000);
    }
}
