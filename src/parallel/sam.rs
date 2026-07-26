//! Multithreaded SAM drivers.
//!
//! Same pattern as [`crate::parallel::fastq`]: one [`SamByteReader`] behind an `Arc<Mutex<..>>`
//! hands out batches of complete records, and `num_threads` scoped workers each parse their
//! batches with a private [`SamReader`]. The scope guarantees all threads finish inside the
//! call, so neither the closure nor the input stream needs `'static`.

use std::sync::{Arc, Mutex};

use crate::{
    parallel::fastq::Merge,
    sam_byte_reader::{SamByteReader, SamHeader, SamReadId},
    sam_reader::SamReader,
    sequence::sam_record::RefSamRecord,
};

/// Runs `f` over every alignment record and returns the sum of its results.
pub fn read_sam_par<G, T>(
    file: T,
    buffer_size: usize,
    num_threads: u32,
    f: G,
) -> Result<usize, std::io::Error>
where
    G: Fn(&RefSamRecord) -> usize + Clone + Send,
    T: std::io::Read + std::marker::Send,
{
    let byte_reader = Arc::new(Mutex::new(SamByteReader::new(file, buffer_size)?));
    let mut total = 0;

    std::thread::scope(|scope| {
        let mut threads = Vec::new();

        for _thread in 0..num_threads {
            let mut sam_reader = SamReader::with_capacity(buffer_size);
            let mut reader_local = byte_reader.clone(); // Need to clone for move
            let f_local = f.clone();

            threads.push(scope.spawn(move || {
                let mut result = 0;
                while let Some(()) = sam_reader
                    .load_batch_par(&mut reader_local)
                    .expect("Batch is invalid")
                {
                    while let Some(record) = sam_reader.next() {
                        result += f_local(&record);
                    }
                }
                result
            }));
        }

        for thread_guard in threads {
            total += thread_guard.join().expect("Worker thread panicked");
        }
    });

    Ok(total)
}

/// Runs `f` over every alignment record with a per-thread `State`, merging the thread states
/// into `global_state` at the end. Counterpart of
/// [`read_fastq_paired_end_state_par`](crate::parallel::fastq::read_fastq_paired_end_state_par).
pub fn read_sam_state_par<G, T, State>(
    file: T,
    buffer_size: usize,
    num_threads: u32,
    global_state: State,
    f: G,
) -> Result<State, std::io::Error>
where
    G: FnMut(&RefSamRecord, &mut State) -> () + Clone + Send,
    T: std::io::Read + std::marker::Send,
    State: Default + Clone + Send + Merge,
{
    let byte_reader = Arc::new(Mutex::new(SamByteReader::new(file, buffer_size)?));

    Ok(read_sam_state_par_shared(
        &byte_reader,
        buffer_size,
        num_threads,
        global_state,
        f,
    ))
}

/// Like [`read_sam_state_par`], but also returns the file's header. The header is consumed by
/// the byte reader before any worker starts, so this costs nothing extra.
pub fn read_sam_state_par_with_header<G, T, State>(
    file: T,
    buffer_size: usize,
    num_threads: u32,
    global_state: State,
    f: G,
) -> Result<(SamHeader, State), std::io::Error>
where
    G: FnMut(&RefSamRecord, &mut State) -> () + Clone + Send,
    T: std::io::Read + std::marker::Send,
    State: Default + Clone + Send + Merge,
{
    let byte_reader = Arc::new(Mutex::new(SamByteReader::new(file, buffer_size)?));

    let state = read_sam_state_par_shared(&byte_reader, buffer_size, num_threads, global_state, f);
    let header = byte_reader
        .lock()
        .expect("Locking ByteReader was unsuccessful")
        .header()
        .clone();

    Ok((header, state))
}

/// The driver behind [`read_sam_state_par`], taking an already constructed byte reader. Use it
/// when the header has to be inspected before deciding how to process the records:
///
/// ```no_run
/// # use std::sync::{Arc, Mutex};
/// # use bioreader::sam_byte_reader::SamByteReader;
/// # use bioreader::parallel::sam::read_sam_state_par_shared;
/// # use bioreader::parallel::fastq::Merge;
/// # use bioreader::sequence::sam_record::RefSamRecord;
/// # #[derive(Default, Clone)] struct Count(usize);
/// # impl Merge for Count { fn merge_from(&mut self, other: &mut Self) { self.0 += other.0; } }
/// # fn main() -> std::io::Result<()> {
/// let file = std::fs::File::open("alignments.sam")?;
/// let byte_reader = Arc::new(Mutex::new(SamByteReader::new(file, 1 << 24)?));
/// let references = byte_reader.lock().unwrap().header().reference_sequences().len();
///
/// let counts = read_sam_state_par_shared(&byte_reader, 1 << 24, 8, Count::default(),
///     |record: &RefSamRecord, state: &mut Count| state.0 += record.is_primary() as usize);
/// # Ok(())
/// # }
/// ```
pub fn read_sam_state_par_shared<G, T, State>(
    byte_reader: &Arc<Mutex<SamByteReader<T>>>,
    buffer_size: usize,
    num_threads: u32,
    mut global_state: State,
    f: G,
) -> State
where
    G: FnMut(&RefSamRecord, &mut State) -> () + Clone + Send,
    T: std::io::Read + std::marker::Send,
    State: Default + Clone + Send + Merge,
{
    std::thread::scope(|scope| {
        let mut threads = Vec::new();

        for _thread in 0..num_threads {
            let mut sam_reader = SamReader::with_capacity(buffer_size);
            let mut reader_local = byte_reader.clone();
            let mut state = global_state.clone();
            let mut f_local = f.clone();

            threads.push(scope.spawn(move || {
                while let Some(()) = sam_reader
                    .load_batch_par(&mut reader_local)
                    .expect("Batch is invalid")
                {
                    while let Some(record) = sam_reader.next() {
                        f_local(&record, &mut state);
                    }
                }
                state
            }));
        }

        for thread_guard in threads {
            let mut state = match thread_guard.join() {
                Ok(state) => state,
                Err(err) => panic!("Error {:?}", err),
            };
            global_state.merge_from(&mut state);
        }
    });

    global_state
}

/// Read-id variant of [`read_sam_state_par`]. Identical driver, except each record is passed to
/// the closure alongside its deterministic [`SamReadId`] `(ordinal, byte_offset, ..)`. See
/// [`SamReadId`] for the determinism guarantees.
///
/// This is additive: the plain variant is untouched, so existing callers pay nothing for the id
/// machinery. Batch decomposition, threading and merge semantics are exactly the same.
pub fn read_sam_state_par_ids<G, T, State>(
    file: T,
    buffer_size: usize,
    num_threads: u32,
    mut global_state: State,
    f: G,
) -> Result<State, std::io::Error>
where
    G: FnMut(&RefSamRecord, SamReadId, &mut State) -> () + Clone + Send,
    T: std::io::Read + std::marker::Send,
    State: Default + Clone + Send + Merge,
{
    let byte_reader = Arc::new(Mutex::new(SamByteReader::new(file, buffer_size)?));

    std::thread::scope(|scope| {
        let mut threads = Vec::new();

        for _thread in 0..num_threads {
            let mut sam_reader = SamReader::with_capacity(buffer_size);
            let reader_local = byte_reader.clone();
            let mut state = global_state.clone();
            let mut f_local = f.clone();

            threads.push(scope.spawn(move || {
                while let Some(()) = sam_reader
                    .load_batch_par_ids(&reader_local)
                    .expect("Batch is invalid")
                {
                    while let Some((record, read_id)) = sam_reader.next_id() {
                        f_local(&record, read_id, &mut state);
                    }
                }
                state
            }));
        }

        for thread_guard in threads {
            let mut state = match thread_guard.join() {
                Ok(state) => state,
                Err(err) => panic!("Error {:?}", err),
            };
            global_state.merge_from(&mut state);
        }
    });

    Ok(global_state)
}
