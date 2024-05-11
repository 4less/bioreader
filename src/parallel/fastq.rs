use std::sync::{mpsc, Arc, Mutex};

use crate::{
    byte_reader::{ByteReader, PairedByteReader},
    fastq_reader::{FastqReader, PairedFastqReader},
    sequence::fastq_record::RefRecord,
};

pub fn read_fastq_par<G, T>(
    file: T,
    buffer_size: usize,
    num_threads: u32,
    f: G,
) -> Result<usize, std::io::Error>
where
    G: Fn(&RefRecord) -> usize + Clone + Send,
    T: std::io::Read + std::marker::Send,
{

    let byte_reader = Arc::new(Mutex::new(ByteReader::new(file, buffer_size)?));
    let mut total = 0;
    std::thread::scope(|scope| { 


        let mut threads = Vec::new();

        let (tx, rx) = mpsc::channel::<usize>();

        // Results thread
        scope.spawn(move || {
            let mut total_result = 0;
            for e in rx {
                total_result += e;
            }
            println!("Totale resultate: {}", total_result);
        });

        for _thread in 0..num_threads {
            let mut fastq_reader = FastqReader::with_capacity(buffer_size);
            let mut reader_clone = byte_reader.clone(); // Need to clone for move
            let mut count = 0;
            let mut result_buffer = 0;
            let f_clone = f.clone();
            let tx: mpsc::Sender<usize> = tx.clone();

            threads.push(scope.spawn(move || {
                while let Some(()) = fastq_reader
                    .load_batch_par(&mut reader_clone)
                    .expect("Batch is invalid")
                {
                    while let Some(record) = fastq_reader.next() {
                        count += 1;
                        if !record.valid() {
                            panic!("Invalid record {}", record)
                        }
                        result_buffer += f_clone(&record);

                        if count > 10_000 {
                            let _res = tx.send(result_buffer);
                            result_buffer = 0;
                            count -= 10_000;
                        }
                    }
                }
                count
            }));
        }

        for thread_guard in threads {
            total += thread_guard.join().unwrap();
        }

    
    });
    Ok(total)
}

pub fn read_fastq_pair_par<G, T>(
    file1: T,
    file2: T,
    buffer_size: usize,
    num_threads: u32,
    f: G,
) -> Result<u64, std::io::Error>
where
    G: Fn(&RefRecord, &RefRecord) -> usize + Clone + Send,
    T: std::io::Read + std::marker::Send,
{
    let mut total = 0;

    // This scope guarantees that all threads finish within the scope. This way, lifetimes of G and T do not need to be 'static
    std::thread::scope(|scope| { 
        let byte_reader = Arc::new(Mutex::new(PairedByteReader::new(file1, file2, buffer_size)));
        let mut threads = Vec::new();

        let (tx, rx) = mpsc::channel::<usize>();

        // Results thread. Clarify what this wants to do
        scope.spawn(move || {
            let mut total_result = 0;
            for e in rx {
                total_result += e;
            }
            println!("Totale result: {}", total_result);
        });

        for _thread in 0..num_threads {
            let mut fastq_reader = PairedFastqReader::new(byte_reader.clone(), buffer_size);
            let mut count = 0;
            let mut count2 = 0;
            let mut result_buffer = 0;
            let f_clone = f.clone();
            let tx: mpsc::Sender<usize> = tx.clone();
            
            threads.push(scope.spawn(move || {
                while let Some((record1, record2)) = fastq_reader.next() {
                    count += 1;
                    if !record1.valid() {
                        panic!("Invalid record {}", record1)
                    }
                    if !record2.valid() {
                        panic!("Invalid record {}", record2)
                    }

                    result_buffer += f_clone(&record1, &record2);
                    count2 += 1;

                    if count > 10_000 {
                        let _res = tx.send(result_buffer);
                        result_buffer = 0;
                        count -= 10_000;
                    }
                }
                let _ = tx.send(result_buffer);
                count2
            }));
        }
        for thread_guard in threads {
            total += thread_guard.join().unwrap();
        }
    });
    Ok(total)
}
