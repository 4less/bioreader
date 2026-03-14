use std::sync::{mpsc, Arc, Mutex};

use crate::{
    fasta_byte_reader::FastaByteReader, fasta_reader::FastaReader, fastq_byte_reader::{FastqByteReader, FastqPairedByteReader}, fastq_reader::{self, FastqReader, PairedFastqReader}, sequence::{fasta_record::OwnedFastaRecord, fastq_record::RefFastqRecord}
};

pub fn read_fastq_par<G, T>(
    file: T,
    buffer_size: usize,
    num_threads: u32,
    f: G,
) -> Result<usize, std::io::Error>
where
    G: Fn(&RefFastqRecord) -> usize + Clone + Send,
    T: std::io::Read + std::marker::Send,
{

    let byte_reader = Arc::new(Mutex::new(FastqByteReader::new(file, buffer_size)?));
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

pub fn read_fastq_pair_par<G, T, O>(
    file1: T,
    file2: T,
    buffer_size: usize,
    num_threads: u32,
    f: G,
) -> Result<u64, std::io::Error>
where
    G: FnMut(&RefFastqRecord, &RefFastqRecord) -> O + Clone + Send,
    T: std::io::Read + std::marker::Send,
{
    let mut total = 0;

    // This scope guarantees that all threads finish within the scope. This way, lifetimes of G and T do not need to be 'static
    std::thread::scope(|scope| { 
        let byte_reader = Arc::new(Mutex::new(FastqPairedByteReader::new(file1, file2, buffer_size)));
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
            let mut f_clone = f.clone();
            let tx: mpsc::Sender<usize> = tx.clone();
            
            threads.push(scope.spawn(move || {
                while let Some((record1, record2)) = fastq_reader.next() {
                    count += 1;
                    if !record1.valid_extended() {
                        panic!("Invalid record {}", record1)
                    }
                    if !record2.valid_extended() {
                        panic!("Invalid record {}", record2)
                    }

                    f_clone(&record1, &record2);
                    count2 += 1;

                    if count > 10_000 {
                        let _res = tx.send(result_buffer);
                        result_buffer = 0;
                        // count -= 10_000;
                    }
                }
                println!("Read count: {}", count);

                let _ = tx.send(result_buffer);
                // count2
                (count2, f_clone)
            }));
        }
        
        for thread_guard in threads {
            let (count, _closure) = thread_guard.join().unwrap();
            total += count;
        }
    });
    Ok(total)
}


pub fn read_fastq_single_end_state_par<G, T, State>(
    file: T,
    buffer_size: usize,
    num_threads: u32,
    f: G,
) -> Result<State, std::io::Error>
where
    G: FnMut(&RefFastqRecord, &mut State) -> () + Clone + Send,
    T: std::io::Read + std::marker::Send,
    State: Default + Send + Merge,
{

    let byte_reader = Arc::new(Mutex::new(FastqByteReader::new(file, buffer_size)?));
    let mut total = 0;

    let mut global_state = State::default();

    std::thread::scope(|scope| { 


        let mut threads = Vec::new();

        // let (tx, rx) = mpsc::channel::<usize>();

        // // Results thread
        // scope.spawn(move || {
        //     let mut total_result = 0;
        //     for e in rx {
        //         total_result += e;
        //     }
        //     println!("Totale resultate: {}", total_result);
        // });

        for _thread in 0..num_threads {
            let mut fastq_reader = FastqReader::with_capacity(buffer_size);
            let mut reader_local = byte_reader.clone(); // Need to clone for move
            let mut count = 0;
            let mut state = State::default();
            let mut result_buffer = 0;
            let mut f_local = f.clone();
            // let tx: mpsc::Sender<usize> = tx.clone();

            threads.push(scope.spawn(move || {
                while let Some(()) = fastq_reader
                    .load_batch_par(&mut reader_local)
                    .expect("Batch is invalid")
                {
                    while let Some(record) = fastq_reader.next() {
                        count += 1;
                        if !record.valid_extended() {
                            panic!("Invalid record {}", record)
                        }
                        
                        f_local(&record, &mut state);

                        // if count > 10_000 {
                        //     let _res = tx.send(result_buffer);
                        //     result_buffer = 0;
                        //     count -= 10_000;
                        // }
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


pub fn read_fastq_paired_end_state_par<G, T, State>(
    file1: T,
    file2: T,
    buffer_size: usize,
    num_threads: u32,
    mut global_state: State,
    f: G,
) -> Result<State, std::io::Error>
where
    G: FnMut(&RefFastqRecord, &RefFastqRecord, &mut State) -> () + Clone + Send,
    T: std::io::Read + std::marker::Send,
    State: Default + Clone + Send + Merge,
{
    let mut total = 0;

    // let mut global_state = State::default();



    // This scope guarantees that all threads finish within the scope. This way, lifetimes of G and T do not need to be 'static
    std::thread::scope(|scope| { 
        let byte_reader = Arc::new(Mutex::new(FastqPairedByteReader::new(file1, file2, buffer_size)));
        let mut threads = Vec::new();



        for _thread in 0..num_threads {
            let mut fastq_reader = PairedFastqReader::new(byte_reader.clone(), buffer_size);
            let result_buffer = 0;
            let mut state = global_state.clone();
            let mut f_local = f.clone();
            // let tx: mpsc::Sender<usize> = tx.clone();
            
            threads.push(scope.spawn(move || {
                while let Some((record1, record2)) = fastq_reader.next() {
                    if !record1.valid_extended() {
                        panic!("Invalid record {}", record1)
                    }
                    if !record2.valid_extended() {
                        panic!("Invalid record {}", record2)
                    }

                    f_local(&record1, &record2, &mut state);
                }

                // let _ = tx.send(result_buffer);

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

// pub fn read_fastq_paired_end_state_par<G, T, State>(
//     file1: T,
//     file2: T,
//     buffer_size: usize,
//     num_threads: u32,
//     f: G,
// ) -> Result<State, std::io::Error>
// where
//     G: FnMut(&RefFastqRecord, &RefFastqRecord, &mut State) -> () + Clone + Send,
//     T: std::io::Read + std::marker::Send,
//     State: Default + Send + Merge,
// {
//     let mut total = 0;

//     let mut global_state = State::default();

//     // This scope guarantees that all threads finish within the scope. This way, lifetimes of G and T do not need to be 'static
//     std::thread::scope(|scope| { 
//         let byte_reader = Arc::new(Mutex::new(FastqPairedByteReader::new(file1, file2, buffer_size)));
//         let mut threads = Vec::new();



//         for _thread in 0..num_threads {
//             let mut fastq_reader = PairedFastqReader::new(byte_reader.clone(), buffer_size);
//             let result_buffer = 0;
//             let mut state = State::default();
//             let mut f_local = f.clone();
//             // let tx: mpsc::Sender<usize> = tx.clone();
            
//             threads.push(scope.spawn(move || {
//                 while let Some((record1, record2)) = fastq_reader.next() {
//                     if !record1.valid_extended() {
//                         panic!("Invalid record {}", record1)
//                     }
//                     if !record2.valid_extended() {
//                         panic!("Invalid record {}", record2)
//                     }

//                     f_local(&record1, &record2, &mut state);
//                 }

//                 // let _ = tx.send(result_buffer);

//                 state
//             }));
//         }
        
//         for thread_guard in threads {
//             let mut state = match thread_guard.join() {
//                 Ok(state) => state,
//                 Err(err) => panic!("Error {:?}", err),
//             };
//             global_state.merge_from(&mut state);
//         }
//     });
//     Ok(global_state)
// }




pub fn read_fasta_par<G, T>(
    file: T,
    buffer_size: usize,
    num_threads: u32,
    f: G,
) -> Result<usize, std::io::Error>
where
    G: Fn(&OwnedFastaRecord) -> usize + Clone + Send,
    T: std::io::Read + std::marker::Send,
{

    let byte_reader = Arc::new(Mutex::new(FastaByteReader::new(file, buffer_size)?));
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
            let mut fasta_reader = FastaReader::with_capacity(buffer_size);
            let mut reader_clone = byte_reader.clone(); // Need to clone for move
            let mut count = 0;
            let mut result_buffer = 0;
            let f_clone = f.clone();
            let tx: mpsc::Sender<usize> = tx.clone();

            let mut record = OwnedFastaRecord::new();
            threads.push(scope.spawn(move || {

                while let Some(()) = fasta_reader
                    .load_batch_par(&mut reader_clone)
                    .expect("Batch is invalid")
                {
                    while let Some(_) = fasta_reader.next(&mut record) {
                        count += 1;
                        if !record.valid_extended() {
                            panic!("Invalid record {}", record.to_string())
                        }
                        result_buffer += f_clone(&record);

                        if count > 10_000 {
                            let _res = tx.send(result_buffer);
                            result_buffer = 0;
                        }
                    }
                }
                let _res = tx.send(result_buffer);
                count
            }));
        }

        for thread_guard in threads {
            let result = thread_guard.join().unwrap();
            total += result;
        }

    
    });
    Ok(total)
}


pub trait Merge {
    fn merge_from(self: &mut Self, other: &mut Self);
}


pub fn read_fasta_par2<G, H, T, B>(
    file: T,
    buffer_size: usize,
    num_threads: u32,
    f: G,
    h: H,
) -> Result<usize, std::io::Error>
where
    G: Fn(&OwnedFastaRecord, &mut B) -> usize + Clone + Send,
    H: Fn(B) -> Option<()> + Clone + Send,
    T: std::io::Read + std::marker::Send,
    B: Clone + Send + Default + Merge,
{

    let byte_reader = Arc::new(Mutex::new(FastaByteReader::new(file, buffer_size)?));
    
    let mut total = 0;
    std::thread::scope(|scope| { 
        let mut threads = Vec::new();

        let (tx, rx) = mpsc::channel::<B>();

        // Results thread
        scope.spawn(move || {
            let mut total_result = 0;
            for e in rx {
                h(e);
            }
            println!("Totale resultate: {}", total_result);
        });

        for _thread in 0..num_threads {
            let mut fasta_reader = FastaReader::with_capacity(buffer_size);
            let mut reader_clone = byte_reader.clone(); // Need to clone for move
            let mut count = 0;
            
            let mut buffer: B = B::default();

            let f_clone = f.clone();
            let tx: mpsc::Sender<B> = tx.clone();

            let mut record = OwnedFastaRecord::new();
            threads.push(scope.spawn(move || {

                while let Some(()) = fasta_reader
                    .load_batch_par(&mut reader_clone)
                    .expect("Batch is invalid")
                {
                    while let Some(_) = fasta_reader.next(&mut record) {
                        count += 1;
                        if !record.valid_extended() {
                            panic!("Invalid record {}", record.to_string())
                        }
                        f_clone(&record, &mut buffer);

                        if count > 10_000 {
                            let _res = tx.send(buffer);
                            buffer = B::default();
                        }
                    }
                }
                let _res = tx.send(buffer);
                count
            }));
        }

        for thread_guard in threads {
            let result = thread_guard.join().unwrap();
            total += result;
        }

    
    });
    Ok(total)
}
