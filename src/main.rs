#![feature(slice_internals)]
use core::slice::memchr::memchr;
use std::{
    fs::File,
    io::{Error, Read},
    path::Path,
    sync::{Arc, Mutex},
};

use bioreader::{fasta_byte_reader::FastaByteReader, parallel::fastq::read_fasta_par, utils};
use bioreader::{
    fastq_byte_reader::FastqByteReader,
    fasta_reader::FastaReader,
    fastq_reader::FastqReader,
    parallel::fastq::read_fastq_pair_par,
    sequence::{fasta_record::OwnedFastaRecord, fastq_record::RefFastqRecord},
};
use flate2::read::GzDecoder;
use memmap2::Mmap;

fn read_fastq<T>(file: T, buffer_size: usize) -> Result<usize, std::io::Error>
where
    T: Read,
{
    let mut byte_reader = FastqByteReader::new(file, buffer_size)?;
    // let mut byte_reader = BufReader::with_capacity(buffer_size, file);
    let mut fastq_reader = FastqReader::with_capacity(buffer_size);

    let mut count = 0;
    let mut total_length = 0;
    while let Some(_ok) = fastq_reader.load_batch(&mut byte_reader)? {
        while let Some(record) = fastq_reader.next() {
            // println!("{}", record);
            count += 1;
            total_length += record.seq().len();
        }
    }
    println!("Count: {count}");

    Ok(total_length)
}

fn read_fastq_mm(path: impl AsRef<Path>) -> Result<usize, std::io::Error> {
    let file: File = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.as_ref().display(), why),
        Ok(file) => file,
    };
    // Create the memory mapped buffer
    let mmap = unsafe { Mmap::map(&file).expect("Error mapping file") };

    let mut count = 0;
    let mut pos = 0;
    while let Some(res) = memchr(b'\n', &mmap[pos..]) {
        pos += res + 1;
        count += 1;
    }

    Ok(count)
}

fn test_fa() -> Result<(), Error> {
    let path: &Path = Path::new("data/fasta/test.fna");

    let buffer_size = 1_000_000;
    let file: File = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why),
        Ok(file) => file,
    };

    let mut byte_reader = Arc::new(Mutex::new(FastaByteReader::new(file, buffer_size)?));
    let mut fasta_reader = FastaReader::with_capacity(buffer_size);

    let mut count = 0;
    let mut total_length = 0;

    let mut record = OwnedFastaRecord::new();

    while let Some(()) = fasta_reader
        .load_batch_par(&mut byte_reader)
        .expect("Batch is invalid")
    {
        while let Some(_) = fasta_reader.next(&mut record) {
            count += 1;
            total_length += record.seq().len();
            // println!("{}", record.to_string());
        }
        // println!("{}", String::from_utf8_lossy(record.head()));
    }
    

    println!("Count: {count}");

    Ok(())
}

fn test_paired_fq_parallel() -> Result<(), Error> {
    let worker = |rec1: &RefFastqRecord, rec2: &RefFastqRecord| {
        if !rec1.valid() || !rec2.valid() {
            println!("problemo");
        }
        let mut kmer_count = 0;

        kmer_count += rec1.perfect() as u32;

        kmer_count as usize
    };

    // let path_1: &Path = Path::new("data/fastq/small_test_gzip_1.fq.gz");
    // let path_2: &Path = Path::new("data/fastq/small_test_gzip_2.fq.gz");

    let path_1: &Path = Path::new("data/large_data/fastq/weird_1.fq.gz");
    let path_2: &Path = Path::new("data/large_data/fastq/weird_2.fq.gz");

    let file_1 = File::open(&path_1)?;
    let file_2 = File::open(&path_2)?;

    let (duration, result) = utils::time(move || {
        read_fastq_pair_par(
            GzDecoder::new(file_1),
            GzDecoder::new(file_2),
            usize::pow(2, 24),
            8,
            worker,
        )
    })
    .expect("No Shit");
    println!(
        "Time elapsed in read_fastq_mm() is: {:?} -> Result: {}",
        duration, result
    );

    Ok(())
}

fn test_fa_parallel() -> Result<(), Error> {
    let worker = |rec: &OwnedFastaRecord| {
        if !rec.valid() {
            println!("problemo");
        }
        let mut kmer_count = 0;

        kmer_count += rec.perfect() as u32;

        kmer_count as usize
    };

    let path: &Path = Path::new("/usr/users/QIB_fr017/fritsche/ProjectsPrivate/bioreader/data/fasta/test.fna");
    let path: &Path = Path::new("/usr/users/QIB_fr017/fritsche/ProjectsPrivate/flexalign/data/large_data/fasta/combined.fna");

    let file = File::open(&path)?;

    let (duration, result) = utils::time(move || {
        read_fasta_par(
            file,
            usize::pow(2, 24),
            4,
            worker,
        )
    })
    .expect("No Shit");
    println!(
        "Time elapsed in read_fasta_par() is: {:?} -> Result: {}",
        duration, result
    );

    Ok(())
}

fn main() {
    let _result = test_fa();
    // let _result = test_paired_fq_parallel();
    let _result = test_fa_parallel();

    // exit(9);

    // let path_fq: &Path = Path::new("data/fastq/SRR7280802_1c.fastq");
    // let path_fq_gz: &Path = Path::new("data/fastq/SRR7280802_1c.fq.gz");

    // println!("Is gz? {}", is_gzip(path_fq_gz).unwrap());
    // println!("Is gz? {}", is_gzip(path_fq).unwrap());

    // let mut file1 = match File::open(&path_fq) {
    //     Err(why) => panic!("couldn't open {}: {}", &path_fq.display(), why),
    //     Ok(file) => file,
    // };
    // // let mut file1_gz = match File::open(&path_fq_gz).map(|file| GzDecoder::new(file)) {
    // //     Err(why) => panic!("couldn't open {}: {}", path_fq.display(), why),
    // //     Ok(file_gz) => file_gz,
    // // };

    // let mut file1_gz = match File::open(&path_fq_gz) {
    //     Err(why) => panic!("couldn't open {}: {}", path_fq.display(), why),
    //     Ok(file) => file,
    // };

    // // let file1b: Box<dyn Read> = Box::new(GzDecoder::new(&file1));
    // // let filegz: Box<dyn Read> = Box::new(GzDecoder::new(&file1_gz));

    // // let newred = BufReader::new(file1b);
    // // let newred2 = GzDecoder::new(file1b);

    // let buffer_size = 1_000_000;
    // let (duration, result) =
    //     utils::time(|| read_fastq(&file1, buffer_size)).expect("No Shit");
    // println!("Time elapsed in read_fastq() is: {:?} -> Result: {}", duration, result);

    // let (duration, result) =
    //     utils::time(|| read_fastq(GzDecoder::new(&file1_gz), buffer_size)).expect("Cannot create GZDecoder");
    // println!("Time elapsed in read_fastq() is: {:?} -> Result: {}", duration, result);

    // let (duration, result) =
    //     utils::time(move || read_fastq_mm(path_fq)).expect("No Shit");
    // println!("Time elapsed in read_fastq_mm() is: {:?} -> Result: {}", duration, result);

    println!("Hello, world!");
}
