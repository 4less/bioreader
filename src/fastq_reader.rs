use memchr::memchr;
use std::{io::{Error, Read}, sync::{Arc, Mutex}};

// use memmap2::Mmap;

use crate::{byte_reader::{FillBufferPair, PairedByteReader}, sequence::fastq_record::{BufferPosition, RefRecord}};

#[cfg(windows)]
const LINE_ENDING: &'static str = "\r\n";
#[cfg(not(windows))]
const LINE_ENDING: &'static str = "\n";


pub struct PairedFastqReader<T> where T: Read{
    reader: Arc<Mutex<PairedByteReader<T>>>,
    pub buffer1: Vec<u8>,
    pub buffer1_fill: usize,
    pub buffer2: Vec<u8>,
    pub buffer2_fill: usize,
    buf1_pos: BufferPosition,
    buf2_pos: BufferPosition,
}

impl<T: Read> PairedFastqReader<T> {
    pub fn new(reader: Arc<Mutex<PairedByteReader<T>>>, capacity: usize) -> Self {
        PairedFastqReader {
            reader: reader,
            buffer1: vec![0; capacity],
            buffer1_fill: 0,
            buffer2: vec![0; capacity],
            buffer2_fill: 0,
            buf1_pos: BufferPosition::default(),
            buf2_pos: BufferPosition::default(),
        }
    }

    #[inline]
    pub fn load_batch_par(&mut self) -> Result<Option<()>, std::io::Error> {
        
        match self.reader.lock().expect("Locking ByteReader was unsuccessful")
                .fill_buf(&mut self.buffer1, &mut self.buffer2)? {
            
            Some((pos1, pos2)) => {
                self.buffer1_fill = pos1;
                self.buffer2_fill = pos2;
                // println!("{} {}", self.buffer1_fill, self.buffer2_fill);

                Ok(Some(()))
            },
            None => Ok(None)
        }
    }

    pub fn find_position(buffer: &mut [u8], buffer_pos: &mut BufferPosition) -> Option<()> {
        if buffer[buffer_pos.pos.1] != b'@' {
            println!("{}", std::str::from_utf8(&buffer[buffer_pos.pos.0..buffer_pos.pos.1+10]).unwrap());
        }
        assert!(buffer[buffer_pos.pos.1] == b'@');


        let pos1: usize = memchr(b'\n', &buffer[buffer_pos.pos.1..]).expect("1: Couldn't find newline") + buffer_pos.pos.1+1;
        let pos2 = memchr(b'\n', &buffer[pos1..]).expect("2: Couldn't find newline ") + pos1+1;
        
        assert!(buffer[pos2] == b'+');
        let pos3  = memchr(b'\n', &buffer[pos2..]).expect("3: Couldn't find newline") + pos2+1;
        let mut pos4 = memchr(b'\n', &buffer[pos3..]).unwrap_or(buffer.len() - pos3);
        pos4 += pos3;

        buffer_pos.pos = (buffer_pos.pos.1, pos4);
        buffer_pos.seq = pos1;
        buffer_pos.sep = pos2;
        buffer_pos.qual = pos3;

        Some(())
    }

    pub fn next_position(&mut self) -> Option<(RefRecord, RefRecord)> {
        self.buf1_pos.pos.1 += (self.buf1_pos.pos.1 > 0) as usize;
        self.buf2_pos.pos.1 += (self.buf2_pos.pos.1 > 0) as usize;

        let at_end1 = self.buf1_pos.pos.1 >= self.buffer1_fill;
        let at_end2 = self.buf2_pos.pos.1 >= self.buffer2_fill;

        assert_eq!(at_end1, at_end2);
        
        if at_end2 && at_end2 {
            let load = self.load_batch_par().expect("Valid filestream");

            self.buf1_pos.reset(0);
            self.buf2_pos.reset(0);
            // println!("Fill {} {}", self.buffer1_fill, self.buffer2_fill);
            if load.is_none() { return None; }
        }


        if self.buffer1[self.buf1_pos.pos.0] != b'@' {
            
            println!("{} {:?}", self.buffer1_fill, self.buf1_pos);
            println!(">>>>>>{}<<<<<<", std::str::from_utf8(&self.buffer1[0..100]).unwrap());
            println!("{}", self.buffer1.len());
        }

        assert!(self.buffer1[self.buf1_pos.pos.0] == b'@');
        assert!(self.buffer2[self.buf2_pos.pos.0] == b'@');
        assert!(self.buf1_pos.pos.0 == 0  || self.buffer1[self.buf1_pos.pos.0-1] == b'\n');
        assert!(self.buf2_pos.pos.0 == 0  || self.buffer2[self.buf2_pos.pos.0-1] == b'\n');

        PairedFastqReader::<T>::find_position(&mut self.buffer1[..self.buffer1_fill], &mut self.buf1_pos);
        PairedFastqReader::<T>::find_position(&mut self.buffer2[..self.buffer2_fill], &mut self.buf2_pos);

        // println!("Record length1: {}", (self.buf1_pos.pos.1 - self.buf1_pos.pos.0));
        // println!("Record length2: {}", (self.buf2_pos.pos.1 - self.buf2_pos.pos.0));

        let r1 = RefRecord {
            buffer: &self.buffer1,
            buf_pos: &self.buf1_pos,
        };
        let r2 = RefRecord {
            buffer: &self.buffer1,
            buf_pos: &self.buf1_pos,
        };

        Some((r1, r2))
    }


}


#[derive(Debug, Clone)]
pub struct FastqReader {
    pub buffer: Vec<u8>,
    pub buffer_size: usize,
    buf_pos: BufferPosition,
    finished: bool,
}




impl<'a> FastqReader { 

    pub fn new() -> FastqReader {
        FastqReader::with_capacity(usize::pow(2, 20))
    }

    pub fn with_capacity(capacity: usize) -> FastqReader {
        FastqReader {
            buffer: vec![0; capacity],
            buffer_size: 0,
            buf_pos: BufferPosition::default(),
            finished: false,
        }
    }

    // #[inline]
    // pub fn load_batch(&mut self, br: &mut impl FillBuffer) -> Result<Option<()>, std::io::Error> {
    //     self.buf_pos.reset(0);
    //     match br.fill_buf(&mut self.buffer)? {
    //         Some(bytes) => {
    //             self.buffer_size = bytes;
    //             Ok(Some(()))
    //         },
    //         None => Ok(None)
    //     }
    // }

    #[inline]
    pub fn load_batch<T>(&mut self, br: &mut T) -> Result<Option<()>, std::io::Error> where 
            T: std::io::Read {

        self.buf_pos.reset(0);
        match br.read(&mut self.buffer)? {
            bytes if bytes > 0 => {
                self.buffer_size = bytes;
                Ok(Some(()))
            },
            _ => Ok(None)
        }
    }

    #[inline]
    pub fn load_batch_par<T>(&mut self, br: &mut Arc<Mutex<T>>) -> Result<Option<()>, std::io::Error> where 
            T: std::io::Read {
        self.buf_pos.reset(0);
        match br.lock().expect("Locking ByteReader was unsuccessful").read(&mut self.buffer)? {
            bytes if bytes > 0 => {
                self.buffer_size = bytes;
                Ok(Some(()))
            },
            _ => Ok(None)
        }
    }

    pub fn next_position(&mut self) -> Option<RefRecord> {
        self.buf_pos.pos.1 += (self.buf_pos.pos.1 > 0) as usize;

        if self.buf_pos.pos.1 >= self.buffer_size {
            return None;
        }

        // println!("Size: {}\nBuffer: {}", self.buffer_size,  std::str::from_utf8(&self.buffer).unwrap());
        // println!("Size: {}", self.buffer_size);

        assert!(self.buf_pos.pos.0 == 0  || self.buffer[self.buf_pos.pos.0-1] == b'\n');
        assert!(self.buffer[self.buf_pos.pos.0] == b'@');
        let pos1: usize = memchr(b'\n', &self.buffer[self.buf_pos.pos.1..self.buffer_size]).expect("Couldn't find newline") + self.buf_pos.pos.1+1;
        let pos2 = memchr(b'\n', &self.buffer[pos1..self.buffer_size]).expect("Couldn't find newline") + pos1+1;
        
        assert!(self.buffer[pos2] == b'+');
        let pos3  = memchr(b'\n', &self.buffer[pos2..self.buffer_size]).expect("Couldn't find newline") + pos2+1;
        let mut pos4 = memchr(b'\n', &self.buffer[pos3..self.buffer_size]).unwrap_or(self.buffer_size - pos3);
        pos4 += pos3;

        self.buf_pos.pos = (self.buf_pos.pos.1, pos4);
        self.buf_pos.seq = pos1;
        self.buf_pos.sep = pos2;
        self.buf_pos.qual = pos3;

        Some(RefRecord {
            buffer: &self.buffer,
            buf_pos: &self.buf_pos,
        })
    }


    // / Searches the next FASTQ record and returns a [RefRecord](struct.RefRecord.html) that
    // / borrows its data from the underlying buffer of this reader.
    // /
    // / # Example:
    // /
    // / ```no_run
    // / use seq_io::fastq::{Reader, Record};
    // /
    // / let mut reader = Reader::from_path("seqs.fastq").unwrap();
    // /
    // / while let Some(record) = reader.next() {
    // /     let record = record.unwrap();
    // /     println!("{}", record.id().unwrap());
    // / }
    // / ```
    #[allow(clippy::should_implement_trait)]
    #[inline]
    pub fn next(&mut self) -> Option<Result<RefRecord, Error>> {
        if self.finished {
            return None;
        }

        match self.next_position() {
            Some(_rec) => {},
            None => { return None }
        };

        Some(Ok(RefRecord {
            buffer: &self.buffer,
            buf_pos: &self.buf_pos,
        }))
    }
}


// // Implement `Iterator` for `Fibonacci`.
// // The `Iterator` trait only requires a method to be defined for the `next` element.
// impl<'a> Iterator for FastqReader<'a> {
//     // We can refer to this type using Self::Item
//     type Item = RefRecord<'a>;

//     fn next(&mut self) -> Option<Self::Item> {
//         self.next().
//         Some()
//     }
// }
