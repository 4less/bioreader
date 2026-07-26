use memchr::memchr;
use std::{io::{Error, Read}, sync::{Arc, Mutex}};

// use memmap2::Mmap;

use crate::{fastq_byte_reader::{BatchBase, FillBufferPair, FastqPairedByteReader, ReadId}, sequence::fastq_record::{BufferPosition, RefFastqRecord}};

#[cfg(windows)]
const LINE_ENDING: &'static str = "\r\n";
#[cfg(not(windows))]
const LINE_ENDING: &'static str = "\n";


pub struct PairedFastqReader<T> where T: Read{
    reader: Arc<Mutex<FastqPairedByteReader<T>>>,
    pub buffer1: Vec<u8>,
    pub buffer1_fill: usize,
    pub buffer2: Vec<u8>,
    pub buffer2_fill: usize,
    buf1_pos: BufferPosition,
    buf2_pos: BufferPosition,
    // Read-id state for the `*_ids` path only; the plain `next()` never reads these.
    // `batch_base` is the id of the current batch's first record (fetched under the reader
    // lock by `load_batch_par_ids`); `idx_in_batch` counts records emitted from this batch.
    batch_base: BatchBase,
    idx_in_batch: u64,
}

impl<T: Read> PairedFastqReader<T> {
    pub fn new(reader: Arc<Mutex<FastqPairedByteReader<T>>>, capacity: usize) -> Self {
        Self {
            reader: reader,
            buffer1: vec![0; capacity],
            buffer1_fill: 0,
            buffer2: vec![0; capacity],
            buffer2_fill: 0,
            buf1_pos: BufferPosition::default(),
            buf2_pos: BufferPosition::default(),
            batch_base: BatchBase::default(),
            idx_in_batch: 0,
        }
    }

    #[inline]
    pub fn load_batch_par(&mut self) -> Result<Option<()>, std::io::Error> {
        
        match self.reader.lock().expect("Locking ByteReader was unsuccessful")
                .fill_buf(&mut self.buffer1, &mut self.buffer2)? {
            
            Some((pos1, pos2)) => {
                self.buffer1_fill = pos1;
                self.buffer2_fill = pos2;

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

    pub fn next(&mut self) -> Option<(RefFastqRecord, RefFastqRecord)> {
        self.buf1_pos.pos.1 += (self.buf1_pos.pos.1 > 0) as usize;
        self.buf2_pos.pos.1 += (self.buf2_pos.pos.1 > 0) as usize;

        let at_end1 = self.buf1_pos.pos.1 >= self.buffer1_fill;
        let at_end2 = self.buf2_pos.pos.1 >= self.buffer2_fill;

        assert_eq!(at_end1, at_end2);
        
        if at_end2 && at_end2 {
            let load = self.load_batch_par().expect("Valid filestream");

            self.buf1_pos.reset(0);
            self.buf2_pos.reset(0);
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

        let r1 = RefFastqRecord {
            buffer: &self.buffer1,
            buf_pos: &self.buf1_pos,
        };
        let r2 = RefFastqRecord {
            buffer: &self.buffer2,
            buf_pos: &self.buf2_pos,
        };

        Some((r1, r2))
    }

    /// Like [`load_batch_par`](Self::load_batch_par), but pulls the batch through
    /// [`fill_buf_ids`](FastqPairedByteReader::fill_buf_ids) and records the batch's id base,
    /// resetting the within-batch record counter. Used only by [`next_id`](Self::next_id).
    #[inline]
    pub fn load_batch_par_ids(&mut self) -> Result<Option<()>, std::io::Error> {
        match self.reader.lock().expect("Locking ByteReader was unsuccessful")
                .fill_buf_ids(&mut self.buffer1, &mut self.buffer2)? {

            Some((pos1, pos2, base)) => {
                self.buffer1_fill = pos1;
                self.buffer2_fill = pos2;
                self.batch_base = base;
                self.idx_in_batch = 0;

                Ok(Some(()))
            },
            None => Ok(None)
        }
    }

    /// The read-id counterpart of [`next`](Self::next): yields each pair together with its
    /// deterministic [`ReadId`]. Behaves identically to `next` for record iteration — this is
    /// deliberately a near-duplicate rather than a wrapper, because `next` returns records
    /// that borrow `self`, leaving no room to also hand back an id from a shared helper. Keep
    /// the two in sync; the only differences are the `_ids` batch load and the id it derives
    /// from the batch base plus each record's in-buffer offset.
    pub fn next_id(&mut self) -> Option<(RefFastqRecord, RefFastqRecord, ReadId)> {
        self.buf1_pos.pos.1 += (self.buf1_pos.pos.1 > 0) as usize;
        self.buf2_pos.pos.1 += (self.buf2_pos.pos.1 > 0) as usize;

        let at_end1 = self.buf1_pos.pos.1 >= self.buffer1_fill;
        let at_end2 = self.buf2_pos.pos.1 >= self.buffer2_fill;

        assert_eq!(at_end1, at_end2);

        if at_end1 && at_end2 {
            let load = self.load_batch_par_ids().expect("Valid filestream");

            self.buf1_pos.reset(0);
            self.buf2_pos.reset(0);
            if load.is_none() { return None; }
        }

        assert!(self.buffer1[self.buf1_pos.pos.0] == b'@');
        assert!(self.buffer2[self.buf2_pos.pos.0] == b'@');
        assert!(self.buf1_pos.pos.0 == 0  || self.buffer1[self.buf1_pos.pos.0-1] == b'\n');
        assert!(self.buf2_pos.pos.0 == 0  || self.buffer2[self.buf2_pos.pos.0-1] == b'\n');

        PairedFastqReader::<T>::find_position(&mut self.buffer1[..self.buffer1_fill], &mut self.buf1_pos);
        PairedFastqReader::<T>::find_position(&mut self.buffer2[..self.buffer2_fill], &mut self.buf2_pos);

        // `pos.0` is the record's start within the batch buffer, which begins at the batch's
        // file offset; adding the two gives the record's absolute byte offset in the file.
        let read_id = ReadId {
            ordinal: self.batch_base.ordinal + self.idx_in_batch,
            byte_offset1: self.batch_base.byte_offset1 + self.buf1_pos.pos.0 as u64,
            byte_offset2: self.batch_base.byte_offset2 + self.buf2_pos.pos.0 as u64,
            batch_id: self.batch_base.batch_id,
            idx_in_batch: self.idx_in_batch,
        };
        self.idx_in_batch += 1;

        let r1 = RefFastqRecord {
            buffer: &self.buffer1,
            buf_pos: &self.buf1_pos,
        };
        let r2 = RefFastqRecord {
            buffer: &self.buffer2,
            buf_pos: &self.buf2_pos,
        };

        Some((r1, r2, read_id))
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

    pub fn next(&mut self) -> Option<RefFastqRecord> {
        self.buf_pos.pos.1 += (self.buf_pos.pos.1 > 0) as usize;

        if self.buf_pos.pos.1 >= self.buffer_size {
            return None;
        }


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

        Some(RefFastqRecord {
            buffer: &self.buffer,
            buf_pos: &self.buf_pos,
        })
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
