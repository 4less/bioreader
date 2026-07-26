use memchr::{memchr_iter, memmem::Finder};
use memmap2::Mmap;
use std::{
    fs::File, io::{Error, Read}
};



pub trait FillBuffer {
    fn fill_buf(&mut self, buffer: &mut Vec<u8>) -> Result<Option<usize>, Error>;
}

pub trait FillBufferPair {
    fn fill_buf(
        &mut self,
        buffer1: &mut Vec<u8>,
        buffer2: &mut Vec<u8>,
    ) -> Result<Option<(usize, usize)>, Error>;
}

pub struct ByteReaderMmap {
    mmap: Mmap,
    buffer_fill: usize,
    position: usize,
    record_finder: Finder<'static>,
}

impl ByteReaderMmap {
    pub fn with_capacity(file: &File, buffer_size: usize) -> Result<ByteReaderMmap, Error> {
        let br = ByteReaderMmap {
            mmap: unsafe { Mmap::map(file).expect("Error mapping file") },
            buffer_fill: buffer_size,
            position: 0,
            record_finder: Finder::new("\n@"),
        };
        Ok(br)
    }
}

impl FillBuffer for ByteReaderMmap {
    fn fill_buf(&mut self, buffer: &mut Vec<u8>) -> Result<Option<usize>, Error> {
        if self.position >= self.mmap.len() {
            return Ok(None);
        }

        if self.position + self.buffer_fill >= self.mmap.len() {
            let chunk_size = self.mmap.len() - self.position;
            buffer[..chunk_size].copy_from_slice(&self.mmap[self.position..]);
            self.position = self.mmap.len();
            return Ok(Some(chunk_size));
        }
        //TODO: If one record does not fit within chunk, what to do?

        let chunk_size = match self
            .record_finder
            .find(&self.mmap[self.position + self.buffer_fill..])
        {
            Some(pos) => self.buffer_fill + pos + 1,
            None => self.mmap.len(),
        };

        if chunk_size > 2 * self.buffer_fill {
            panic!("Record does not fit within buffer size for ByteReaderMMap");
        }

        while chunk_size > buffer.len() {
            buffer.resize(buffer.len() * 2, 0);
        }

        buffer[..chunk_size].copy_from_slice(&self.mmap[self.position..self.position + chunk_size]);
        self.position += chunk_size;

        Ok(Some(chunk_size))
    }
}

// impl<T: std::io::Read> FillBuffer for ByteReader<T> {
//     fn fill_buf(&mut self, buffer: &mut Vec<u8>) -> Result<Option<usize>, Error> {
//         // Implies that file has been read to the end, but there still might be data in the buffer
//         if self.finished {
//             if self.buffer_fill == 0 {
//                 return Ok(None);
//             } else if self.buffer_fill <= buffer.len() {
//                 buffer[..self.buffer_fill].copy_from_slice(&self.buffer[..self.buffer_fill]);
//                 let fill = self.buffer_fill;
//                 self.buffer_fill = 0;
//                 return Ok(Some(fill));
//             } // File has been read and buffer is empty
//         }

//         // Find record in slice, at this point buffer must be filled to the brim
//         assert!(self.buffer_fill >= buffer.len());

//         let buffer_slice = &self.buffer[..buffer.len()];

//         // Find end of last complete Fastq record in local buffer
//         let mut index = 0;
//         for (i, &c) in buffer_slice.iter().enumerate().rev() {
//             if i > 1 && c == b'@' && buffer_slice[i - 1] == b'\n' && buffer_slice[i - 2] != b'+' {
//                 index = i;
//                 break;
//             }
//         }

//         // Copy local buffer of complete Fastq records into external buffer
//         buffer[..index].copy_from_slice(&buffer_slice[..index]);

//         self.buffer.drain(0..index);
//         self.buffer_fill -= index; // unprocessed bytes
//         self.buffer.resize(self.buffer.capacity(), 0); // restore buffer fill

//         if self.buffer_fill < self.buffer.len() {
//             self.read()?;
//         }

//         Ok(Some(index))
//     }
// }

pub struct FastqByteReader<T>
where
    T: std::io::Read,
{
    file: T,
    buffer_fill: usize,
    buffer: Vec<u8>,
    finished: bool,
}

impl<T: std::io::Read> std::io::Read for FastqByteReader<T> {
    
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        // Implies that file has been read to the end, but there still might be data in the buffer
        if self.finished {
            if self.buffer_fill == 0 {
                return Ok(0);
            } else if self.buffer_fill <= buf.len() {
                buf[..self.buffer_fill].copy_from_slice(&self.buffer[..self.buffer_fill]);
                let fill = self.buffer_fill;
                self.buffer_fill = 0;
                return Ok(fill);
            } // File has been read and buffer is empty
        }

        // Find record in slice, at this point buffer must be filled to the brim
        // assert!(self.buffer_fill >= buf.len());

        let buffer_slice = &self.buffer[..self.buffer_fill];

        // Find end of last complete Fastq record in local buffer
        let mut index = 0;
        for (i, &c) in buffer_slice.iter().enumerate().rev() {
            if i > 1 && c == b'@' && buffer_slice[i - 1] == b'\n' && buffer_slice[i - 2] != b'+' {
                index = i;
                break;
            }
        }

        // Copy local buffer of complete Fastq records into external buffer
        buf[..index].copy_from_slice(&buffer_slice[..index]);
        // println!("- Bytes: {index} {}", buf.len());
        // println!("- End: {}", std::str::from_utf8(&buf[index-200..index]).unwrap());

        self.buffer.drain(0..index);
        self.buffer_fill -= index; // unprocessed bytes
        self.buffer.resize(self.buffer.capacity(), 0); // restore buffer fill


        if self.buffer_fill < self.buffer.len() {
            self.read()?;
        }

        Ok(index)
    }
}



impl<T: std::io::Read> FastqByteReader<T> {
    pub fn new(reader: T, chunk_size: usize) -> Result<FastqByteReader<T>, Error> {
        let mut br = FastqByteReader {
            file: reader,
            buffer_fill: 0,
            buffer: vec![0; chunk_size],
            finished: false,
        };
        br.read()?;
        Ok(br)
    }

    pub fn read(&mut self) -> Result<Option<usize>, Error> {
        if self.finished {
            return Ok(None);
        }

        loop {
            let n_bytes = self.file.read(&mut self.buffer[self.buffer_fill..])?;
            self.buffer_fill += n_bytes;
            if n_bytes == 0 || self.buffer_fill == self.buffer.capacity() {
                self.finished = n_bytes == 0;
                break
            }
        }

        Ok(Some(self.buffer_fill))
    }
}




/// A deterministic, file-derived identifier for a paired read.
///
/// Every field is a pure function of the input byte stream, independent of thread count,
/// scheduling, and `buffer_size`:
/// - `ordinal` — the read pair's 0-based position in file order.
/// - `byte_offset1` / `byte_offset2` — the byte offset of the record's first byte (`@`) in
///   file1 / file2 respectively. On gzip input these are offsets into the *decompressed*
///   stream (still deterministic, but not directly seekable in the compressed file).
///
/// These are assigned from cumulative counters maintained under the shared reader's lock
/// (see [`FastqPairedByteReader::fill_buf_ids`]), so which worker thread processes a batch
/// never changes a read's id.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct ReadId {
    pub ordinal: u64,
    pub byte_offset1: u64,
    pub byte_offset2: u64,
    /// Dense, deterministic batch number (the chunk key for on-disk evidence). Batches are cut by
    /// the shared reader under its lock, so this is a pure function of the byte stream.
    pub batch_id: u64,
    /// Position within the batch (the group key inside a chunk). `ordinal == batch_base + idx_in_batch`.
    pub idx_in_batch: u64,
}

/// The id of a batch's first record, captured before the batch is drained. Per-record ids
/// are `base + within-batch offset`, computed by the worker thread (see
/// [`crate::fastq_reader::PairedFastqReader::next_id`]).
#[derive(Debug, Clone, Copy, Default)]
pub struct BatchBase {
    pub byte_offset1: u64,
    pub byte_offset2: u64,
    pub ordinal: u64,
    pub batch_id: u64,
}

pub struct FastqPairedByteReader<T>
where
    T: std::io::Read, {
    file1: T,
    file2: T,
    buffer2: Vec<u8>,
    buffer1: Vec<u8>,
    buffer1_fill: usize,
    buffer2_fill: usize,
    finished1: bool,
    finished2: bool,
    // Cumulative counters for read-id assignment, advanced only by `fill_buf_ids` (the plain
    // `fill_buf` path leaves them untouched, so it pays nothing). Each is the base for the
    // *next* batch: total bytes drained from file1/file2, and total records emitted, so far.
    file1_offset: u64,
    file2_offset: u64,
    record_offset: u64,
    batch_counter: u64,
}

impl<T: Read> FillBufferPair for FastqPairedByteReader<T> {
    fn fill_buf(
        &mut self,
        buffer1: &mut Vec<u8>,
        buffer2: &mut Vec<u8>,
    ) -> Result<Option<(usize, usize)>, Error> {
        Ok(self
            .cut_batch(buffer1, buffer2)?
            .map(|(pos1, pos2, _records)| (pos1, pos2)))
    }
}

impl<T: Read> FastqPairedByteReader<T> {
    /// Drains the next batch into `buffer1`/`buffer2` and reports how many read pairs it
    /// contains. This is the single copy of the batch-cutting logic; both the plain
    /// [`fill_buf`](FillBufferPair::fill_buf) and [`fill_buf_ids`](Self::fill_buf_ids) call
    /// it, so the two paths can never drift on where batches are cut. It does **not** touch
    /// the read-id counters — the plain path pays nothing for the id machinery.
    fn cut_batch(
        &mut self,
        buffer1: &mut Vec<u8>,
        buffer2: &mut Vec<u8>,
    ) -> Result<Option<(usize, usize, u64)>, Error> {
        self.read()?;

        if self.done() {
            return Ok(None);
        };

        let (lines, pos1, pos2) = self.byte_pos().expect("msg");

        assert!(lines > 0);

        assert_eq!(self.buffer1[0], b'@');
        buffer1[..pos1].copy_from_slice(&self.buffer1[..pos1]);
        self.buffer1.drain(0..pos1);
        self.buffer1_fill -= pos1;
        self.buffer1.resize(self.buffer1.capacity(), 0); // restore buffer fill
        assert_eq!(buffer1[0], b'@');


        assert_eq!(self.buffer2[0], b'@');
        buffer2[..pos2].copy_from_slice(&self.buffer2[..pos2]);
        self.buffer2.drain(0..pos2);
        self.buffer2_fill -= pos2;
        self.buffer2.resize(self.buffer2.capacity(), 0); // restore buffer fill
        assert_eq!(buffer2[0], b'@');

        // `lines` is the per-file line count for this batch (min across the two files); a
        // FASTQ record is exactly 4 lines, so both files emit `lines / 4` records here.
        Ok(Some((pos1, pos2, (lines / 4) as u64)))
    }

    /// Like [`fill_buf`](FillBufferPair::fill_buf), but also returns the [`BatchBase`] — the
    /// byte offsets and record ordinal of this batch's first record — and advances the
    /// cumulative counters. The base is captured *before* draining, and all counter updates
    /// happen here under the shared reader lock, so a batch's base is a deterministic
    /// function of the byte stream regardless of which thread calls this or how many do.
    pub fn fill_buf_ids(
        &mut self,
        buffer1: &mut Vec<u8>,
        buffer2: &mut Vec<u8>,
    ) -> Result<Option<(usize, usize, BatchBase)>, Error> {
        let base = BatchBase {
            byte_offset1: self.file1_offset,
            byte_offset2: self.file2_offset,
            ordinal: self.record_offset,
            batch_id: self.batch_counter,
        };

        match self.cut_batch(buffer1, buffer2)? {
            Some((pos1, pos2, records)) => {
                self.file1_offset += pos1 as u64;
                self.file2_offset += pos2 as u64;
                self.record_offset += records;
                self.batch_counter += 1;
                Ok(Some((pos1, pos2, base)))
            }
            None => Ok(None),
        }
    }
}

impl<T: Read> FastqPairedByteReader<T> {
    fn done(&mut self) -> bool {
        self.finished1 && self.finished2 && self.buffer1_fill == 0 && self.buffer2_fill == 0
    }

    fn byte_pos(&mut self) -> Option<(usize, usize, usize)> {
        let search = self.buffer1_fill > 0
            || self.buffer2_fill > 0;

        let mut lines = 0;
        let pos1;
        let pos2;


        if search {
            lines = self.load_lines()?;
            pos1 = memchr_iter(b'\n', &self.buffer1).take(lines).last()?;
            pos2 = memchr_iter(b'\n', &self.buffer2).take(lines).last()?;
        } else {
            assert!(self.buffer1_fill > 0);
            assert!(self.buffer2_fill > 0);
            pos1 = self.buffer1_fill - 1;
            pos2 = self.buffer2_fill - 1;
        }
        Some((lines, pos1 + 1, pos2 + 1))
    }

    fn load_lines(&mut self) -> Option<usize> {
        let lines1 = memchr_iter(b'\n', &self.buffer1).count() & !0b11;
        let lines2 = memchr_iter(b'\n', &self.buffer2).count() & !0b11;

        // println!("Lines: {} == {}", lines1, lines2);
        let min = std::cmp::min(lines1, lines2);
        if min > 0 {
            Some(min)
        } else {
            None
        }
    }

    pub fn new(
        file1: T,
        file2: T,
        buff_capacity: usize,
    ) -> Self {
        FastqPairedByteReader {
            file1: file1,
            file2: file2,
            buffer1_fill: 0,
            buffer2_fill: 0,
            buffer1: vec![0; buff_capacity],
            buffer2: vec![0; buff_capacity],
            finished1: false,
            finished2: false,
            file1_offset: 0,
            file2_offset: 0,
            record_offset: 0,
            batch_counter: 0,
        }
    }

    pub fn invalid(&mut self) -> bool {
        (self.buffer1_fill == 0 || self.buffer2_fill == 0) && self.buffer1_fill != self.buffer2_fill
    }

    pub fn fill_both_buffs(&mut self) -> Result<(), Error> {
        // println!("Fill buffs: {} {} {} {}", self.buffer1_fill, self.buffer2_fill, self.finished1, self.finished2);

        if self.buffer1_fill < self.buffer1.capacity() {
            // GzDecoder does not read full buffer but only chunks so we need to loop
            // to Fill the buffer
            loop {
                let n_bytes1 = self.file1.read(&mut self.buffer1[self.buffer1_fill..])?;
                self.buffer1_fill += n_bytes1;
                if n_bytes1 == 0 || self.buffer1_fill == self.buffer1.capacity() {
                    self.finished1 = n_bytes1 == 0;
                    break
                }
            }
            let n_bytes1 = self.file1.read(&mut self.buffer1[self.buffer1_fill..])?;
            self.buffer1_fill += n_bytes1;
        }

        if self.buffer2_fill < self.buffer2.capacity() {
            // GzDecoder does not read full buffer but only chunks so we need to loop
            // to Fill the buffer
            loop {
                let n_bytes2 = self.file2.read(&mut self.buffer2[self.buffer2_fill..])?;
                self.buffer2_fill += n_bytes2;
                if n_bytes2 == 0 || self.buffer2_fill == self.buffer2.capacity() {
                    self.finished2 = n_bytes2 == 0;
                    break
                }
            }
        }

        if self.done() { return Ok(()); };

        assert_eq!(self.buffer1[0], b'@');
        assert_eq!(self.buffer2[0], b'@');

        assert!(self.done() || (self.buffer1[0] == b'@' && self.buffer2[0] == b'@'));
        Ok(())
    }

    pub fn read(&mut self) -> Result<Option<()>, Error> {
        if self.finished1 && self.finished2 {
            return Ok(None);
        }

        let _ = self.fill_both_buffs();

        if self.invalid() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "Fastq files of different length.",
            ));
        }

        Ok(Some(()))
    }
}
