use flate2::bufread::GzDecoder;
use memchr::{memchr_iter, memmem::Finder};
use memmap2::Mmap;
use std::{
    fs::File,
    io::{BufReader, Error, Read},
    path::Path,
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

pub struct ByteReader<T>
where
    T: std::io::Read,
{
    file: T,
    buffer_fill: usize,
    buffer: Vec<u8>,
    finished: bool,
}

impl<T: std::io::Read> std::io::Read for ByteReader<T> {
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

impl<T: std::io::Read> ByteReader<T> {
    pub fn new(reader: T, chunk_size: usize) -> Result<ByteReader<T>, Error> {
        let mut br = ByteReader {
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

pub struct PairedByteReader<T>
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
}

impl<T: Read> FillBufferPair for PairedByteReader<T> {
    fn fill_buf(
        &mut self,
        buffer1: &mut Vec<u8>,
        buffer2: &mut Vec<u8>,
    ) -> Result<Option<(usize, usize)>, Error> {

        // let fr1 = self.buffer1_fill as f64/ self.buffer1.capacity() as f64;
        // let fr2 = self.buffer2_fill as f64/ self.buffer2.capacity() as f64;
        // if fr1 > fr2*2f64 {
        //     self.buffer2.resize(std::cmp::min(self.buffer2.capacity()*2, 536870912), 0);
        //     println!("Resize Buffer 2 {}", self.buffer2.capacity());
        //     println!("{} {}", fr1, fr2);
        //     println!("{} {}", self.buffer1_fill, self.buffer2_fill);
        //     println!("{} {}", self.buffer1.capacity(), self.buffer2.capacity());
        // }
        // if fr2*2f64 > fr1 {
        //     self.buffer1.resize(std::cmp::min(self.buffer1.capacity()*2, 536870912), 0);
        //     println!("Resize Buffer 1 {}", std::cmp::min(self.buffer1.capacity(), 536870912));
        //     println!("{} {}", fr1, fr2);
        //     println!("{} {}", self.buffer1_fill, self.buffer2_fill);
        // }

        self.read()?;


        if self.done() {
            return Ok(None);
        };

        let (_lines, pos1, pos2) = self.byte_pos().expect("msg");


        // println!("Lines: {}, pos1 {}, pos2 {}", _lines, pos1, pos2);
        // println!("Fill1 {}, Fill2 {}", self.buffer1_fill, self.buffer2_fill);
        assert!(_lines > 0);
        // let f1 = self.buffer1_fill;
        // println!("1After---{}\n---", std::str::from_utf8(&self.buffer1[f1-200..f1]).unwrap());

        // let f2 = self.buffer2_fill;
        // println!("2After---{}\n---", std::str::from_utf8(&self.buffer2[f2-200..f2]).unwrap());


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


        Ok(Some((pos1, pos2)))
    }
}

impl<T: Read> PairedByteReader<T> {
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
        PairedByteReader {
            file1: file1,
            file2: file2,
            buffer1_fill: 0,
            buffer2_fill: 0,
            buffer1: vec![0; buff_capacity],
            buffer2: vec![0; buff_capacity],
            finished1: false,
            finished2: false,
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
