// pub enum FastaLoadBatch {
//     Resize(usize),
//     Good(usize),
//     Done(usize),
// }

// pub trait FastaFillBuffer {
//     fn fill_buf(&mut self, buffer: &mut Vec<u8>) -> Result<FastaLoadBatch, std::io::Error>;
// }

use crate::fastq_byte_reader::FillBuffer;


pub struct FastaByteReader<T>
where
    T: std::io::Read,
{
    file: T,
    buffer_fill: usize,
    buffer: Vec<u8>,
    finished: bool,
}

impl<T: std::io::Read> FillBuffer for FastaByteReader<T> {
    fn fill_buf(&mut self, buf: &mut Vec<u8>) -> Result<Option<usize>, std::io::Error> {
        // Implies that file has been read to the end, but there still might be data in the buffer
        if self.finished {
            if self.buffer_fill == 0 {
                return Ok(None);
            } else if self.buffer_fill <= buf.len() {
                buf[..self.buffer_fill].copy_from_slice(&self.buffer[..self.buffer_fill]);
                let fill = self.buffer_fill;
                self.buffer_fill = 0;
                return Ok(Some(fill));
            } // File has been read and buffer is empty
        }

        // Find record in slice, at this point buffer must be filled to the brim
        // assert!(self.buffer_fill >= buf.len());


        // Find end of last complete Fastq record in local buffer
        let mut index = self.find_next(&self.buffer[..self.buffer_fill]);

        while index == 0 && !self.finished {
            self.buffer.resize(self.buffer.len() * 2, 0);
            self.read_file()?;
            index = self.find_next(&self.buffer[..self.buffer_fill]);
        }
        buf.resize(self.buffer.len(), 0);


        let buffer_slice = &self.buffer[..self.buffer_fill];

        assert!(index < buf.len());

        // Copy local buffer of complete Fastq records into external buffer
        buf[..index].copy_from_slice(&buffer_slice[..index]);
        // println!("- Bytes: {index} {}", buf.len());
        // println!("- End: {}", std::str::from_utf8(&buf[index-200..index]).unwrap());

        self.buffer.drain(0..index);
        self.buffer_fill -= index; // unprocessed bytes
        self.buffer.resize(self.buffer.capacity(), 0); // restore buffer fill


        if self.buffer_fill < self.buffer.len() {
            self.read_file()?;
        }
        
        return Ok(Some(index))
    }
}

impl<T: std::io::Read> FastaByteReader<T> {
    pub fn new(reader: T, chunk_size: usize) -> Result<Self, std::io::Error> {
        let mut br = Self {
            file: reader,
            buffer_fill: 0,
            buffer: vec![0; chunk_size],
            finished: false,
        };
        br.read_file()?;

        Ok(br)
    }

    pub fn find_next(&self, buffer_slice: &[u8]) -> usize {
        for (i, &c) in buffer_slice.iter().enumerate().rev() {
            if i > 1 && c == b'>' && buffer_slice[i - 1] == b'\n' {
                return i
            }
        }
        0
    }

    pub fn read_file(&mut self) -> std::io::Result<usize> {
        let read_bytes = self.file.read(&mut self.buffer[self.buffer_fill..])?;

        self.finished = read_bytes == 0;
        self.buffer_fill += read_bytes;
        Ok(read_bytes)
    }
}