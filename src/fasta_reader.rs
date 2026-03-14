use std::sync::{Arc, Mutex};
use crate::{fastq_byte_reader::FillBuffer, sequence::fasta_record::OwnedFastaRecord};
use memchr::memchr;


#[derive(Debug, Clone)]
pub struct FastaReader {
    pub buffer: Vec<u8>,
    pub buffer_pos: usize,
    pub buffer_fill: usize,
}

impl FastaReader {
    pub fn new(capacity: usize) -> Self {
        Self {
            buffer: vec![0; capacity],
            buffer_pos: 0,
            buffer_fill: 0,
        }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self::new(capacity)
    }

    #[inline]
    pub fn load_batch_par(&mut self, br: &mut Arc<Mutex<impl FillBuffer>>) -> Result<Option<()>, std::io::Error> {
        self.buffer_pos = 0;

        let bytes = match br.lock().expect("Locking ByteReader was unsuccessful").fill_buf(&mut self.buffer)? {
            Some(bytes) => bytes,
            None => 0,
        };

        self.buffer_fill = bytes;

        if bytes == 0 {
            return Ok(None)
        }

        Ok(Some(()))
     }

    pub fn next(&mut self, record: &mut OwnedFastaRecord) -> Option<()> {
        if self.buffer_pos >= self.buffer_fill {
            return None
        }

        assert!(self.buffer[self.buffer_pos] == b'>');
        let header_start: usize = self.buffer_pos;
        let header_length: usize = memchr(b'\n', &self.buffer[self.buffer_pos..]).expect("1: Couldn't find newline");
        self.buffer_pos += header_length + 1;

        let mut length = 0;

        record.clear();
        record.header.extend_from_slice(&self.buffer[header_start..header_start + header_length]);

        while self.buffer_pos < self.buffer_fill - 1 && self.buffer[self.buffer_pos] != b'>' {
            let newline_pos = match memchr(b'\n', &self.buffer[self.buffer_pos..]) {
                Some(pos) => pos,
                None => {
                    self.buffer_fill - self.buffer_pos
                },
            };

            length = newline_pos;

            if length > record.sequence.len() {
                record.sequence.reserve(length);
            }

            record.sequence.extend_from_slice(&self.buffer[self.buffer_pos..self.buffer_pos+newline_pos]);

            self.buffer_pos += newline_pos + 1;
        }

        Some(())
    }
}