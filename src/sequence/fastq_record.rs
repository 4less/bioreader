
use core::slice::memchr;
use std::fmt::{self, Display};
use memchr::memchr;

/// Represents the position of a record within a buffer
#[derive(Debug, Clone, Default)] //, Serialize, Deserialize
pub struct BufferPosition {
    // (start, stop), but might include \r at the end
    pub pos: (usize, usize),
    pub seq: usize,
    pub sep: usize,
    pub qual: usize,
}

#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd)]
pub enum SearchPosition {
    Id,
    Sequence,
    Separator,
    Quality,
}

/// Remove a final '\r' from a byte slice
#[inline]
fn trim_cr(line: &[u8]) -> &[u8] {
    if let Some((&b'\r', remaining)) = line.split_last() {
        remaining
    } else {
        line
    }
}

/// Whether it uses \r\n or only \n
#[derive(Debug, PartialEq, Eq, Hash, Copy, Clone)]
pub enum LineEnding {
    Windows,
    Unix,
}

impl LineEnding {
    pub fn to_bytes(&self) -> Vec<u8> {
        match self {
            Self::Windows => vec![b'\r', b'\n'],
            Self::Unix => vec![b'\n'],
        }
    }
}

pub fn find_line_ending(bytes: &[u8]) -> Option<LineEnding> {
    if !bytes.is_empty() {
        if let Some(idx) = memchr(b'\n', bytes) {
            if idx > 0 && bytes[idx - 1] == b'\r' {
                return Some(LineEnding::Windows);
            }

            return Some(LineEnding::Unix);
        }
    }
    None
}

impl BufferPosition {
    #[inline]
    pub fn is_new(&self) -> bool {
        self.pos.1 == 0
    }

    #[inline]
    pub fn reset(&mut self, start: usize) {
        self.pos.0 = start;
        self.pos.1 = 0;
    }

    #[inline]
    pub fn head<'a>(&'a self, buffer: &'a [u8]) -> &'a [u8] {
        trim_cr(&buffer[self.pos.0 + 1..self.seq - 1])
    }

    #[inline]
    pub fn seq<'a>(&'a self, buffer: &'a [u8]) -> &'a [u8] {
        trim_cr(&buffer[self.seq..self.sep - 1])
    }

    #[inline]
    pub fn qual<'a>(&'a self, buffer: &'a [u8]) -> &'a [u8] {
        trim_cr(&buffer[self.qual..self.pos.1])
    }
}


/// A FASTQ record that borrows data from a buffer
#[derive(Debug, Clone)]
pub struct RefRecord<'a> {
    pub buffer: &'a [u8],
    pub buf_pos: &'a BufferPosition,
}

impl<'a> RefRecord<'a> { // Record for 
    #[inline]
    pub fn head(&self) -> &[u8] {
        self.buf_pos.head(self.buffer)
    }

    #[inline]
    pub fn seq(&self) -> &[u8] {
        self.buf_pos.seq(self.buffer)
    }

    #[inline]
    pub fn qual(&self) -> &[u8] {
        self.buf_pos.qual(self.buffer)
    }

    #[inline]
    pub fn valid(&self) -> bool {
        let mut valid = self.seq().len() == self.qual().len();
        valid &= self.seq().iter().all(|&c| c == b'A' || c == b'C' || c == b'G' || c == b'T' || c == b'N');
        return valid
    }

    #[inline]
    pub fn perfect(&self) -> bool {
        let mut valid = self.seq().len() == self.qual().len();
        valid &= self.seq().iter().all(|&c| c == b'A' || c == b'C' || c == b'G' || c == b'T');
        return valid
    }
}

impl<'a> Display for RefRecord<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}\n{}\n+\n{}", 
            std::str::from_utf8(self.head()).expect("Expect printable string"),
            std::str::from_utf8(self.seq()).expect("Expect printable string"), 
            std::str::from_utf8(self.qual()).expect("Expect printable string"))
    }
}