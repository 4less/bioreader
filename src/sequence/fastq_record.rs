
use core::slice::memchr;
use std::{fmt::{self, Display}, u8};
use memchr::memchr;
use colored::{Color, ColoredString, Colorize, CustomColor};

use super::utils::reverse_complement_into_vec;

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

#[derive(Clone)]
pub struct OwnedFastqRecord {
    pub header: Vec<u8>,
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
}

impl OwnedFastqRecord {
    pub fn new() -> Self {
        Self {
            header: Vec::new(),
            sequence: Vec::new(),
            quality: Vec::new(),
        }
    }

    pub fn head(&self) -> &[u8] {
        &self.header
    }

    pub fn seq(&self) -> &[u8] {
        &self.sequence
    }

    pub fn qual(&self) -> &[u8] {
        &self.quality
    }
}

impl Display for OwnedFastqRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}\n{}\n+\n{}", 
            std::str::from_utf8(self.head()).expect("Expect printable string"),
            std::str::from_utf8(self.seq()).expect("Expect printable string"), 
            std::str::from_utf8(self.qual()).expect("Expect printable string"))
    }
}

pub fn qual_to_color(q: u8) -> CustomColor {
    match q {
        0..10 => CustomColor::new(255,0,255 - q*25),
        10..20 => CustomColor::new(255,0 + (q-10)*25,0),
        20..30 => CustomColor::new(255 - (q-20)*25,255,0),
        _ => CustomColor::new(0,255,0),
    }
}

pub fn print_color_qualities(qual: &[u8], qual_offset: Option<u8>) {
    let offset = qual_offset.unwrap_or(0);
    
    for q in qual {
        let char = *q as char;
        eprint!("{}", char.to_string().custom_color(qual_to_color(q - offset)));
    }
    eprintln!();
}

/// A FASTQ record that borrows data from a buffer
#[derive(Debug, Clone)]
pub struct RefFastqRecord<'a> {
    pub buffer: &'a [u8],
    pub buf_pos: &'a BufferPosition,
}

impl<'a> RefFastqRecord<'a> { // Record for 
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
    pub fn valid_extended(&self) -> bool {
        let valid = self.seq().iter().all(
            |&c| { 
                c == b'A' || c == b'C' || c == b'G' || c == b'T' || c == b'N' || 
                c == b'R' || c == b'Y' || c == b'K' || c == b'M' || c == b'S' ||
                c == b'W' || c == b'B' || c == b'D' || c == b'H' || c == b'V'
            }

        );
        return valid
    }

    #[inline]
    pub fn perfect(&self) -> bool {
        let mut valid = self.seq().len() == self.qual().len();
        valid &= self.seq().iter().all(|&c| c == b'A' || c == b'C' || c == b'G' || c == b'T');
        return valid
    }

    #[inline]
    pub fn reverse_complement(&self, rec: &mut OwnedFastqRecord) -> () {
        rec.header.clear();
        rec.sequence.clear();
        rec.quality.clear();

        rec.header.extend_from_slice(self.head());
        rec.quality.extend_from_slice(self.qual());
        reverse_complement_into_vec(self.seq(), &mut rec.sequence);
        rec.quality.reverse();

    }
}

impl<'a> Display for RefFastqRecord<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}\n{}\n+\n{}", 
            std::str::from_utf8(self.head()).expect("Expect printable string"),
            std::str::from_utf8(self.seq()).expect("Expect printable string"), 
            std::str::from_utf8(self.qual()).expect("Expect printable string"))
    }
}