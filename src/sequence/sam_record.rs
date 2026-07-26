//! SAM alignment records that borrow from a batch buffer.
//!
//! A SAM record is a single line of tab separated fields: 11 mandatory ones
//! (QNAME..QUAL) followed by an arbitrary number of optional `TAG:TYPE:VALUE` fields.
//! Like [`RefFastqRecord`](crate::sequence::fastq_record::RefFastqRecord), a record never
//! owns its bytes — [`SamBufferPosition`] stores the field offsets within the reader's batch
//! buffer and every accessor is a slice into that buffer.

use std::fmt::{self, Display};

use memchr::memchr;

/// Number of mandatory SAM fields (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT,
/// TLEN, SEQ, QUAL).
pub const MANDATORY_FIELDS: usize = 11;

/// FLAG field bits (SAM spec section 1.4.2).
pub mod flags {
    pub const PAIRED: u16 = 0x1;
    pub const PROPER_PAIR: u16 = 0x2;
    pub const UNMAPPED: u16 = 0x4;
    pub const MATE_UNMAPPED: u16 = 0x8;
    pub const REVERSE: u16 = 0x10;
    pub const MATE_REVERSE: u16 = 0x20;
    pub const FIRST_IN_PAIR: u16 = 0x40;
    pub const LAST_IN_PAIR: u16 = 0x80;
    pub const SECONDARY: u16 = 0x100;
    pub const QC_FAIL: u16 = 0x200;
    pub const DUPLICATE: u16 = 0x400;
    pub const SUPPLEMENTARY: u16 = 0x800;
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

/// Parse an unsigned decimal number. `None` on empty input or any non-digit byte, so callers
/// can distinguish a malformed field from a zero.
#[inline]
pub fn parse_u64(bytes: &[u8]) -> Option<u64> {
    if bytes.is_empty() {
        return None;
    }
    let mut value: u64 = 0;
    for &c in bytes {
        if !c.is_ascii_digit() {
            return None;
        }
        value = value.wrapping_mul(10).wrapping_add((c - b'0') as u64);
    }
    Some(value)
}

/// Parse a signed decimal number (TLEN is the only signed mandatory field).
#[inline]
pub fn parse_i64(bytes: &[u8]) -> Option<i64> {
    match bytes.split_first() {
        Some((&b'-', rest)) => parse_u64(rest).map(|v| -(v as i64)),
        Some((&b'+', rest)) => parse_u64(rest).map(|v| v as i64),
        _ => parse_u64(bytes).map(|v| v as i64),
    }
}

/// Position of a record within a batch buffer.
///
/// `pos` is the half-open byte range of the record line, with the line terminator (and a
/// Windows `\r`) excluded. `fields[i]` is the start offset of mandatory field `i`; the extra
/// trailing entry `fields[MANDATORY_FIELDS]` is the end offset of QUAL — either the position
/// of the tab that introduces the optional fields, or `pos.1` when there are none.
#[derive(Debug, Clone)]
pub struct SamBufferPosition {
    pub pos: (usize, usize),
    pub fields: [usize; MANDATORY_FIELDS + 1],
}

impl Default for SamBufferPosition {
    fn default() -> Self {
        Self {
            pos: (0, 0),
            fields: [0; MANDATORY_FIELDS + 1],
        }
    }
}

impl SamBufferPosition {
    /// Record the field offsets of the line `buffer[start..end]` (terminator already trimmed).
    ///
    /// Returns `false` if the line has fewer than [`MANDATORY_FIELDS`] fields, leaving `self`
    /// in an unspecified state — the caller decides whether that is an error or a line to skip.
    pub fn parse(&mut self, buffer: &[u8], start: usize, end: usize) -> bool {
        self.pos = (start, end);

        let mut cursor = start;
        for i in 0..MANDATORY_FIELDS {
            self.fields[i] = cursor;
            match memchr(b'\t', &buffer[cursor..end]) {
                Some(tab) => cursor += tab + 1,
                None if i + 1 == MANDATORY_FIELDS => {
                    // QUAL is the last field on the line: no optional fields follow.
                    cursor = end + 1;
                }
                None => return false,
            }
        }
        // `cursor` sits one past the tab that ended QUAL, so `cursor - 1` is that tab — or
        // `end` when QUAL ran to the end of the line.
        self.fields[MANDATORY_FIELDS] = cursor - 1;

        true
    }

    /// The `i`-th mandatory field, without its separator.
    #[inline]
    pub fn field<'a>(&self, buffer: &'a [u8], i: usize) -> &'a [u8] {
        debug_assert!(i < MANDATORY_FIELDS);
        let end = if i + 1 == MANDATORY_FIELDS {
            self.fields[MANDATORY_FIELDS]
        } else {
            self.fields[i + 1] - 1
        };
        &buffer[self.fields[i]..end]
    }

    /// The raw optional-field region (`TAG:TYPE:VALUE` separated by tabs), empty if the record
    /// carries no optional fields.
    #[inline]
    pub fn optional<'a>(&self, buffer: &'a [u8]) -> &'a [u8] {
        let qual_end = self.fields[MANDATORY_FIELDS];
        if qual_end < self.pos.1 {
            &buffer[qual_end + 1..self.pos.1]
        } else {
            &[]
        }
    }

    /// The full record line.
    #[inline]
    pub fn line<'a>(&self, buffer: &'a [u8]) -> &'a [u8] {
        &buffer[self.pos.0..self.pos.1]
    }
}

/// A SAM record that borrows data from a buffer.
#[derive(Debug, Clone)]
pub struct RefSamRecord<'a> {
    pub buffer: &'a [u8],
    pub buf_pos: &'a SamBufferPosition,
}

impl<'a> RefSamRecord<'a> {
    #[inline]
    pub fn field(&self, i: usize) -> &[u8] {
        self.buf_pos.field(self.buffer, i)
    }

    #[inline]
    pub fn qname(&self) -> &[u8] {
        self.field(0)
    }

    #[inline]
    pub fn rname(&self) -> &[u8] {
        self.field(2)
    }

    #[inline]
    pub fn cigar(&self) -> &[u8] {
        self.field(5)
    }

    #[inline]
    pub fn rnext(&self) -> &[u8] {
        self.field(6)
    }

    #[inline]
    pub fn seq(&self) -> &[u8] {
        self.field(9)
    }

    #[inline]
    pub fn qual(&self) -> &[u8] {
        trim_cr(self.field(10))
    }

    /// The optional fields as one raw, tab separated slice. See [`tags`](Self::tags) to walk
    /// them one by one.
    #[inline]
    pub fn optional(&self) -> &[u8] {
        trim_cr(self.buf_pos.optional(self.buffer))
    }

    /// The record as it appears in the file, without the line terminator.
    #[inline]
    pub fn as_bytes(&self) -> &[u8] {
        self.buf_pos.line(self.buffer)
    }

    #[inline]
    pub fn try_flag(&self) -> Option<u16> {
        parse_u64(self.field(1)).filter(|&v| v <= u16::MAX as u64).map(|v| v as u16)
    }

    #[inline]
    pub fn try_pos(&self) -> Option<u64> {
        parse_u64(self.field(3))
    }

    #[inline]
    pub fn try_mapq(&self) -> Option<u8> {
        parse_u64(self.field(4)).filter(|&v| v <= u8::MAX as u64).map(|v| v as u8)
    }

    #[inline]
    pub fn try_pnext(&self) -> Option<u64> {
        parse_u64(self.field(7))
    }

    #[inline]
    pub fn try_tlen(&self) -> Option<i64> {
        parse_i64(self.field(8))
    }

    #[inline]
    pub fn flag(&self) -> u16 {
        self.try_flag().expect("FLAG is not a valid number")
    }

    /// 1-based leftmost mapping position; 0 for an unmapped read.
    #[inline]
    pub fn pos(&self) -> u64 {
        self.try_pos().expect("POS is not a valid number")
    }

    #[inline]
    pub fn mapq(&self) -> u8 {
        self.try_mapq().expect("MAPQ is not a valid number")
    }

    #[inline]
    pub fn pnext(&self) -> u64 {
        self.try_pnext().expect("PNEXT is not a valid number")
    }

    #[inline]
    pub fn tlen(&self) -> i64 {
        self.try_tlen().expect("TLEN is not a valid number")
    }

    #[inline]
    pub fn has_flag(&self, mask: u16) -> bool {
        self.flag() & mask != 0
    }

    #[inline]
    pub fn is_paired(&self) -> bool {
        self.has_flag(flags::PAIRED)
    }

    #[inline]
    pub fn is_proper_pair(&self) -> bool {
        self.has_flag(flags::PROPER_PAIR)
    }

    #[inline]
    pub fn is_unmapped(&self) -> bool {
        self.has_flag(flags::UNMAPPED)
    }

    #[inline]
    pub fn is_mate_unmapped(&self) -> bool {
        self.has_flag(flags::MATE_UNMAPPED)
    }

    #[inline]
    pub fn is_reverse(&self) -> bool {
        self.has_flag(flags::REVERSE)
    }

    #[inline]
    pub fn is_first_in_pair(&self) -> bool {
        self.has_flag(flags::FIRST_IN_PAIR)
    }

    #[inline]
    pub fn is_last_in_pair(&self) -> bool {
        self.has_flag(flags::LAST_IN_PAIR)
    }

    #[inline]
    pub fn is_secondary(&self) -> bool {
        self.has_flag(flags::SECONDARY)
    }

    #[inline]
    pub fn is_supplementary(&self) -> bool {
        self.has_flag(flags::SUPPLEMENTARY)
    }

    #[inline]
    pub fn is_duplicate(&self) -> bool {
        self.has_flag(flags::DUPLICATE)
    }

    /// A primary alignment: neither secondary nor supplementary.
    #[inline]
    pub fn is_primary(&self) -> bool {
        self.flag() & (flags::SECONDARY | flags::SUPPLEMENTARY) == 0
    }

    /// Iterator over the optional `TAG:TYPE:VALUE` fields.
    #[inline]
    pub fn tags(&self) -> SamTagIter<'_> {
        SamTagIter {
            rest: self.optional(),
        }
    }

    /// The first optional field with the given two-character tag, e.g. `record.tag(b"NM")`.
    #[inline]
    pub fn tag(&self, tag: &[u8; 2]) -> Option<SamTag<'_>> {
        self.tags().find(|t| t.tag == *tag)
    }

    /// All numeric fields parse, QNAME is non-empty, and SEQ/QUAL agree in length unless
    /// either is absent (`*`).
    pub fn valid(&self) -> bool {
        if self.qname().is_empty() || self.seq().is_empty() || self.qual().is_empty() {
            return false;
        }
        if self.try_flag().is_none()
            || self.try_pos().is_none()
            || self.try_mapq().is_none()
            || self.try_pnext().is_none()
            || self.try_tlen().is_none()
        {
            return false;
        }
        self.seq() == b"*" || self.qual() == b"*" || self.seq().len() == self.qual().len()
    }

    /// [`valid`](Self::valid) plus a check that SEQ contains only IUPAC nucleotide codes
    /// (or is absent). Mirrors `RefFastqRecord::valid_extended`.
    pub fn valid_extended(&self) -> bool {
        if !self.valid() {
            return false;
        }
        let seq = self.seq();
        seq == b"*"
            || seq.iter().all(|&c| {
                c == b'A' || c == b'C' || c == b'G' || c == b'T' || c == b'N' ||
                c == b'R' || c == b'Y' || c == b'K' || c == b'M' || c == b'S' ||
                c == b'W' || c == b'B' || c == b'D' || c == b'H' || c == b'V' ||
                c == b'='
            })
    }
}

impl<'a> Display for RefSamRecord<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", String::from_utf8_lossy(self.as_bytes()))
    }
}

/// One optional SAM field, split into its `TAG`, single-character `TYPE` and raw `VALUE`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SamTag<'a> {
    pub tag: [u8; 2],
    pub tag_type: u8,
    pub value: &'a [u8],
}

impl<'a> SamTag<'a> {
    /// The value of an integer tag (`i`).
    #[inline]
    pub fn as_int(&self) -> Option<i64> {
        (self.tag_type == b'i').then(|| parse_i64(self.value)).flatten()
    }

    /// The value of a character (`A`), string (`Z`) or hex (`H`) tag.
    #[inline]
    pub fn as_bytes(&self) -> &'a [u8] {
        self.value
    }
}

/// Iterator over the optional fields of a [`RefSamRecord`]. Fields that are not shaped like
/// `TAG:TYPE:VALUE` are skipped.
pub struct SamTagIter<'a> {
    rest: &'a [u8],
}

impl<'a> Iterator for SamTagIter<'a> {
    type Item = SamTag<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        while !self.rest.is_empty() {
            let (field, remainder) = match memchr(b'\t', self.rest) {
                Some(tab) => (&self.rest[..tab], &self.rest[tab + 1..]),
                None => (self.rest, &self.rest[self.rest.len()..]),
            };
            self.rest = remainder;

            if field.len() >= 5 && field[2] == b':' && field[4] == b':' {
                return Some(SamTag {
                    tag: [field[0], field[1]],
                    tag_type: field[3],
                    value: &field[5..],
                });
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const LINE: &[u8] =
        b"read1\t99\tchr1\t1000\t60\t8M\t=\t1200\t208\tACGTACGT\tIIIIIIII\tNM:i:2\tMD:Z:8\tAS:i:-3";

    fn parse(line: &[u8]) -> SamBufferPosition {
        let mut buf_pos = SamBufferPosition::default();
        assert!(buf_pos.parse(line, 0, line.len()));
        buf_pos
    }

    #[test]
    fn parses_mandatory_fields() {
        let buf_pos = parse(LINE);
        let record = RefSamRecord { buffer: LINE, buf_pos: &buf_pos };

        assert_eq!(record.qname(), b"read1");
        assert_eq!(record.flag(), 99);
        assert_eq!(record.rname(), b"chr1");
        assert_eq!(record.pos(), 1000);
        assert_eq!(record.mapq(), 60);
        assert_eq!(record.cigar(), b"8M");
        assert_eq!(record.rnext(), b"=");
        assert_eq!(record.pnext(), 1200);
        assert_eq!(record.tlen(), 208);
        assert_eq!(record.seq(), b"ACGTACGT");
        assert_eq!(record.qual(), b"IIIIIIII");
        assert!(record.valid_extended());
        assert!(record.is_paired() && record.is_first_in_pair() && !record.is_unmapped());
    }

    #[test]
    fn parses_optional_fields() {
        let buf_pos = parse(LINE);
        let record = RefSamRecord { buffer: LINE, buf_pos: &buf_pos };

        let tags: Vec<_> = record.tags().collect();
        assert_eq!(tags.len(), 3);
        assert_eq!(record.tag(b"NM").unwrap().as_int(), Some(2));
        assert_eq!(record.tag(b"MD").unwrap().as_bytes(), b"8");
        assert_eq!(record.tag(b"AS").unwrap().as_int(), Some(-3));
        assert!(record.tag(b"XX").is_none());
    }

    #[test]
    fn record_without_optional_fields() {
        let line: &[u8] = b"read2\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*";
        let buf_pos = parse(line);
        let record = RefSamRecord { buffer: line, buf_pos: &buf_pos };

        assert_eq!(record.qname(), b"read2");
        assert_eq!(record.qual(), b"*");
        assert_eq!(record.optional(), b"");
        assert_eq!(record.tags().count(), 0);
        assert!(record.is_unmapped());
        assert!(record.valid());
    }

    #[test]
    fn truncated_record_is_rejected() {
        let line: &[u8] = b"read3\t0\tchr1\t1\t0\t*";
        let mut buf_pos = SamBufferPosition::default();
        assert!(!buf_pos.parse(line, 0, line.len()));
    }

    #[test]
    fn seq_and_qual_length_mismatch_is_invalid() {
        let line: &[u8] = b"read4\t0\tchr1\t1\t0\t4M\t*\t0\t0\tACGT\tIII";
        let buf_pos = parse(line);
        let record = RefSamRecord { buffer: line, buf_pos: &buf_pos };
        assert!(!record.valid());
    }
}
