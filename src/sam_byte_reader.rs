//! Shared, chunk-cutting byte reader for SAM files.
//!
//! Same shape as [`FastqByteReader`](crate::fastq_byte_reader::FastqByteReader) and
//! [`FastaByteReader`](crate::fasta_byte_reader::FastaByteReader): one instance sits behind an
//! `Arc<Mutex<..>>`, every worker thread pulls a batch of *complete records* out of it and then
//! parses that batch on its own, lock-free. SAM is line based, so a batch simply ends at the
//! last newline that fits in the buffer.
//!
//! The `@`-prefixed header lines are consumed once, at construction, and kept in
//! [`SamHeader`] — worker batches therefore contain alignment lines only.
//!
//! **Batch size matters more than anything else here.** Every byte is copied through a buffer
//! of `chunk_size`, so whether that buffer stays in cache dominates throughput: `benchmarks/`
//! measures 256 KiB–1 MiB as the sweet spot, with 16 MiB costing 25% on short reads and 2.3x
//! on long reads. Prefer `1 << 20`; bigger is not better.

use memchr::{memchr, memchr_iter, memrchr};
use std::io::{Error, Read};

use crate::fastq_byte_reader::FillBuffer;

/// The raw header block of a SAM file (all leading `@` lines, terminators included).
#[derive(Debug, Clone, Default)]
pub struct SamHeader {
    raw: Vec<u8>,
}

impl SamHeader {
    #[inline]
    pub fn raw(&self) -> &[u8] {
        &self.raw
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.raw.is_empty()
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.raw.len()
    }

    pub fn text(&self) -> std::borrow::Cow<'_, str> {
        String::from_utf8_lossy(&self.raw)
    }

    /// The header lines, without their terminators.
    pub fn lines(&self) -> impl Iterator<Item = &[u8]> {
        self.raw.split(|&c| c == b'\n').filter_map(|line| {
            let line = match line.split_last() {
                Some((&b'\r', rest)) => rest,
                _ => line,
            };
            (!line.is_empty()).then_some(line)
        })
    }

    /// The reference sequences declared by the `@SQ` lines, as `(SN, LN)` pairs. Lines missing
    /// either tag are skipped.
    pub fn reference_sequences(&self) -> Vec<(&[u8], u64)> {
        self.lines()
            .filter(|line| line.starts_with(b"@SQ"))
            .filter_map(|line| {
                let mut name = None;
                let mut length = None;
                for field in line.split(|&c| c == b'\t') {
                    if let Some(value) = field.strip_prefix(b"SN:".as_slice()) {
                        name = Some(value);
                    } else if let Some(value) = field.strip_prefix(b"LN:".as_slice()) {
                        length = crate::sequence::sam_record::parse_u64(value);
                    }
                }
                Some((name?, length?))
            })
            .collect()
    }
}

/// A deterministic, file-derived identifier for a SAM record.
///
/// Every field is a pure function of the input byte stream, independent of thread count,
/// scheduling and `buffer_size`:
/// - `ordinal` — the record's 0-based position in file order (header lines excluded).
/// - `byte_offset` — the byte offset of the record's first byte in the file, header
///   included. On gzip input this is an offset into the *decompressed* stream (still
///   deterministic, but not directly seekable in the compressed file).
///
/// These are assigned from cumulative counters maintained under the shared reader's lock
/// (see [`SamByteReader::fill_buf_ids`]), so which worker thread processes a batch never
/// changes a record's id.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct SamReadId {
    pub ordinal: u64,
    pub byte_offset: u64,
    /// Dense, deterministic batch number. Batches are cut by the shared reader under its
    /// lock, so this is a pure function of the byte stream.
    pub batch_id: u64,
    /// Position within the batch. `ordinal == batch_base + idx_in_batch`.
    pub idx_in_batch: u64,
}

/// The id of a batch's first record, captured before the batch is drained. Per-record ids are
/// `base + within-batch offset`, computed by the worker thread (see
/// [`crate::sam_reader::SamReader::next_id`]).
#[derive(Debug, Clone, Copy, Default)]
pub struct SamBatchBase {
    pub byte_offset: u64,
    pub ordinal: u64,
    pub batch_id: u64,
}

pub struct SamByteReader<T>
where
    T: Read,
{
    file: T,
    buffer: Vec<u8>,
    buffer_fill: usize,
    finished: bool,
    header: SamHeader,
    // Cumulative counters for read-id assignment, advanced only by `fill_buf_ids` (the plain
    // `fill_buf` path leaves them untouched, so it pays nothing). Each is the base for the
    // *next* batch: total bytes drained from the file (header included, hence the header-length
    // initialisation), and total records emitted, so far.
    file_offset: u64,
    record_offset: u64,
    batch_counter: u64,
}

impl<T: Read> FillBuffer for SamByteReader<T> {
    fn fill_buf(&mut self, buffer: &mut Vec<u8>) -> Result<Option<usize>, Error> {
        self.cut_batch(buffer)
    }
}

impl<T: Read> SamByteReader<T> {
    /// Reads the header and leaves the reader positioned at the first alignment line.
    pub fn new(reader: T, chunk_size: usize) -> Result<Self, Error> {
        let mut br = Self {
            file: reader,
            buffer: vec![0; chunk_size],
            buffer_fill: 0,
            finished: false,
            header: SamHeader::default(),
            file_offset: 0,
            record_offset: 0,
            batch_counter: 0,
        };
        br.read_header()?;
        br.file_offset = br.header.len() as u64;

        Ok(br)
    }

    #[inline]
    pub fn header(&self) -> &SamHeader {
        &self.header
    }

    /// Drains the next batch of complete records into `buffer` and returns its size in bytes.
    /// This is the single copy of the batch-cutting logic; both the plain
    /// [`fill_buf`](FillBuffer::fill_buf) and [`fill_buf_ids`](Self::fill_buf_ids) call it, so
    /// the two paths can never drift on where batches are cut. It does **not** touch the
    /// read-id counters — the plain path pays nothing for the id machinery.
    fn cut_batch(&mut self, buffer: &mut Vec<u8>) -> Result<Option<usize>, Error> {
        self.fill_buffer()?;

        if self.buffer_fill == 0 {
            return Ok(None);
        }

        let bytes = loop {
            match memrchr(b'\n', &self.buffer[..self.buffer_fill]) {
                Some(pos) => break pos + 1,
                // No complete line in the buffer: either the file ends without a trailing
                // newline, or a single record is longer than the buffer and it has to grow.
                None if self.finished => break self.buffer_fill,
                None => {
                    self.buffer.resize(self.buffer.len() * 2, 0);
                    self.fill_buffer()?;
                }
            }
        };

        if buffer.len() < bytes {
            buffer.resize(bytes, 0);
        }
        buffer[..bytes].copy_from_slice(&self.buffer[..bytes]);
        self.consume(bytes);

        Ok(Some(bytes))
    }

    /// Like [`fill_buf`](FillBuffer::fill_buf), but also returns the [`SamBatchBase`] — the byte
    /// offset and record ordinal of this batch's first record — and advances the cumulative
    /// counters. The base is captured *before* draining, and all counter updates happen here
    /// under the shared reader lock, so a batch's base is a deterministic function of the byte
    /// stream regardless of which thread calls this or how many do.
    pub fn fill_buf_ids(
        &mut self,
        buffer: &mut Vec<u8>,
    ) -> Result<Option<(usize, SamBatchBase)>, Error> {
        let base = SamBatchBase {
            byte_offset: self.file_offset,
            ordinal: self.record_offset,
            batch_id: self.batch_counter,
        };

        match self.cut_batch(buffer)? {
            Some(bytes) => {
                self.file_offset += bytes as u64;
                self.record_offset += count_records(&buffer[..bytes]);
                self.batch_counter += 1;
                Ok(Some((bytes, base)))
            }
            None => Ok(None),
        }
    }

    /// Consumes the leading `@` lines. The header can be far larger than one buffer (a
    /// reference with many contigs), so this drains whole buffers' worth of header lines at a
    /// time, growing the buffer only if a single header line does not fit.
    fn read_header(&mut self) -> Result<(), Error> {
        loop {
            self.fill_buffer()?;
            if self.buffer_fill == 0 {
                return Ok(());
            }

            // Walk complete `@` lines in one pass, so the buffer is drained once per fill
            // rather than once per header line.
            let mut pos = 0;
            while pos < self.buffer_fill && self.buffer[pos] == b'@' {
                match memchr(b'\n', &self.buffer[pos..self.buffer_fill]) {
                    Some(newline) => pos += newline + 1,
                    None => break,
                }
            }

            self.header.raw.extend_from_slice(&self.buffer[..pos]);
            let drained_all = pos == self.buffer_fill;
            self.consume(pos);

            if !drained_all {
                if self.buffer[0] != b'@' {
                    return Ok(()); // First alignment line reached.
                }
                if self.finished {
                    // Final header line without a trailing newline.
                    self.header.raw.extend_from_slice(&self.buffer[..self.buffer_fill]);
                    let rest = self.buffer_fill;
                    self.consume(rest);
                    return Ok(());
                }
                if self.buffer_fill == self.buffer.len() {
                    self.buffer.resize(self.buffer.len() * 2, 0); // Header line exceeds the buffer.
                }
            } else if self.finished {
                return Ok(()); // Header-only file.
            }
        }
    }

    /// Fills the local buffer to the brim. Loops because decompressing readers (`GzDecoder`)
    /// hand back one chunk at a time rather than filling the whole slice.
    fn fill_buffer(&mut self) -> Result<(), Error> {
        while !self.finished && self.buffer_fill < self.buffer.len() {
            let n_bytes = self.file.read(&mut self.buffer[self.buffer_fill..])?;
            self.buffer_fill += n_bytes;
            self.finished = n_bytes == 0;
        }
        Ok(())
    }

    /// Drops the first `n` bytes, shifting whatever is left to the front. The buffer keeps its
    /// length throughout, so this is one `memmove` of the (small) remainder — no re-zeroing of
    /// the drained region, which on a 16 MiB buffer would cost more than the parsing does.
    /// Bytes past `buffer_fill` are stale and must never be read.
    fn consume(&mut self, n: usize) {
        self.buffer.copy_within(n..self.buffer_fill, 0);
        self.buffer_fill -= n;
    }
}

/// Number of records in a cut batch — non-empty lines, matching exactly what
/// [`SamReader::next`](crate::sam_reader::SamReader::next) will yield from it, including a
/// final line that has no trailing newline.
fn count_records(batch: &[u8]) -> u64 {
    let mut records = 0;
    let mut start = 0;

    for newline in memchr_iter(b'\n', batch) {
        let mut end = newline;
        if end > start && batch[end - 1] == b'\r' {
            end -= 1;
        }
        records += (end > start) as u64;
        start = newline + 1;
    }
    records += (start < batch.len()) as u64;

    records
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAM: &[u8] = b"@HD\tVN:1.6\tSO:coordinate\n\
                         @SQ\tSN:chr1\tLN:248956422\n\
                         @SQ\tSN:chr2\tLN:242193529\n\
                         @PG\tID:bwa\tPN:bwa\n\
                         r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
                         r2\t16\tchr2\t7\t60\t4M\t*\t0\t0\tTGCA\tIIII\n";

    #[test]
    fn header_is_consumed_up_front() {
        let reader = SamByteReader::new(SAM, 1 << 16).unwrap();

        assert_eq!(reader.header().lines().count(), 4);
        assert_eq!(
            reader.header().reference_sequences(),
            vec![(b"chr1".as_slice(), 248956422), (b"chr2".as_slice(), 242193529)]
        );
    }

    /// The header alone is bigger than the buffer, so `read_header` has to drain and refill.
    #[test]
    fn header_larger_than_buffer() {
        let mut sam = Vec::new();
        for i in 0..500 {
            sam.extend_from_slice(format!("@SQ\tSN:contig{}\tLN:{}\n", i, 1000 + i).as_bytes());
        }
        sam.extend_from_slice(b"r1\t0\tcontig0\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n");

        let mut reader = SamByteReader::new(sam.as_slice(), 64).unwrap();
        assert_eq!(reader.header().reference_sequences().len(), 500);

        let mut batch = vec![0; 64];
        let bytes = reader.fill_buf(&mut batch).unwrap().unwrap();
        assert!(batch[..bytes].starts_with(b"r1\t"));
    }

    /// A record longer than the buffer forces the buffer to grow instead of cutting mid-line.
    #[test]
    fn record_larger_than_buffer() {
        let seq = vec![b'A'; 4096];
        let qual = vec![b'I'; 4096];
        let mut sam = Vec::from(b"@HD\tVN:1.6\n".as_slice());
        sam.extend_from_slice(b"long\t0\tchr1\t1\t60\t4096M\t*\t0\t0\t");
        sam.extend_from_slice(&seq);
        sam.push(b'\t');
        sam.extend_from_slice(&qual);
        sam.push(b'\n');

        let mut reader = SamByteReader::new(sam.as_slice(), 128).unwrap();
        let mut batch = vec![0; 128];
        let bytes = reader.fill_buf(&mut batch).unwrap().unwrap();

        assert_eq!(bytes, sam.len() - reader.header().len());
        assert!(reader.fill_buf(&mut batch).unwrap().is_none());
    }

    #[test]
    fn batches_end_on_record_boundaries() {
        let mut reader = SamByteReader::new(SAM, 80).unwrap();
        let mut batch = vec![0; 80];

        let mut total = 0;
        while let Some(bytes) = reader.fill_buf(&mut batch).unwrap() {
            assert_eq!(batch[bytes - 1], b'\n');
            total += count_records(&batch[..bytes]);
        }
        assert_eq!(total, 2);
    }

    #[test]
    fn counts_final_line_without_newline() {
        assert_eq!(count_records(b"a\nb\nc"), 3);
        assert_eq!(count_records(b"a\nb\nc\n"), 3);
        assert_eq!(count_records(b"a\r\nb\r\n"), 2);
        assert_eq!(count_records(b""), 0);
    }
}
