//! Record-level SAM reader.
//!
//! Mirrors [`FastqReader`](crate::fastq_reader::FastqReader) / [`FastaReader`](crate::fasta_reader::FastaReader):
//! each worker thread owns one `SamReader` with its own batch buffer, pulls a batch of complete
//! records from the shared [`SamByteReader`] (the only synchronised step) and then parses that
//! batch on its own.

use std::sync::{Arc, Mutex};

use memchr::memchr;

use crate::{
    fastq_byte_reader::FillBuffer,
    sam_byte_reader::{SamBatchBase, SamByteReader, SamReadId},
    sequence::sam_record::{RefSamRecord, SamBufferPosition},
};

#[derive(Debug, Clone)]
pub struct SamReader {
    pub buffer: Vec<u8>,
    pub buffer_fill: usize,
    buffer_pos: usize,
    buf_pos: SamBufferPosition,
    // Read-id state for the `*_ids` path only; the plain `next()` never reads these.
    // `batch_base` is the id of the current batch's first record (fetched under the reader lock
    // by `load_batch_par_ids`); `idx_in_batch` counts records emitted from this batch.
    batch_base: SamBatchBase,
    idx_in_batch: u64,
}

impl SamReader {
    pub fn new(capacity: usize) -> Self {
        Self {
            buffer: vec![0; capacity],
            buffer_fill: 0,
            buffer_pos: 0,
            buf_pos: SamBufferPosition::default(),
            batch_base: SamBatchBase::default(),
            idx_in_batch: 0,
        }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self::new(capacity)
    }

    /// Single-threaded batch load, for readers that are not shared.
    #[inline]
    /// Point the parser at a batch someone else filled.
    ///
    /// The prefetching reader hands a worker a whole buffer rather than filling the
    /// worker's own, so the parser has to be told how much of it is live. Bytes past
    /// `bytes` are stale and are never read, exactly as after `load_batch`.
    #[inline]
    pub fn set_batch(&mut self, bytes: usize) {
        self.buffer_pos = 0;
        self.buffer_fill = bytes;
    }

    pub fn load_batch(&mut self, br: &mut impl FillBuffer) -> Result<Option<()>, std::io::Error> {
        self.buffer_pos = 0;
        match br.fill_buf(&mut self.buffer)? {
            Some(bytes) if bytes > 0 => {
                self.buffer_fill = bytes;
                Ok(Some(()))
            }
            _ => {
                self.buffer_fill = 0;
                Ok(None)
            }
        }
    }

    #[inline]
    pub fn load_batch_par(
        &mut self,
        br: &mut Arc<Mutex<impl FillBuffer>>,
    ) -> Result<Option<()>, std::io::Error> {
        self.buffer_pos = 0;
        match br
            .lock()
            .expect("Locking ByteReader was unsuccessful")
            .fill_buf(&mut self.buffer)?
        {
            Some(bytes) if bytes > 0 => {
                self.buffer_fill = bytes;
                Ok(Some(()))
            }
            _ => {
                self.buffer_fill = 0;
                Ok(None)
            }
        }
    }

    /// Like [`load_batch_par`](Self::load_batch_par), but pulls the batch through
    /// [`fill_buf_ids`](SamByteReader::fill_buf_ids) and records the batch's id base, resetting
    /// the within-batch record counter. Used only by [`next_id`](Self::next_id).
    #[inline]
    pub fn load_batch_par_ids<T: std::io::Read>(
        &mut self,
        br: &Arc<Mutex<SamByteReader<T>>>,
    ) -> Result<Option<()>, std::io::Error> {
        self.buffer_pos = 0;
        match br
            .lock()
            .expect("Locking ByteReader was unsuccessful")
            .fill_buf_ids(&mut self.buffer)?
        {
            Some((bytes, base)) if bytes > 0 => {
                self.buffer_fill = bytes;
                self.batch_base = base;
                self.idx_in_batch = 0;
                Ok(Some(()))
            }
            _ => {
                self.buffer_fill = 0;
                Ok(None)
            }
        }
    }

    /// Advances to the next record in the current batch, `None` at the end of the batch.
    pub fn next(&mut self) -> Option<RefSamRecord<'_>> {
        self.parse_next()?;

        Some(RefSamRecord {
            buffer: &self.buffer,
            buf_pos: &self.buf_pos,
        })
    }

    /// The read-id counterpart of [`next`](Self::next): yields each record together with its
    /// deterministic [`SamReadId`]. Requires the batch to have been loaded through
    /// [`load_batch_par_ids`](Self::load_batch_par_ids).
    pub fn next_id(&mut self) -> Option<(RefSamRecord<'_>, SamReadId)> {
        self.parse_next()?;

        // `pos.0` is the record's start within the batch buffer, which begins at the batch's
        // file offset; adding the two gives the record's absolute byte offset in the file.
        let read_id = SamReadId {
            ordinal: self.batch_base.ordinal + self.idx_in_batch,
            byte_offset: self.batch_base.byte_offset + self.buf_pos.pos.0 as u64,
            batch_id: self.batch_base.batch_id,
            idx_in_batch: self.idx_in_batch,
        };
        self.idx_in_batch += 1;

        Some((
            RefSamRecord {
                buffer: &self.buffer,
                buf_pos: &self.buf_pos,
            },
            read_id,
        ))
    }

    /// Locates the next record line and fills `buf_pos` with its field offsets. Empty lines are
    /// skipped (they carry no record, and the byte reader's record count agrees).
    fn parse_next(&mut self) -> Option<()> {
        loop {
            if self.buffer_pos >= self.buffer_fill {
                return None;
            }

            let start = self.buffer_pos;
            let line_end = match memchr(b'\n', &self.buffer[start..self.buffer_fill]) {
                Some(newline) => start + newline,
                None => self.buffer_fill, // Final line of the file, no trailing newline.
            };
            self.buffer_pos = line_end + 1;

            let mut end = line_end;
            if end > start && self.buffer[end - 1] == b'\r' {
                end -= 1;
            }
            if end == start {
                continue;
            }

            if !self.buf_pos.parse(&self.buffer, start, end) {
                panic!(
                    "Invalid SAM record, expected 11 tab separated fields: {}",
                    String::from_utf8_lossy(&self.buffer[start..end])
                );
            }
            return Some(());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAM: &[u8] = b"@HD\tVN:1.6\n\
                         @SQ\tSN:chr1\tLN:1000\n\
                         r1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\n\
                         r2\t16\tchr1\t20\t30\t4M\t*\t0\t0\tTTTT\t!!!!\n\
                         r3\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";

    fn collect(sam: &[u8], buffer_size: usize) -> Vec<(Vec<u8>, u64)> {
        let mut byte_reader = SamByteReader::new(sam, buffer_size).unwrap();
        let mut reader = SamReader::with_capacity(buffer_size);
        let mut records = Vec::new();

        while let Some(()) = reader.load_batch(&mut byte_reader).unwrap() {
            while let Some(record) = reader.next() {
                assert!(record.valid_extended(), "invalid record {}", record);
                records.push((record.qname().to_vec(), record.pos()));
            }
        }
        records
    }

    #[test]
    fn reads_every_record() {
        let records = collect(SAM, 1 << 16);
        assert_eq!(
            records,
            vec![
                (b"r1".to_vec(), 10),
                (b"r2".to_vec(), 20),
                (b"r3".to_vec(), 0),
            ]
        );
    }

    #[test]
    fn small_buffers_do_not_change_the_records() {
        // Several batches instead of one; the record stream must be identical.
        for buffer_size in [32, 48, 64, 100, 4096] {
            assert_eq!(collect(SAM, buffer_size), collect(SAM, 1 << 16));
        }
    }

    #[test]
    fn handles_missing_trailing_newline() {
        let sam = &SAM[..SAM.len() - 1];
        assert_eq!(collect(sam, 1 << 16).len(), 3);
        assert_eq!(collect(sam, 40).len(), 3);
    }

    #[test]
    fn handles_windows_line_endings() {
        let sam: Vec<u8> = String::from_utf8_lossy(SAM).replace('\n', "\r\n").into_bytes();
        let records = collect(&sam, 1 << 16);

        assert_eq!(records.len(), 3);
        assert_eq!(records[1].0, b"r2".to_vec());
    }
}
