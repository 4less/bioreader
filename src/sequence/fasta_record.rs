/// A FASTQ record that borrows data from a buffer
#[derive(Debug, Clone)]


pub struct OwnedFastaRecord {
    pub header: Vec<u8>,
    pub sequence: Vec<u8>,
}

impl OwnedFastaRecord { // Record for 
    pub fn new() -> Self {
        Self {
            header: Vec::new(),
            sequence: Vec::new()
        }
    }

    pub fn clear(&mut self) {
        self.header.clear();
        self.sequence.clear();
    }

    #[inline]
    pub fn head(&self) -> &[u8] {
        &self.header
    }

    #[inline]
    pub fn seq(&self) -> &[u8] {
        &self.sequence
    }

    #[inline]
    pub fn valid(&self) -> bool {
        let valid = self.seq().iter().all(|&c| c == b'A' || c == b'C' || c == b'G' || c == b'T' || c == b'N');
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
        let valid = self.seq().iter().all(|&c| c == b'A' || c == b'C' || c == b'G' || c == b'T');
        return valid
    }

    pub fn to_string(&self) -> String {
        let header = String::from_utf8(self.header.clone()).expect("Invalid UTF-8");
        let sequence = String::from_utf8(self.sequence.clone()).expect("Invalid UTF-8");
        format!("{}\n{}", header, sequence)
    }
}