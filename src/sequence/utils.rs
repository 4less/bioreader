
static COMPLEMENT: &'static [u8] = &[
    b'T', // A -> T
    b'V', // B -> V
    b'G', // C -> G
    b'H', // D -> H
    255u8, // E 
    255u8, // F
    b'C', // G -> C
    b'D', // H -> D
    255u8, // I
    255u8, // J
    b'M', // K -> M
    255u8, // L
    b'K', // M -> K
    b'N', // N -> N
    255u8, // O
    255u8, // P
    255u8, // Q
    b'Y', // R -> Y
    b'S', // S -> S
    b'A', // T -> A
    255u8, // U
    b'B', // V -> B
    b'W', // W -> W
    255u8, // X
    b'R', // Y -> R
    b'Z', // Z -> Z
];

// A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
// T V G H - - C D - - M - K N - - - Y S A - B W - R Z

// Offset 65
// A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z
// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21  

pub fn complement(base: u8) -> u8 {
    assert!(base > 64 && base < 91);
    return COMPLEMENT[base as usize - 65]
}

pub fn reverse_complement_into_vec(from: &[u8], to: &mut Vec<u8>) {
    for &c in from.iter().rev() {
        to.push(complement(c));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complement() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'N'), b'N');
        assert_eq!(complement(b'C'), b'G');
        assert_eq!(complement(b'B'), b'V');
        assert_eq!(complement(b'G'), b'C');
        assert_eq!(complement(b'Z'), b'Z');
        assert_eq!(complement(b'D'), b'H');
        assert_eq!(complement(b'H'), b'D');
    }

    #[test]
    fn test_reverse_complement_into_vec() {
        let seq =            "CGACTGATGTCGACYBVBBZNZNNNCGACTGATC";
        let rev_complement = "GATCAGTCGNNNZNZVVBVRGTCGACATCAGTCG";

        let mut vec = Vec::new();
        reverse_complement_into_vec(seq.as_bytes(), &mut vec);
        assert_eq!(String::from_utf8_lossy(&vec), rev_complement.to_string());
    }
}