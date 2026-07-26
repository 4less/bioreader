//! Generate a synthetic SAM file: `gen <path> <n_records> <read_len>`.
use std::io::{BufWriter, Write};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = &args[1];
    let n: usize = args[2].parse().unwrap();
    let read_len: usize = args[3].parse().unwrap();

    let refs = 5usize;
    let mut out = BufWriter::with_capacity(1 << 22, std::fs::File::create(path).unwrap());
    writeln!(out, "@HD\tVN:1.6\tSO:unsorted").unwrap();
    for r in 0..refs {
        writeln!(out, "@SQ\tSN:chr{}\tLN:{}", r + 1, 250_000_000 - r * 10_000_000).unwrap();
    }
    writeln!(out, "@RG\tID:sample1\tSM:sample1\tPL:ILLUMINA").unwrap();
    writeln!(out, "@PG\tID:bwa\tPN:bwa\tVN:0.7.17").unwrap();

    let bases = b"ACGT";
    let mut seq = vec![b'A'; read_len];
    let qual = vec![b'I'; read_len];

    for i in 0..n {
        for (j, b) in seq.iter_mut().enumerate() {
            *b = bases[(i + j) % 4];
        }
        // ~5% unmapped, as in a real alignment file.
        if i % 20 == 19 {
            out.write_all(b"read").unwrap();
            write!(out, "{}\t4\t*\t0\t0\t*\t*\t0\t0\t", i).unwrap();
            out.write_all(&seq).unwrap();
            out.write_all(b"\t").unwrap();
            out.write_all(&qual).unwrap();
            out.write_all(b"\n").unwrap();
            continue;
        }
        write!(
            out,
            "read{}\t{}\tchr{}\t{}\t60\t{}M\t=\t{}\t{}\t",
            i,
            if i % 2 == 0 { 99 } else { 147 },
            i % refs + 1,
            i % 200_000_000 + 1,
            read_len,
            i % 200_000_000 + 300,
            read_len * 2
        )
        .unwrap();
        out.write_all(&seq).unwrap();
        out.write_all(b"\t").unwrap();
        out.write_all(&qual).unwrap();
        write!(
            out,
            "\tNM:i:{}\tMD:Z:{}\tAS:i:{}\tXS:i:{}\tRG:Z:sample1\n",
            i % 5,
            read_len,
            read_len - i % 7,
            i % 23
        )
        .unwrap();
    }
    out.flush().unwrap();
}
