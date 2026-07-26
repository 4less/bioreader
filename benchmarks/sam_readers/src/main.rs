//! Head-to-head read benchmark: bioreader's SAM reader vs. noodles-sam vs. a naive
//! `BufReader` + `split('\t')` baseline.
//!
//! Every contender does the *same* per-record work — sum SEQ length + POS + the `NM` tag —
//! and returns (checksum, record count) so the runs can be verified to agree.
//!
//! Run: `cargo run --release -- <file.sam> [runs]`

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use bioreader::sam_byte_reader::SamByteReader;
use bioreader::sam_reader::SamReader;
use bioreader::sequence::sam_record::RefSamRecord;

use noodles_sam::{self as sam, alignment::record::data::field::Tag};

static BUFFER_SIZE_CELL: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(1 << 24);

/// Batch size for the chunking readers; overridable so it can be swept.
#[allow(non_snake_case)]
#[inline]
fn BUFFER_SIZE() -> usize {
    BUFFER_SIZE_CELL.load(std::sync::atomic::Ordering::Relaxed)
}

static BUFREAD_SIZE_CELL: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(1 << 20);

/// `BufReader` capacity for the line-at-a-time readers; swept the same way.
#[allow(non_snake_case)]
#[inline]
fn BUFREAD_SIZE() -> usize {
    BUFREAD_SIZE_CELL.load(std::sync::atomic::Ordering::Relaxed)
}

type Outcome = (u64, u64); // (checksum, records)

// ---------------------------------------------------------------- bioreader

#[inline]
fn work_bioreader(record: &RefSamRecord) -> u64 {
    let nm = record.tag(b"NM").and_then(|tag| tag.as_int()).unwrap_or(0) as u64;
    record.seq().len() as u64 + record.pos() + nm
}

/// Single threaded: byte reader and record reader without any lock or thread.
fn bioreader_single(path: &str) -> Outcome {
    let mut byte_reader = SamByteReader::new(File::open(path).unwrap(), BUFFER_SIZE()).unwrap();
    let mut reader = SamReader::with_capacity(BUFFER_SIZE());

    let mut checksum = 0;
    let mut records = 0;
    while let Some(()) = reader.load_batch(&mut byte_reader).unwrap() {
        while let Some(record) = reader.next() {
            checksum += work_bioreader(&record);
            records += 1;
        }
    }
    (checksum, records)
}

#[derive(Default, Clone)]
struct Sum {
    checksum: u64,
    records: u64,
}

impl bioreader::parallel::fastq::Merge for Sum {
    fn merge_from(&mut self, other: &mut Self) {
        self.checksum += other.checksum;
        self.records += other.records;
    }
}

fn bioreader_par(path: &str, threads: u32) -> Outcome {
    let sum = bioreader::parallel::sam::read_sam_state_par(
        File::open(path).unwrap(),
        BUFFER_SIZE(),
        threads,
        Sum::default(),
        |record: &RefSamRecord, state: &mut Sum| {
            state.checksum += work_bioreader(record);
            state.records += 1;
        },
    )
    .unwrap();
    (sum.checksum, sum.records)
}

/// The read-id variant, to price the id plumbing.
fn bioreader_par_ids(path: &str, threads: u32) -> Outcome {
    let sum = bioreader::parallel::sam::read_sam_state_par_ids(
        File::open(path).unwrap(),
        BUFFER_SIZE(),
        threads,
        Sum::default(),
        |record: &RefSamRecord, id: bioreader::sam_byte_reader::SamReadId, state: &mut Sum| {
            state.checksum += work_bioreader(record) + (id.ordinal & 0);
            state.records += 1;
        },
    )
    .unwrap();
    (sum.checksum, sum.records)
}

// ------------------------------------------------------------------ noodles

#[inline]
fn work_noodles(record: &sam::Record) -> u64 {
    let pos = match record.alignment_start() {
        Some(start) => usize::from(start.unwrap()) as u64,
        None => 0,
    };
    let nm = match record.data().get(&Tag::EDIT_DISTANCE) {
        Some(value) => value.unwrap().as_int().unwrap_or(0) as u64,
        None => 0,
    };
    record.sequence().len() as u64 + pos + nm
}

fn noodles_single(path: &str) -> Outcome {
    let file = BufReader::with_capacity(BUFREAD_SIZE(), File::open(path).unwrap());
    let mut reader = sam::io::Reader::new(file);
    reader.read_header().unwrap();

    let mut record = sam::Record::default();
    let mut checksum = 0;
    let mut records = 0;
    while reader.read_record(&mut record).unwrap() != 0 {
        checksum += work_noodles(&record);
        records += 1;
    }
    (checksum, records)
}

/// noodles' record parser fed by bioreader's parallel chunker: each worker runs a noodles
/// reader over an in-memory batch. Isolates parser cost from I/O and gives noodles — which
/// ships no parallel SAM reader — a multithreaded showing.
fn noodles_par(path: &str, threads: u32) -> Outcome {
    let byte_reader = Arc::new(Mutex::new(
        SamByteReader::new(File::open(path).unwrap(), BUFFER_SIZE()).unwrap(),
    ));

    let mut checksum = 0;
    let mut records = 0;
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for _ in 0..threads {
            let mut reader = SamReader::with_capacity(BUFFER_SIZE());
            let mut byte_reader = byte_reader.clone();

            handles.push(scope.spawn(move || {
                let mut record = sam::Record::default();
                let mut checksum = 0u64;
                let mut records = 0u64;

                while let Some(()) = reader.load_batch_par(&mut byte_reader).unwrap() {
                    let batch = &reader.buffer[..reader.buffer_fill];
                    let mut inner = sam::io::Reader::new(batch);
                    while inner.read_record(&mut record).unwrap() != 0 {
                        checksum += work_noodles(&record);
                        records += 1;
                    }
                }
                (checksum, records)
            }));
        }
        for handle in handles {
            let (c, r) = handle.join().unwrap();
            checksum += c;
            records += r;
        }
    });
    (checksum, records)
}

// -------------------------------------------------------------------- naive

/// What a straightforward implementation looks like: read a line, split on tabs, parse.
fn naive_single(path: &str) -> Outcome {
    let mut reader = BufReader::with_capacity(BUFREAD_SIZE(), File::open(path).unwrap());

    let mut line = Vec::new();
    let mut checksum = 0;
    let mut records = 0;
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line).unwrap() == 0 {
            break;
        }
        let trimmed = line.strip_suffix(b"\n").unwrap_or(&line);
        if trimmed.is_empty() || trimmed[0] == b'@' {
            continue;
        }

        let fields: Vec<&[u8]> = trimmed.split(|&c| c == b'\t').collect();
        let pos: u64 = std::str::from_utf8(fields[3]).unwrap().parse().unwrap();
        let nm = fields[11..]
            .iter()
            .find(|field| field.starts_with(b"NM:i:"))
            .map(|field| std::str::from_utf8(&field[5..]).unwrap().parse::<u64>().unwrap())
            .unwrap_or(0);

        checksum += fields[9].len() as u64 + pos + nm;
        records += 1;
    }
    (checksum, records)
}

// ------------------------------------------------------- floors and gzip

/// The I/O floor: pull the file through `read()` and touch every byte, no parsing at all.
fn raw_read(path: &str) -> Outcome {
    let mut file = File::open(path).unwrap();
    let mut buffer = vec![0; BUFFER_SIZE()];
    let mut checksum = 0u64;
    loop {
        let bytes = file.read(&mut buffer).unwrap();
        if bytes == 0 {
            break;
        }
        checksum += memchr::memchr_iter(b'\n', &buffer[..bytes]).count() as u64;
    }
    (checksum, checksum)
}

fn bioreader_gz(path: &str, threads: u32) -> Outcome {
    let sum = bioreader::parallel::sam::read_sam_state_par(
        flate2::read::MultiGzDecoder::new(File::open(path).unwrap()),
        BUFFER_SIZE(),
        threads,
        Sum::default(),
        |record: &RefSamRecord, state: &mut Sum| {
            state.checksum += work_bioreader(record);
            state.records += 1;
        },
    )
    .unwrap();
    (sum.checksum, sum.records)
}

fn noodles_gz(path: &str) -> Outcome {
    let file = BufReader::with_capacity(
        BUFREAD_SIZE(),
        flate2::read::MultiGzDecoder::new(File::open(path).unwrap()),
    );
    let mut reader = sam::io::Reader::new(file);
    reader.read_header().unwrap();

    let mut record = sam::Record::default();
    let mut checksum = 0;
    let mut records = 0;
    while reader.read_record(&mut record).unwrap() != 0 {
        checksum += work_noodles(&record);
        records += 1;
    }
    (checksum, records)
}

// ----------------------------------------------------------------- harness

fn warm_cache(path: &str) {
    let mut file = File::open(path).unwrap();
    let mut buffer = vec![0; 1 << 22];
    while file.read(&mut buffer).unwrap() != 0 {}
}

/// Median wall time of `runs` runs, plus the throughput it works out to.
fn measure(runs: usize, size: u64, f: impl Fn() -> Outcome) -> (Outcome, Duration, f64, f64) {
    let mut timings = Vec::new();
    let mut outcome = (0, 0);
    for _ in 0..runs {
        let start = Instant::now();
        outcome = f();
        timings.push(start.elapsed());
    }
    timings.sort();
    let median: Duration = timings[timings.len() / 2];

    let seconds = median.as_secs_f64();
    let gb_per_s = size as f64 / seconds / 1e9;
    let mrec_per_s = outcome.1 as f64 / seconds / 1e6;
    (outcome, median, gb_per_s, mrec_per_s)
}

fn bench(name: &str, runs: usize, size: u64, f: impl Fn() -> Outcome) -> Outcome {
    let (outcome, median, gb_per_s, mrec_per_s) = measure(runs, size, f);
    println!(
        "{:<28} {:>8.0} ms   {:>5.2} GB/s   {:>6.2} M rec/s   [{} records, checksum {}]",
        name,
        median.as_secs_f64() * 1000.0,
        gb_per_s,
        mrec_per_s,
        outcome.1,
        outcome.0
    );
    outcome
}

// ---------------------------------------------------------- regression floor

#[derive(Clone, Copy, PartialEq)]
enum Unit {
    MRecPerSec,
    GbPerSec,
    /// Speedup over noodles-sam on the same file — the portable signal, since it cancels out
    /// most of the hardware difference.
    RelativeToNoodles,
}

impl Unit {
    fn label(&self) -> &'static str {
        match self {
            Unit::MRecPerSec => "M rec/s",
            Unit::GbPerSec => "GB/s",
            Unit::RelativeToNoodles => "x noodles",
        }
    }
}

/// A floor is "this must not get slower". Values are 85% of the medians measured on the
/// baseline machine (see `README.md`), leaving headroom for run-to-run noise while still
/// catching a real regression. Re-baseline deliberately, never to make a red run go green.
struct Floor {
    scenario: &'static str,
    unit: Unit,
    floor: f64,
    baseline: f64,
}

/// Baseline: AMD Ryzen AI 9 HX 370 (12C/24T), 93 GB RAM, rustc nightly, 2026-07-26.
/// Data: `gen short.sam 4000000 150` (1.58 GB, 4M records).
const SHORT_FLOORS: &[Floor] = &[
    Floor { scenario: "bioreader 1 thread", unit: Unit::MRecPerSec, floor: 9.7, baseline: 11.45 },
    Floor { scenario: "bioreader 4 threads", unit: Unit::MRecPerSec, floor: 21.0, baseline: 24.96 },
    Floor { scenario: "bioreader 8 threads", unit: Unit::MRecPerSec, floor: 23.4, baseline: 27.61 },
    Floor { scenario: "bioreader+ids 8 threads", unit: Unit::MRecPerSec, floor: 15.1, baseline: 17.86 },
    Floor { scenario: "bioreader 1 thread", unit: Unit::RelativeToNoodles, floor: 1.45, baseline: 1.73 },
];

/// Data: `gen long.sam 60000 10000` (1.21 GB, 60k records of 10 kb).
const LONG_FLOORS: &[Floor] = &[
    Floor { scenario: "bioreader 1 thread", unit: Unit::GbPerSec, floor: 13.6, baseline: 16.09 },
    Floor { scenario: "bioreader 8 threads", unit: Unit::GbPerSec, floor: 8.4, baseline: 9.94 },
    Floor { scenario: "bioreader 1 thread", unit: Unit::RelativeToNoodles, floor: 0.74, baseline: 0.87 },
];

/// Runs the floor scenarios and reports PASS/FAIL per row. Returns false if any floor is missed.
fn check_floors(path: &str, runs: usize, size: u64, floors: &[Floor]) -> bool {
    // Both readers at the settings the report uses.
    BUFFER_SIZE_CELL.store(1 << 20, std::sync::atomic::Ordering::Relaxed);
    BUFREAD_SIZE_CELL.store(1 << 18, std::sync::atomic::Ordering::Relaxed);

    let noodles = measure(runs, size, || noodles_single(path));
    let mut passed = true;

    for floor in floors {
        let (_, _, gb_per_s, mrec_per_s) = match floor.scenario {
            "bioreader 1 thread" => measure(runs, size, || bioreader_single(path)),
            "bioreader 4 threads" => measure(runs, size, || bioreader_par(path, 4)),
            "bioreader 8 threads" => measure(runs, size, || bioreader_par(path, 8)),
            "bioreader+ids 8 threads" => measure(runs, size, || bioreader_par_ids(path, 8)),
            other => panic!("unknown floor scenario: {}", other),
        };

        let value = match floor.unit {
            Unit::MRecPerSec => mrec_per_s,
            Unit::GbPerSec => gb_per_s,
            Unit::RelativeToNoodles => gb_per_s / noodles.2,
        };
        let ok = value >= floor.floor;
        passed &= ok;

        println!(
            "{}  {:<26} {:>7.2} {:<10} floor {:>6.2}  (baseline {:.2})",
            if ok { "PASS" } else { "FAIL" },
            floor.scenario,
            value,
            floor.unit.label(),
            floor.floor,
            floor.baseline
        );
    }
    passed
}

fn main() {
    let args: Vec<String> = std::env::args().collect();

    // `--check short|long <file> [runs]` runs only the regression floor.
    if args.get(1).map(String::as_str) == Some("--check") {
        let profile = args[2].as_str();
        let path = &args[3];
        let runs: usize = args.get(4).map(|r| r.parse().unwrap()).unwrap_or(3);
        let size = std::fs::metadata(path).unwrap().len();

        println!(
            "\nregression floor: {} profile, {} ({:.2} GB), {} runs, median\n",
            profile,
            path,
            size as f64 / 1e9,
            runs
        );
        warm_cache(path);

        let floors = match profile {
            "short" => SHORT_FLOORS,
            "long" => LONG_FLOORS,
            other => panic!("unknown profile {}, expected `short` or `long`", other),
        };
        let passed = check_floors(path, runs, size, floors);
        println!();
        std::process::exit(if passed { 0 } else { 1 });
    }

    let path = &args[1];
    let runs: usize = args.get(2).map(|r| r.parse().unwrap()).unwrap_or(3);
    let size = std::fs::metadata(path).unwrap().len();

    println!(
        "\nfile: {} ({:.2} GB), {} runs each, median reported, page cache warm\n",
        path,
        size as f64 / 1e9,
        runs
    );
    warm_cache(path);

    let mut expected: Option<Outcome> = None;
    let mut check = |name: &str, outcome: Outcome| match expected {
        None => expected = Some(outcome),
        Some(first) => assert_eq!(first, outcome, "{} disagrees with the first run", name),
    };

    // Batch size sweep: the chunking readers copy every byte through a buffer of this size, so
    // whether that buffer stays in cache dominates their throughput.
    println!("-- batch size sweep --");
    for shift in [16u32, 18, 20, 22, 24, 26] {
        BUFFER_SIZE_CELL.store(1 << shift, std::sync::atomic::Ordering::Relaxed);
        let label = format!("{} KiB", (1usize << shift) / 1024);
        bench(&format!("  floor       {:>9}", label), runs, size, || raw_read(path));
        check(
            "sweep 1t",
            bench(&format!("  bioreader 1t{:>9}", label), runs, size, || bioreader_single(path)),
        );
        check(
            "sweep 4t",
            bench(&format!("  bioreader 4t{:>9}", label), runs, size, || bioreader_par(path, 4)),
        );
    }

    println!();
    for shift in [16u32, 18, 20, 22] {
        BUFREAD_SIZE_CELL.store(1 << shift, std::sync::atomic::Ordering::Relaxed);
        let label = format!("{} KiB", (1usize << shift) / 1024);
        check(
            "sweep noodles",
            bench(&format!("  noodles-sam 1t{:>8}", label), runs, size, || noodles_single(path)),
        );
    }
    BUFREAD_SIZE_CELL.store(1 << 18, std::sync::atomic::Ordering::Relaxed);
    BUFFER_SIZE_CELL.store(1 << 20, std::sync::atomic::Ordering::Relaxed);
    println!("\n-- comparison at best settings (bioreader 1 MiB batches, noodles 256 KiB BufReader) --");

    bench("read+memchr (I/O floor)", runs, size, || raw_read(path));
    check("bioreader 1t", bench("bioreader 1 thread", runs, size, || bioreader_single(path)));
    check("noodles 1t", bench("noodles-sam 1 thread", runs, size, || noodles_single(path)));
    check("naive 1t", bench("naive split 1 thread", runs, size, || naive_single(path)));

    println!();
    for threads in [1u32, 2, 4, 8, 16] {
        let name = format!("bioreader {} threads", threads);
        check(&name, bench(&name, runs, size, || bioreader_par(path, threads)));
    }

    println!();
    for threads in [1u32, 2, 4, 8, 16] {
        let name = format!("noodles-sam {} threads*", threads);
        check(&name, bench(&name, runs, size, || noodles_par(path, threads)));
    }

    println!();
    for threads in [1u32, 8] {
        let name = format!("bioreader+ids {} threads", threads);
        check(&name, bench(&name, runs, size, || bioreader_par_ids(path, threads)));
    }

    // Gzipped input, if a `<path>.gz` sits next to the plain file. Uncompressed size is used
    // for the throughput columns so the rows stay comparable with the ones above.
    let gz_path = format!("{}.gz", path);
    if std::fs::metadata(&gz_path).is_ok() {
        println!();
        for threads in [1u32, 4, 8] {
            let name = format!("bioreader gz {} threads", threads);
            check(&name, bench(&name, runs, size, || bioreader_gz(&gz_path, threads)));
        }
        check("noodles gz", bench("noodles-sam gz 1 thread", runs, size, || noodles_gz(&gz_path)));
    }

    println!("\n* noodles has no parallel SAM reader; its parser is fed by bioreader's chunker.\n");
}
