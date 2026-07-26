# SAM reader benchmark

Head-to-head read throughput for bioreader's SAM reader against
[noodles-sam](https://crates.io/crates/noodles-sam) and a naive `BufReader` + `split('\t')`
baseline, plus a regression floor to run after touching the reader.

Every contender does the *same* per-record work — sum `SEQ` length + `POS` + the `NM` tag — and
returns `(checksum, record count)`. The harness asserts all implementations produce the same
checksum, so a "fast" run that skips work fails instead of scoring well.

## Generating data

```sh
cargo build --release
./target/release/gen short.sam 4000000 150     # 1.58 GB, 4M Illumina-shaped records
./target/release/gen long.sam    60000 10000   # 1.21 GB, 60k long reads
pigz -1 -k short.sam                           # optional: adds the gzip rows
```

Both files carry a 8-line header, ~5% unmapped records and five optional fields
(`NM MD AS XS RG`) per mapped record.

## Running

```sh
./target/release/sambench short.sam 5     # full report: batch-size sweep + comparison
./target/release/sambench --check short short.sam 3    # regression floor, exits 1 on failure
./target/release/sambench --check long  long.sam  3
```

Everything is measured with the page cache warm (the harness reads the file once first), so the
numbers are parser throughput, not disk throughput.

## Regression floor

`--check` compares against `SHORT_FLOORS` / `LONG_FLOORS` in `src/main.rs`. The floors are 85%
of the medians measured on the baseline machine, which leaves room for run-to-run noise while
still catching a real regression. Each floor also records the baseline value it came from.

Baseline: AMD Ryzen AI 9 HX 370 (12C/24T), 93 GB RAM, Linux 7.0.11, rustc nightly, 2026-07-26.

| profile | scenario | baseline | floor |
| --- | --- | --- | --- |
| short | bioreader 1 thread | 11.45 M rec/s | 9.7 |
| short | bioreader 4 threads | 24.96 M rec/s | 21.0 |
| short | bioreader 8 threads | 27.61 M rec/s | 23.4 |
| short | bioreader+ids 8 threads | 17.86 M rec/s | 15.1 |
| short | bioreader vs noodles, 1 thread | 1.73x | 1.45x |
| long | bioreader 1 thread | 16.09 GB/s | 13.6 |
| long | bioreader 8 threads | 9.94 GB/s | 8.4 |
| long | bioreader vs noodles, 1 thread | 0.87x | 0.74x |

The `x noodles` rows are the portable signal: absolute throughput moves with the machine, the
ratio against another reader on the same machine mostly does not. On different hardware, expect
the absolute floors to need re-baselining and the ratios to hold.

## Buffer size matters more than anything else

Both readers copy every byte through a buffer, so whether that buffer stays in cache dominates.
The sweep at the top of the report exists to keep that honest — 16 MiB batches cost bioreader
~35% on short reads and ~2.3x on long reads compared with 1 MiB. `1 << 20` is the setting the
comparison and the floors use.
