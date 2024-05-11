# BioReader

Bioreader is a rust project for reading bioinformatics files but currently only supports reading fastq file. Its design philosophy is to parallelize reading files by moving all sequence based checks into the worker threads and keep the IO thread clean. Bioreader is currently faster than Noodles and needletail, two common rust libraries for reading fastq files. Bioreader is still in its early development stage and will only be released as a crate once enough unit and integration tests are in place. Benchmarks against noodles and needletail will follow.

# Roadmap
- Fasta reading
- Fastq writing
- Fasta writing
- Sam reading
- Sam writing
