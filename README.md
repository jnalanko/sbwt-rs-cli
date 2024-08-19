# sbwt-rs command-line interface

A command line interface for the [sbwt crate](https://crates.io/crates/sbwt). To compile, run `cargo build -r`. Then, run `./target/release/sbwt` for instructions. At the time of writing this README.md, the CLI supports the provides the following subdommands:

```
Command line tools using for the sbwt library.

Usage: sbwt [OPTIONS] [COMMAND]

Commands:
  build                
  benchmark            Benchmark queries on the index
  lookup               Look up the colex ranks of all k-mers. Printed as space-separated ascii integers, one line per query sequence. If a k-mer does not exist in the index, the colex rank is -1.
  matching-statistics  
  dump-kmers           Prints all k-mer strings in colex order, one per line.
  dump-unitigs         Print all unitigs in arbitrary order in FASTA format
  help                 Print this message or the help of the given subcommand(s)

Options:
  -t, --threads <threads>  Number of threads to use for parallelized operations [default: 4]
  -v, --verbose            Print more information when running.
  -h, --help               Print help
```
