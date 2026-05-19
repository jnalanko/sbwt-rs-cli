# Set Operations on SBWT Indexes

Command-line interface for building, querying, and performing set operations on [SBWT](https://crates.io/crates/sbwt) (*k*-mer Spectral Burrows-Wheeler Transform) indexes. This is an extended version of the original sbwt-rs CLI with added support for **union/merge**, **intersection**, and **set difference** directly on SBWT indexes, producing a valid compressed SBWT index for the result *k*-mer set.

## Quick start

```bash
cargo build --release
# then use ./target/release/sbwt --help

# Build two indexes
sbwt build input1.fna -o index1
sbwt build input2.fna -o index2

# Set operations
sbwt merge     index1.sbwt index2.sbwt -o union        
sbwt intersect index1.sbwt index2.sbwt -o intersection 
sbwt difference index1.sbwt index2.sbwt -o difference  
```

Add `--low-ram` to any set-operation command to trade some runtime for lower peak memory usage.

The binary is placed at `./target/release/sbwt`. Run `sbwt --help` or `sbwt <subcommand> --help` for full usage.

## Subcommands

| Subcommand | Description |
|---|---|
| `build` | Build an SBWT index from FASTA/FASTQ input |
| `build-lcs` | Build the LCS array for an existing SBWT index |
| `merge` | Compute the union of two SBWT indexes |
| `intersect` | Compute the intersection of two SBWT indexes |
| `difference` | Compute the set difference (index1 \ index2) of two SBWT indexes |
| `jaccard` | Compute the Jaccard index of the *k*-mer sets in two SBWT indexes |
| `lookup` | Look up the colex ranks of query *k*-mers |
| `matching-statistics` | Compute matching statistics against an SBWT index |
| `dump-kmers` | Print all *k*-mer strings in colexicographic order |
| `dump-unitigs` | Print all unitigs in FASTA format |
| `kmer-at-colex` | Print the *k*-mer at a given colexicographic rank |
| `stats` | Print summary statistics of an SBWT index |
| `check` | Verify the structural integrity of an SBWT index |
| `benchmark` | Benchmark query performance on an index |

### Set operation flags

All three set-operation subcommands (`merge`, `intersect`, `difference`) accept:

- `-t / --threads <N>` — number of parallel threads (default: 4)
- `--low-ram` — compact 3-bit character encoding during merge plan construction instead of one byte per character; reduces peak RAM at the cost of some runtime

## Large-scale construction via the binomial merge pipeline

The `binomial_merge_pipeline/` directory contains the script used to construct a 661K bacterial SBWT index from the [ENA2018-bacteria-661k](https://ftp.ebi.ac.uk/pub/databases/ENA2018-bacteria-661k) dataset. The pipeline builds batch indexes and merges them incrementally using a binomial merge stack.

### Before running `build_661k.sh`

Open `binomial_merge_pipeline/build_661k.sh` and set the following variables at the top of the **Configuration** section:

| Variable | What to set |
|---|---|
| `MANIFEST` | Path to `sampleid_assembly_paths.txt` (see download command in the script header) |
| `SBWT` | Path to the compiled `sbwt` binary, e.g. `./target/release/sbwt` |
| `WORKDIR` | Directory where all output, intermediate indexes, and logs will be written |
| `THREADS` | Number of CPU threads to use for building and merging |
| `MEM_GB` | Memory budget in GB passed to `sbwt build --mem-gb` |
| `BUILD_MEM_LIMIT` | systemd cgroup memory cap for build steps — set higher than `MEM_GB` |
| `MERGE_MEM_LIMIT` | systemd cgroup memory cap for merge steps |

`DOWNLOAD_PARALLELISM` controls how many FASTA files are downloaded simultaneously per batch; the default of 8 is a safe value for the EBI FTP server.

## Credits

This repository extends [jnalanko/sbwt-rs-cli](https://github.com/jnalanko/sbwt-rs-cli) with set-operation support. The set-operation algorithms were implemented in collaboration with [Jarno Alanko](https://github.com/jnalanko), under the supervision of Camille Marchet and Simon Puglisi.

The underlying SBWT data structure is described in:

> Jarno N. Alanko, Simon J. Puglisi, and Jaakko Vuohtoniemi.
> "Small Searchable κ-Spectra via Subset Rank Queries on the Spectral Burrows-Wheeler Transform."
> *SIAM Conference on Applied and Computational Discrete Algorithms (ACDA23)*, pp. 225–236.
> DOI: [10.1137/1.9781611977714.20](https://doi.org/10.1137/1.9781611977714.20)
