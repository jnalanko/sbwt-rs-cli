# SBWT command-line interface

Command-line interface for building, querying, and performing set operations on [SBWT](https://crates.io/crates/sbwt) (*k*-mer Spectral Burrows-Wheeler Transform) indexes.

# Installation

```bash
cargo build --release
```

The program will be compiled to `./target/release/sbwt`. Run the program without any command-line arguments for instructions on how to use it.

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

# Citation

If you use this code, please cite the following paper:

```
@inproceedings{alanko2023small,
  title={Small searchable k-spectra via subset rank queries on the spectral Burrows-Wheeler transform},
  author={Alanko, Jarno N and Puglisi, Simon J and Vuohtoniemi, Jaakko},
  booktitle={SIAM Conference on Applied and Computational Discrete Algorithms (ACDA23)},
  pages={225--236},
  year={2023},
  organization={SIAM}
}
```


