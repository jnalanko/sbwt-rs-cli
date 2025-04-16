# Changelog

## v0.3.10

- Fix `push_labels_forward` behavior at root: now fills in dollars, as documentation says.

## v0.3.9

- Make `build_last_column` public.

## v0.3.8

- Add disk usage reporting during construction.
- Make `push_labels_forward` public.

## v0.3.7

- Parallelise the subroutine `push_labels_forward`. This speeds up the k-mer dump, LCS array construction and DBG construction.
- Add a CLI command `build-lcs` to build the LCS array from the SBWT on the command line.
- Make CLI commands `lookup`, `matching-statistics`, `dump-kmers` and `dump-unitigs` all print to stdout if an output file is not given.

## v0.3.6

- Fixed a corner case when importing an index with an empty precalc table from the C++ format. This caused a crash at query time, which is fixed in this patch.


## v0.3.5

- Add `load_from_cpp_plain_matrix_format` to load the index from the C++ API
  or CLI format: https://github.com/algbio/SBWT.

## v0.3.4

- Add an optional in-memory construction implementation by Tommi Mäklin: https://github.com/jnalanko/sbwt-rs-cli/pull/4. This can be built with by enabling the feature `bpks-mem`.
- Fixes to make the crate work with WebAssembly.
- Parallel C-array construction.

## v0.3.3

- Documentation: Add an example using Needletail for FASTX parsing.
- Change usize to u32 in `IS_DNA` for platform-independence.

## v0.3.2

- Fixes to documentation.

## v0.3.1

- Added function `from_subset_seq` to SbwtIndex.

## v0.3.0

- Added serialization magic string to the start of the serialization format. This breaks compatibility with files generated using the previous versions.
- Made some functions generic over the particular implementation of SubsetSeq in preparation for the future when we will have multiple implementations.
- Added to the documentation a link to the Github repository of the CLI.

## v0.2.1

- Added CHANGELOG.md
- Updated jseqio to version 0.1.3. This makes it so that lower-case characters are made upper-case during parsing. Previously k-mers containing lower-case characters were ignored.
