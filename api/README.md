This crate provides an implementation of the Bit Matrix SBWT data structure, as described in [Small Searchable k-Spectra via Subset Rank Queries on the Spectral Burrows-Wheeler Transform](https://epubs.siam.org/doi/abs/10.1137/1.9781611977714.20), for the DNA alphabet ACGT. The data structure uses a variant of the Burrows-Wheeler transform to compress a set of k-mers in way that allows fast lookup queries. If the input k-mers are consecutive k-mers from longer underlying sequences, the index takes typically around 5 bits per distinct k-mer, supporting k-mer lookup queries at a speed of around 1 μs / k-mer on modern hardware.

Queries can be further sped up by using the Longest common suffix array, (see [here](https://link.springer.com/chapter/10.1007/978-3-031-43980-3_1)) taking roughly log(k) bits of space per k-mer. The LCS array also enables the computation of *k-bounded matching statistics* and shortest frequency-bounded suffixes (both described [here](https://www.biorxiv.org/content/10.1101/2024.02.19.580943v1)). Finally, the crate provides an interface for traversing the node-centric de Bruijn graph of the k-mers.

API documentation available at [https://docs.rs/sbwt/latest/sbwt/](https://docs.rs/sbwt/latest/sbwt/).

A CLI for key features of the API can be found at [https://github.com/jnalanko/sbwt-rs-cli](https://github.com/jnalanko/sbwt-rs-cli). The API is developed together with the CLI at [https://github.com/jnalanko/sbwt-rs-cli/tree/master/api](https://github.com/jnalanko/sbwt-rs-cli/tree/master/api).

Changes between releases are documented in the [changelog](https://github.com/jnalanko/sbwt-rs-cli/blob/master/api/CHANGELOG.md).
