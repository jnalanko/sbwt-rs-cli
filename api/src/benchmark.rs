//! Benchmarking queries on an existing [SbwtIndex].

use rand::{rngs::StdRng, Rng};
use rand::SeedableRng;

use crate::streaming_index::StreamingIndex;
use crate::{sbwt::SbwtIndex, streaming_index::LcsArray, subsetseq::SubsetSeq};

fn sample_kmers<SS: SubsetSeq>(sbwt: &SbwtIndex<SS>, n_queries: usize) -> Vec<Vec<u8>> {
    let mut rng = StdRng::from_entropy();
    let mut queries = Vec::<Vec::<u8>>::new();
    for _ in 0..n_queries {
        loop {
            let colex = rng.gen_range(0, sbwt.n_kmers());
            let kmer = sbwt.access_kmer(colex);
            if kmer.iter().all(|&c| c != b'$') {
                // Not a dummy k-mer
                queries.push(kmer);
                break;
            }
        }
    }
    queries
}

fn generate_random_kmers(n_queries: usize, k: usize, rng: &mut StdRng) -> Vec<Vec<u8>> {
    let mut queries = Vec::<Vec::<u8>>::new();
    for _ in 0..n_queries {
        let mut kmer = Vec::<u8>::with_capacity(k);
        for _ in 0..k {
            let c = match rng.gen_range(0, 4) {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => unreachable!(),
            };
            kmer.push(c);
        }
        queries.push(kmer);
    }
    queries
}


fn benchmark_single_positive_query<SS: SubsetSeq>(sbwt: &SbwtIndex<SS>) {

    let n_queries = 1_000_000;

    log::info!("Sampling {} k-mers", n_queries);
    let queries = sample_kmers(sbwt, n_queries);

    log::info!("Querying {} positive k-mers", n_queries);
    let start_time = std::time::Instant::now();
    let mut checksum = 0_usize;
    for kmer in queries.iter(){
        checksum += sbwt.search(kmer).unwrap().start;
    }
    let end_time = std::time::Instant::now();
    log::info!("Sum of answers: {}", checksum);
    println!("Elapsed time: {:.2} seconds", (end_time - start_time).as_secs_f64());
    println!("{:.2} nanoseconds / k-mer", (end_time - start_time).as_nanos() as f64 / queries.len() as f64);
}

fn benchmark_streaming_random_query<SS: SubsetSeq>(sbwt: &SbwtIndex<SS>, lcs: &LcsArray) {
    let n_queries = 100;
    let query_len = 10_000;

    let mut rng = StdRng::from_entropy();
    log::info!("Generating {} random sequences of length {}", n_queries, query_len);
    let queries = generate_random_kmers(n_queries, query_len, &mut rng);

    let index = StreamingIndex::new(sbwt, lcs);

    log::info!("Querying {} random sequences of length {}", n_queries, query_len);
    let start_time = std::time::Instant::now();
    let mut checksum = 0_usize;
    let mut n_kmers_queried = 0_usize;
    for query in queries.iter(){
        for x in index.matching_statistics(query).iter(){
            checksum += x.0;
        }
        n_kmers_queried += std::cmp::max(query.len() as isize - sbwt.k() as isize + 1, 0) as usize;
    }
    let end_time = std::time::Instant::now();
    log::info!("Sum of answers: {}", checksum);
    println!("Elapsed time: {:.2} seconds", (end_time - start_time).as_secs_f64());
    println!("{:.2} nanoseconds / k-mer", (end_time - start_time).as_nanos() as f64 / n_kmers_queried as f64);
}

fn benchmark_single_random_query<SS: SubsetSeq>(sbwt: &SbwtIndex<SS>) {
    
    let n_queries = 1_000_000;
    let mut rng = StdRng::from_entropy();

    log::info!("Generating {} random k-mers", n_queries);
    let queries = generate_random_kmers(n_queries, sbwt.k(), &mut rng);

    log::info!("Querying {} random k-mers", n_queries);
    let start_time = std::time::Instant::now();
    let mut checksum = 0_usize;
    for kmer in queries.iter(){
        checksum += sbwt.search(kmer).is_some() as usize;
    }
    let end_time = std::time::Instant::now();
    log::info!("Sum of answers: {}", checksum);
    println!("Elapsed time: {:.2} seconds", (end_time - start_time).as_secs_f64());
    println!("{:.2} nanoseconds / k-mer", (end_time - start_time).as_nanos() as f64 / queries.len() as f64);

}

// Assumes sbwt has select support
fn benchmark_access<SS: SubsetSeq>(sbwt: &SbwtIndex<SS>){
    let n_queries = 1_000_000;
    let mut rng = StdRng::from_entropy();
    let queries = (0..n_queries).map(|_| rng.gen_range(0, sbwt.n_kmers())).collect::<Vec<usize>>();

    log::info!("Accessing {} k-mers in random order", n_queries);
    let mut buf = Vec::<u8>::with_capacity(sbwt.k());
    let start_time = std::time::Instant::now();
    let mut checksum = 0_usize;
    for colex in queries.iter(){
        buf.clear();
        sbwt.push_kmer_to_vec(*colex, &mut buf);
        checksum += buf[0] as usize;
    }
    let end_time = std::time::Instant::now();
    log::info!("Sum of answers: {}", checksum);
    println!("Elapsed time: {:.2} seconds", (end_time - start_time).as_secs_f64());
    println!("{:.2} nanoseconds / k-mer", (end_time - start_time).as_nanos() as f64 / queries.len() as f64);

}

/// Bechmark various queries on the sbwt. If the LCS array is given, also benchmarks
/// streaming search.
pub fn benchmark_all<SS: SubsetSeq>(mut sbwt: SbwtIndex<SS>, lcs: Option<LcsArray>) {

    log::info!("Preparing for benchmarks");
    sbwt.build_select();
    benchmark_single_positive_query(&sbwt);
    benchmark_single_random_query(&sbwt);

    if let Some(lcs) = lcs {
        // Todo: streaming positive. Where do we get the queries?
        benchmark_streaming_random_query(&sbwt, &lcs);
    }

    benchmark_access(&sbwt);
}