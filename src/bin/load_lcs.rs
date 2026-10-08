//! Stand-alone binary that loads an LCS array (as written by `sbwt build --build-lcs`
//! or `sbwt build-lcs`) into memory and counts the number of distinct k'-mers
//! for each k' from 1 to k. Dummy ($-padded) nodes are not treated specially.

use std::io::BufReader;
use sbwt::LcsArray;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <file.lcs> <k>", args[0]);
        std::process::exit(1);
    }
    let path = &args[1];
    let k: usize = args[2].parse().unwrap_or_else(|e| {
        eprintln!("Invalid k {}: {}", args[2], e);
        std::process::exit(1);
    });

    let file = std::fs::File::open(path).unwrap_or_else(|e| {
        eprintln!("Could not open {}: {}", path, e);
        std::process::exit(1);
    });
    let mut reader = BufReader::new(file);

    let start = std::time::Instant::now();
    let lcs = LcsArray::load(&mut reader).unwrap_or_else(|e| {
        eprintln!("Could not load LCS array from {}: {}", path, e);
        std::process::exit(1);
    });

    eprintln!("Loaded LCS array from {} in {:.2?}", path, start.elapsed());
    eprintln!("Length: {}", lcs.len());
    eprintln!("Size in memory: {} bytes", lcs.size_in_bytes());

    // Nodes are in colex order, so a new distinct k'-suffix starts at i iff LCS[i] < k'.
    // Histogram the LCS values, then count(k') = number of LCS values less than k'.
    let mut hist = vec![0usize; k];
    for i in 0..lcs.len() {
        let v = lcs.access(i);
        if v >= k {
            eprintln!("LCS[{}] = {} is not less than k = {}. Wrong k?", i, v, k);
            std::process::exit(1);
        }
        hist[v] += 1;
    }

    let mut count = 0usize;
    for (kp_minus_1, h) in hist.iter().enumerate() {
        count += h;
        println!("{}\t{}", kp_minus_1 + 1, count);
    }
}
