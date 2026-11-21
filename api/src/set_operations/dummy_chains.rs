//! The dummy nodes that the result of an intersection, a difference or a filter needs for its
//! sources: the nodes of the result that lost all their in-neighbors. Such a source `s` needs the
//! dummy chain `$^{k-1} s[0]`, `$^{k-2} s[0..2]`, ..., `$ s[0..k-1]` from the root. A source can
//! be a dummy node `$^j β` itself, which needs the shorter dummy nodes of the chain of β.
//!
//! A dummy node `$^i α` (α a prefix of length k-i of a source) comes in colex order right before
//! the first node that ends with α: from the last character on, they are equal on α, and then `$`
//! is smaller than any character. The nodes that end with α are a colex range that `search(α)`
//! finds, and α ends a node of index1, which has the source with its in-neighbors (k-mers or its
//! own dummy nodes). So the position of each dummy node in a merge plan of the inputs comes from
//! searching the prefixes of the sources, and the dummy nodes that already exist (in index1 or
//! index2) are found by following the edges of the source from the root, as far as they exist.

use rayon::iter::{Either, IntoParallelRefIterator, ParallelIterator};
use rayon::slice::ParallelSliceMut;

use crate::sbwt::SbwtIndex;
use crate::subsetseq::SubsetSeq;
use crate::util::ACGT_TO_0123;

type BitVec = bitvec::vec::BitVec<u64, bitvec::order::Lsb0>;

/// A new dummy node that the sources need: inserted before the node at position `pos` of the
/// merge plan. `edges` has bit c set if the node has an out-edge with character c: to the next node
/// of the dummy chain of a source, or to the source itself.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct NewDummy {
    pub(crate) pos: usize,
    pub(crate) n_dollars: u32,
    pub(crate) edges: u8,
}

/// A needed dummy node that exists, at position `pos` of the merge plan. `edges` has bit c set if
/// a source goes through it with character c.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ExistingDummy {
    pub(crate) pos: usize,
    pub(crate) edges: u8,
}

/// The dummy nodes that the dummy chains of the sources consist of (see [needed_dummies]).
pub(crate) struct NeededDummies {
    /// The new nodes, sorted by position, and for equal positions by decreasing number of dollars
    /// (their colex order before the same node), each node once
    pub(crate) new: Vec<NewDummy>,
    /// The needed nodes that exist, the root included, sorted by position, each node once
    pub(crate) existing: Vec<ExistingDummy>,
}

/// The dummy nodes that the dummy chains of `sources` (nodes of index1 in ascii, with `$` for the
/// dummy characters) consist of, with the edges of all sources through them, in the merge plan whose
/// positions of index1 are `s1` and of index2 are `s2`. Without index2, `s1` is all positions. The
/// root is always needed: the difference gets it only as a needed node. Runs in the current rayon
/// thread pool.
///
/// The existing nodes are listed, not marked in bit vectors over the merge plan: a source needs at
/// most k of them, and there are few sources (a few thousand in a filter of 300 million k-mers),
/// while bit vectors would take 5 bits per position of the merge plan (one for the node, and one
/// per character for its edges), whatever the data.
pub(crate) fn needed_dummies<SS: SubsetSeq + Sync>(index1: &SbwtIndex<SS>, index2: Option<&SbwtIndex<SS>>, s1: &BitVec, s2: Option<&BitVec>, sources: &[Vec<u8>]) -> NeededDummies {
    let k = index1.k();
    let select1 = Select::new(s1);
    let select2 = s2.map(Select::new);
    let existing_node = |pos: usize, c: usize| Either::Right(ExistingDummy { pos, edges: 1 << c });

    let (mut new, mut existing): (Vec<NewDummy>, Vec<ExistingDummy>) = sources.par_iter().flat_map_iter(|source| {
        assert_eq!(source.len(), k);
        // The characters after the dollars of a dummy source
        let chars: Vec<usize> = source.iter().skip_while(|&&c| c == b'$').map(|&c| ACGT_TO_0123[c as usize] as usize).collect();
        assert!(!chars.is_empty(), "the root is not a source");
        let mut new = Vec::<Either<NewDummy, ExistingDummy>>::new();
        // The root, with the edge to the first dummy node of the chain
        new.push(existing_node(0, chars[0]));
        let mut range1 = 0..index1.n_sets();
        let mut range2 = index2.map(|index2| 0..index2.n_sets());
        // The existing dummy nodes of the chain in each input, followed from the root
        let mut node1 = Some(0);
        let mut node2 = index2.map(|_| 0);
        for m in 1..chars.len() {
            // The dummy node with α = the first m characters, and k - m dollars
            let c = chars[m - 1];
            range1 = index1.lf_step(range1.start, c)..index1.lf_step(range1.end, c);
            // A source is a node of index1 (the result has only nodes of index1), and every node of
            // an SBWT but the root has an in-neighbor in it, so index1 has the whole chain of
            // in-neighbors of the source: nodes that end with each of its prefixes
            assert!(!range1.is_empty(), "a prefix of a source does not end a node of index1");
            if let (Some(index2), Some(range)) = (index2, range2.as_mut()) {
                *range = index2.lf_step(range.start, c)..index2.lf_step(range.end, c);
            }
            node1 = node1.filter(|&v| index1.sbwt.set_contains(v, c as u8)).map(|v| index1.lf_step(v, c));
            node2 = match (index2, node2) {
                (Some(index2), Some(v)) => index2.sbwt.set_contains(v, c as u8).then(|| index2.lf_step(v, c)),
                _ => None,
            };
            match (node1, node2) {
                (Some(v), _) => new.push(existing_node(select1.select(v), chars[m])),
                (None, Some(v)) => new.push(existing_node(select2.as_ref().unwrap().select(v), chars[m])),
                (None, None) => {
                    // Before the first node that ends with α, in index1 or in index2
                    let pos1 = select1.select(range1.start);
                    let pos2 = range2.as_ref().filter(|r| !r.is_empty()).map(|r| select2.as_ref().unwrap().select(r.start));
                    new.push(Either::Left(NewDummy { pos: pos2.map_or(pos1, |pos2| pos1.min(pos2)), n_dollars: (k - m) as u32, edges: 1 << chars[m] }));
                },
            }
        }
        new
    }).partition_map(|node| node);

    // Each existing node once, with the edges of all sources through it. The root is always needed.
    existing.push(ExistingDummy { pos: 0, edges: 0 });
    existing.par_sort_unstable_by_key(|e| e.pos);
    let mut n_distinct = 0;
    for i in 0..existing.len() {
        let e = existing[i];
        if n_distinct > 0 && existing[n_distinct - 1].pos == e.pos {
            existing[n_distinct - 1].edges |= e.edges;
            continue;
        }
        existing[n_distinct] = e;
        n_distinct += 1;
    }
    existing.truncate(n_distinct);
    existing.shrink_to_fit();

    // Each new node once, with the edges of all sources through it
    new.par_sort_unstable_by_key(|d| (d.pos, std::cmp::Reverse(d.n_dollars)));
    let mut n_distinct = 0;
    for i in 0..new.len() {
        let d = new[i];
        if n_distinct > 0 {
            let last = &mut new[n_distinct - 1];
            if last.pos == d.pos && last.n_dollars == d.n_dollars {
                last.edges |= d.edges;
                continue;
            }
        }
        new[n_distinct] = d;
        n_distinct += 1;
    }
    new.truncate(n_distinct);
    NeededDummies { new, existing }
}

/// Select on a bit vector: the position of the j-th set bit, from the word of every 1024-th set bit.
/// The select of simple-sds would need a copy of the bit vector (its conversion consumes it), and an
/// index; this keeps only the samples.
struct Select<'a> {
    bits: &'a BitVec,
    samples: Vec<(usize, usize)>, // For every 1024-th set bit: its word, and the number of set bits before that word
}

impl<'a> Select<'a> {
    const SAMPLE: usize = 1024;

    fn new(bits: &'a BitVec) -> Self {
        let mut samples = Vec::new();
        let mut ones_before = 0;
        for (w, &word) in bits.as_raw_slice().iter().enumerate() {
            let ones = Self::word(bits, w, word).count_ones() as usize;
            // The sampled set bits in this word
            while samples.len() * Self::SAMPLE < ones_before + ones {
                samples.push((w, ones_before));
            }
            ones_before += ones;
        }
        Self { bits, samples }
    }

    /// The word, without the bits after the end of the bit vector
    fn word(bits: &BitVec, w: usize, word: u64) -> u64 {
        let end = bits.len() - w * 64;
        if end >= 64 { word } else { word & ((1 << end) - 1) }
    }

    fn select(&self, j: usize) -> usize {
        let (mut w, mut ones_before) = self.samples[j / Self::SAMPLE];
        let words = self.bits.as_raw_slice();
        loop {
            let mut word = Self::word(self.bits, w, words[w]);
            let ones = word.count_ones() as usize;
            if ones_before + ones > j {
                for _ in 0..j - ones_before {
                    word &= word - 1; // Clear the lowest set bit
                }
                return w * 64 + word.trailing_zeros() as usize;
            }
            ones_before += ones;
            w += 1;
        }
    }
}
