use std::ops::Range;

use bitvec::bitvec;
use bitvec::order::Lsb0;
use bitvec::vec::BitVec;

use crate::pino::Pred8vPino;
use crate::quad_rank::Base4RankVector;
use crate::subsetseq::SubsetSeq;

#[derive(Clone, Eq, PartialEq, Debug)]
pub struct SubsetCorrectionSets {
    concat: Base4RankVector,
    correction_set_a_pred: Option<Pred8vPino>,
    correction_set_c_pred: Option<Pred8vPino>,
    correction_set_g_pred: Option<Pred8vPino>,
    correction_set_t_pred: Option<Pred8vPino>,
}

impl SubsetCorrectionSets {
    fn rank_with_contains(&self, c: u8, i: usize) -> (usize, bool) {
        let base = self.concat.rank_with_contains(i, c);

        let corr_rank_contains =
            |corr: &Option<Pred8vPino>| corr.as_ref().map_or((0, false), |p| p.rank_contains(i));

        match c {
            0 => {
                let corr = corr_rank_contains(&self.correction_set_a_pred);

                (base.0 - corr.0, base.1 && !corr.1)
            }
            1 => {
                let corr = corr_rank_contains(&self.correction_set_c_pred);

                (base.0 + corr.0, base.1 || corr.1)
            }
            2 => {
                let corr = corr_rank_contains(&self.correction_set_g_pred);

                (base.0 + corr.0, base.1 || corr.1)
            }
            3 => {
                let corr = corr_rank_contains(&self.correction_set_t_pred);

                (base.0 + corr.0, base.1 || corr.1)
            }
            _ => unreachable!(),
        }
    }

    pub fn size_in_bytes(&self) -> usize {
        self.concat.size_in_bytes()
            + self
                .correction_set_a_pred
                .as_ref()
                .map(|p| p.size_in_bytes())
                .unwrap_or(0)
            + self
                .correction_set_c_pred
                .as_ref()
                .map(|p| p.size_in_bytes())
                .unwrap_or(0)
            + self
                .correction_set_g_pred
                .as_ref()
                .map(|p| p.size_in_bytes())
                .unwrap_or(0)
            + self
                .correction_set_t_pred
                .as_ref()
                .map(|p| p.size_in_bytes())
                .unwrap_or(0)
    }
}

impl SubsetSeq for SubsetCorrectionSets {
    fn new(subset_seq: Vec<Vec<u8>>, sigma: usize) -> Self {
        assert_eq!(sigma, 4);
        let n = subset_seq.len();

        let mut rows = vec![
            bitvec![u64, Lsb0; 0; n],
            bitvec![u64, Lsb0; 0; n],
            bitvec![u64, Lsb0; 0; n],
            bitvec![u64, Lsb0; 0; n],
        ];

        for (i, set) in subset_seq.iter().enumerate() {
            for &c in set {
                rows[c as usize].set(i, true);
            }
        }

        Self::new_from_bit_vectors(rows)
    }

    fn new_from_bit_vectors(vecs: Vec<BitVec<u64, Lsb0>>) -> Self {
        /* assert_eq!(vecs.len(), 4);
               log::info!("Inside new from bit vectors");
               let a = &vecs[0];
               let c = &vecs[1];
               let g = &vecs[2];
               let t = &vecs[3];

               assert_eq!(a.len(), c.len());
               assert_eq!(c.len(), g.len());
               assert_eq!(g.len(), t.len());

               let n = a.len();

               let mut concat = Vec::<u8>::with_capacity(n);

               let mut correction_a = Vec::<u64>::with_capacity(n/100);
               let mut correction_c = Vec::<u64>::with_capacity(n/50);
               let mut correction_g = Vec::<u64>::with_capacity(n/20);
               let mut correction_t = Vec::<u64>::with_capacity(n/20);

               for i in 0..n {
                   if i % (1 << 26) == 0 {
                       log::info!("at {i} of {n}. Gone through {:.3}% of the bitvectors", ((i * 100) as f64) / (n as f64));
                   }
                   let a_bit = a[i];
                   let c_bit = c[i];
                   let g_bit = g[i];
                   let t_bit = t[i];

                   let count = (a_bit as u8) + (c_bit as u8) + (g_bit as u8) + (t_bit as u8);

                   if count == 1 {
                       if a_bit {
                           concat.push(0);
                       } else if c_bit {
                           concat.push(1);
                       } else if g_bit {
                           concat.push(2);
                       } else {
                           concat.push(3);
                       }
                   } else {
                       if count == 0 {
                           concat.push(0);
                           correction_a.push(i as u64);
                       } else if a_bit {
                           concat.push(0);

                           if c_bit {
                               correction_c.push(i as u64);
                           }
                           if g_bit {
                               correction_g.push(i as u64);
                           }
                           if t_bit {
                               correction_t.push(i as u64);
                           }
                       } else if c_bit {
                           concat.push(1);

                           if g_bit {
                               correction_g.push(i as u64);
                           }
                           if t_bit {
                               correction_t.push(i as u64);
                           }
                       } else if g_bit {
                           concat.push(2);

                           if t_bit {
                               correction_t.push(i as u64);
                           }
                       }
                   }
               }
        */

        assert_eq!(vecs.len(), 4);

        let a = &vecs[0];
        let c = &vecs[1];
        let g = &vecs[2];
        let t = &vecs[3];

        assert_eq!(a.len(), c.len());
        assert_eq!(c.len(), g.len());
        assert_eq!(g.len(), t.len());

        let n = a.len();

        let mut concat = Vec::<u8>::with_capacity(n);

        // These are estimates only. Vec will automatically grow if necessary.
        let mut correction_a = Vec::<u64>::with_capacity(n / 100);
        let mut correction_c = Vec::<u64>::with_capacity(n / 300);
        let mut correction_g = Vec::<u64>::with_capacity(n / 200);
        let mut correction_t = Vec::<u64>::with_capacity(n / 150);

        // Access the underlying u64 storage. Since these are Lsb0 bitvecs,
        // bit 0 of each u64 corresponds to the first position in that word.
        let a_words = a.as_raw_slice();
        let c_words = c.as_raw_slice();
        let g_words = g.as_raw_slice();
        let t_words = t.as_raw_slice();

        let num_words = (n + 63) / 64;

        for word_idx in 0..num_words {
            let base = word_idx * 64;

            if base % (1 << 30) == 0 {
                log::info!(
                    "at {base} of {n}. Gone through {:.3}% of the bitvectors",
                    if n == 0 {
                        100.0
                    } else {
                        (base as f64 * 100.0) / n as f64
                    }
                );
            }

            // The last u64 may contain unused bits beyond n.
            let bits_in_word = (n - base).min(64);
            let valid_mask = if bits_in_word == 64 {
                u64::MAX
            } else {
                (1u64 << bits_in_word) - 1
            };

            let a_word = a_words[word_idx] & valid_mask;
            let c_word = c_words[word_idx] & valid_mask;
            let g_word = g_words[word_idx] & valid_mask;
            let t_word = t_words[word_idx] & valid_mask;

            // A has priority, then C, then G, then T.
            //
            // correction_a: no base is present
            // correction_c: A and C are present, so A wins
            // correction_g: A/C and G are present, so A/C wins
            // correction_t: A/C/G and T are present, so A/C/G wins
            let correction_a_word = !(a_word | c_word | g_word | t_word) & valid_mask;
            let correction_c_word = a_word & c_word;
            let correction_g_word = (a_word | c_word) & g_word;
            let correction_t_word = (a_word | c_word | g_word) & t_word;

            // Store correction positions.
            //
            // trailing_zeros() finds the next set bit and
            // mask &= mask - 1 removes it.
            let mut mask = correction_a_word;
            while mask != 0 {
                let bit = mask.trailing_zeros() as usize;
                correction_a.push((base + bit) as u64);
                mask &= mask - 1;
            }

            let mut mask = correction_c_word;
            while mask != 0 {
                let bit = mask.trailing_zeros() as usize;
                correction_c.push((base + bit) as u64);
                mask &= mask - 1;
            }

            let mut mask = correction_g_word;
            while mask != 0 {
                let bit = mask.trailing_zeros() as usize;
                correction_g.push((base + bit) as u64);
                mask &= mask - 1;
            }

            let mut mask = correction_t_word;
            while mask != 0 {
                let bit = mask.trailing_zeros() as usize;
                correction_t.push((base + bit) as u64);
                mask &= mask - 1;
            }

            // Determine which base wins at each position.
            //
            // A -> 0
            // C -> 1
            // G -> 2
            // T -> 3
            let selected_c = (!a_word) & c_word;
            let selected_g = (!(a_word | c_word)) & g_word;
            let selected_t = (!(a_word | c_word | g_word)) & t_word;

            // Encode each selected base as a 2-bit value:
            //
            // A = 00
            // C = 01
            // G = 10
            // T = 11
            let low_bits = selected_c | selected_t;
            let high_bits = selected_g | selected_t;

            // Unpack the 64 positions into concat.
            for bit in 0..bits_in_word {
                let value = ((low_bits >> bit) & 1) as u8 | ((((high_bits >> bit) & 1) as u8) << 1);

                concat.push(value);
            }
        }

        // Keep your existing Self construction here.
        // ...

        log::info!(
            "Constructed auxilary data structures needed for construction of the needed structures",
        );

        Self {
            concat: Base4RankVector::from_symbols(&concat),
            correction_set_a_pred: if correction_a.is_empty() {
                None
            } else {
                Some(Pred8vPino::from_sorted(&correction_a))
            },
            correction_set_c_pred: if correction_c.is_empty() {
                None
            } else {
                Some(Pred8vPino::from_sorted(&correction_c))
            },
            correction_set_g_pred: if correction_g.is_empty() {
                None
            } else {
                Some(Pred8vPino::from_sorted(&correction_g))
            },
            correction_set_t_pred: if correction_t.is_empty() {
                None
            } else {
                Some(Pred8vPino::from_sorted(&correction_t))
            },
        }
    }

    fn len(&self) -> usize {
        self.concat.len()
    }

    fn rank(&self, c: u8, i: usize) -> usize {
        let base = self.concat.rank_with_contains(i, c).0;

        let corr_rank =
            |corr: &Option<Pred8vPino>| corr.as_ref().map_or(0, |p| p.rank_contains(i).0);

        match c {
            0 => base - corr_rank(&self.correction_set_a_pred),
            1 => base + corr_rank(&self.correction_set_c_pred),
            2 => base + corr_rank(&self.correction_set_g_pred),
            3 => base + corr_rank(&self.correction_set_t_pred),
            _ => unreachable!(),
        }
    }

    // Returns the bit offset of the set bit of specified rank.
    fn select(&self, c: u8, i: usize) -> Option<usize> {
        if i >= self.rank(c, self.len()) {
            return None;
        }
        // set with rank i cannot occur before i
        let mut lb = i;
        let rank_at_lb = self.rank_with_contains(c, lb);
        if rank_at_lb.0 == i && rank_at_lb.1 {
            return Some(i);
        }
        let mut ub = self.concat.len();
        let mut mid = lb + (ub - lb) / 2;
        while lb < ub {
            let rank_at_mid = self.rank_with_contains(c, mid);
            if rank_at_mid.0 == i {
                // Then the next occurrence of char will be the first set with rank i
                let occurrence = self.next_set_with_char(mid, c);
                return occurrence;
            }
            if rank_at_mid.0 < i {
                lb = mid + 1;
            } else {
                ub = mid;
            }
            mid = lb + (ub - lb) / 2;
        }
        None
    }

    // Returns the bit offset of the set bit of specified rank.
    /* fn select(&self, c: u8, i: usize) -> Option<usize> {
        if i >= self.rank(c, self.len()) {
            return None;
        }
        // set with rank i cannot occur before i
        let mut lb = i;
        let mut rank_at_lb = self.rank(c, lb);
        dbg!(lb, rank_at_lb);
        while rank_at_lb < i {
            // This is a safe skip, we must skip at least this amount of sets
            lb = self.next_set_with_char(lb + (i.saturating_sub(rank_at_lb)), c)?;
            rank_at_lb = self.rank(c, lb);
            dbg!(lb, rank_at_lb);
        }
        dbg!(lb, rank_at_lb);
        Some(lb)
    } */

    fn append_set_to_buf(&self, i: usize, buf: &mut Vec<u8>) {
        for c in 0..4 {
            if self.set_contains(i, c) {
                buf.push(c);
            }
        }
    }

    fn subset_size(&self, i: usize) -> usize {
        let mut result = 0;

        for c in 0..4 {
            if self.set_contains(i, c) {
                result += 1;
            }
        }

        result
    }

    fn set_contains(&self, set_idx: usize, character: u8) -> bool {
        let contains = self.concat.contains(set_idx, character);
        let corr_contains =
            |corr: &Option<Pred8vPino>| corr.as_ref().map_or(false, |p| p.rank_contains(set_idx).1);

        match character {
            0 => contains && !corr_contains(&self.correction_set_a_pred),
            1 => contains || corr_contains(&self.correction_set_c_pred),
            2 => contains || corr_contains(&self.correction_set_g_pred),
            3 => contains || corr_contains(&self.correction_set_t_pred),
            _ => unreachable!(),
        }
    }

    fn build_rank(&mut self) {}

    fn build_select(&mut self) {}

    fn has_select_support(&self) -> bool {
        true
    }

    fn has_rank_support(&self) -> bool {
        true
    }

    fn next_set_with_char(&self, set_idx: usize, c: u8) -> Option<usize> {
        assert!(c < 4);
        assert!(set_idx <= self.len());

        let correction = match c {
            0 => &self.correction_set_a_pred,
            1 => &self.correction_set_c_pred,
            2 => &self.correction_set_g_pred,
            3 => &self.correction_set_t_pred,
            _ => unreachable!(),
        };

        let mut pos = set_idx;

        while pos < self.len() {
            let window_start = pos & !63;

            let concat_bits = self
                .concat
                .window256_starting_from_pos_with_sym_word_aligned(pos, c);
            let correction_bits = correction
                .as_ref()
                .map(|p| p.window256_starting_from_key_word_aligned(pos))
                .unwrap_or([0; 4]);
            //window256_starting_from_key_word_aligned(pos);

            let mut words = [0u64; 4];

            for i in 0..4 {
                words[i] = match c {
                    0 => concat_bits[i] & !correction_bits[i],
                    1 | 2 | 3 => concat_bits[i] | correction_bits[i],
                    _ => unreachable!(),
                };
            }

            // Remove positions before pos in the first word.
            let shift = pos & 63;
            if shift != 0 {
                words[0] &= !((1u64 << shift) - 1);
            }

            for i in 0..4 {
                if words[i] != 0 {
                    if window_start + i * 64 + words[i].trailing_zeros() as usize
                        >= self.concat.len()
                    {
                        return None;
                    }
                    return Some(window_start + i * 64 + words[i].trailing_zeros() as usize);
                }
            }

            // The window covered [window_start, window_start + 255].
            pos = window_start + 256;
        }

        None
    }

    fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize> {
        let mut written = 0;

        written += self.concat.serialize(&mut *out)?;

        written += serialize_pred(&self.correction_set_a_pred, &mut *out)?;
        written += serialize_pred(&self.correction_set_c_pred, &mut *out)?;
        written += serialize_pred(&self.correction_set_g_pred, &mut *out)?;
        written += serialize_pred(&self.correction_set_t_pred, &mut *out)?;

        Ok(written)
    }

    fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self> {
        Ok(Self {
            concat: Base4RankVector::load(&mut *input)?,

            correction_set_a_pred: load_pred(&mut *input)?,
            correction_set_c_pred: load_pred(&mut *input)?,
            correction_set_g_pred: load_pred(&mut *input)?,
            correction_set_t_pred: load_pred(&mut *input)?,
        })
    }

    fn call_on_char_occurrences<F: FnMut(usize)>(
        &self,
        range: Range<usize>,
        c: u8,
        mut callback: F,
    ) {
        assert!(range.len() > 0);
        let mut concat_range = self
            .concat
            .get_words_in_range_for_sym(range.clone(), c as u8)
            .to_vec();
        let correction_range = match c {
            0 => self
                .correction_set_a_pred
                .as_ref()
                .map(|p| p.get_words_in_range(range.clone())),
            1 => self
                .correction_set_c_pred
                .as_ref()
                .map(|p| p.get_words_in_range(range.clone())),
            2 => self
                .correction_set_g_pred
                .as_ref()
                .map(|p| p.get_words_in_range(range.clone())),
            3 => self
                .correction_set_t_pred
                .as_ref()
                .map(|p| p.get_words_in_range(range.clone())),
            _ => unreachable!(),
        };

        let concat_bv = bitvec::slice::BitSlice::<u64, Lsb0>::from_slice_mut(&mut concat_range);
        if let Some(correction_range) = correction_range {
            let correction_bv = bitvec::slice::BitSlice::<u64, Lsb0>::from_slice(&correction_range);

            for i in correction_bv.iter_ones() {
                concat_bv.set(i, !concat_bv[i]);
            }
        }
        /* let concat_bv = bitvec::slice::BitSlice::<u64, Lsb0>::from_slice(&concat_range); */

        for i in concat_bv.iter_ones() {
            if i < (range.start & 63) {
                continue;
            }
            if i - (range.start & 63) >= range.len() {
                break;
            }
            callback(range.start.saturating_sub(range.start & 63) + i);
        }
    }

    fn push_labels_forward(
        &self,
        labels: &[u8],
        sbwt_input_range: Range<usize>,
        output_ranges: Vec<&mut [u8]>,
    ) {
        assert_eq!(labels.len(), sbwt_input_range.len());
        assert_eq!(output_ranges.len(), 4);
        if sbwt_input_range.is_empty() {
            return;
        }

        for (char_idx, output_slice) in output_ranges.into_iter().enumerate() {
            let mut output_offset = 0_usize;

            let mut concat_range = self
                .concat
                .get_words_in_range_for_sym(sbwt_input_range.clone(), char_idx as u8)
                .to_vec();
            let correction_range = match char_idx {
                0 => self
                    .correction_set_a_pred
                    .as_ref()
                    .map(|p| p.get_words_in_range(sbwt_input_range.clone())),
                1 => self
                    .correction_set_c_pred
                    .as_ref()
                    .map(|p| p.get_words_in_range(sbwt_input_range.clone())),
                2 => self
                    .correction_set_g_pred
                    .as_ref()
                    .map(|p| p.get_words_in_range(sbwt_input_range.clone())),
                3 => self
                    .correction_set_t_pred
                    .as_ref()
                    .map(|p| p.get_words_in_range(sbwt_input_range.clone())),
                _ => unreachable!(),
            };

            let concat_bv = bitvec::slice::BitSlice::<u64, Lsb0>::from_slice_mut(&mut concat_range);

            if let Some(correction_range) = correction_range {
                let correction_bv =
                    bitvec::slice::BitSlice::<u64, Lsb0>::from_slice(&correction_range);

                for i in correction_bv.iter_ones() {
                    concat_bv.set(i, !concat_bv[i]);
                }
            }

            for i in concat_bv.iter_ones() {
                if i < (sbwt_input_range.start & 63) {
                    continue;
                }
                if i - (sbwt_input_range.start & 63) >= sbwt_input_range.len() {
                    break;
                }
                output_slice[output_offset] = labels[i - sbwt_input_range
                    .start
                    .saturating_sub(sbwt_input_range.start & !63)];
                output_offset += 1;
            }
            assert_eq!(output_slice.len(), output_offset);
        }
    }
}

fn serialize_pred<W: std::io::Write>(
    pred: &Option<Pred8vPino>,
    out: &mut W,
) -> std::io::Result<usize> {
    let mut written = 1;

    match pred {
        None => out.write_all(&[0])?,
        Some(p) => {
            out.write_all(&[1])?;
            written += p.serialize(out)?;
        }
    }

    Ok(written)
}
fn load_pred<R: std::io::Read>(r: &mut R) -> std::io::Result<Option<Pred8vPino>> {
    let mut flag = [0u8; 1];
    r.read_exact(&mut flag)?;

    if flag[0] == 0 {
        Ok(None)
    } else {
        Ok(Some(Pred8vPino::load(r)?))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_next_set_with_char() {
        // Sets:
        // 0: {A}
        // 1: {C}
        // 2: {A,C}
        // 3: {}
        // 4: {G,T}
        // 5: {A,C,G,T}
        let sets = vec![
            vec![0],
            vec![1],
            vec![0, 1],
            vec![],
            vec![2, 3],
            vec![0, 1, 2, 3],
        ];

        let ss = SubsetCorrectionSets::new(sets, 4);

        // Direct concat hits
        assert_eq!(ss.next_set_with_char(0, 0), Some(0)); // A at 0
        assert_eq!(ss.next_set_with_char(0, 1), Some(1)); // C at 1
        assert_eq!(ss.next_set_with_char(0, 2), Some(4)); // G at 4

        // Start after the first occurrence
        assert_eq!(ss.next_set_with_char(1, 0), Some(2));
        assert_eq!(ss.next_set_with_char(2, 1), Some(2));

        // Empty set corrected by representation:
        // set 3 is stored as A correction and should NOT appear as A.
        assert_eq!(ss.next_set_with_char(3, 0), Some(5));

        // Multi-character correction additions:
        // set 4 has G/T
        assert_eq!(ss.next_set_with_char(3, 2), Some(4));
        assert_eq!(ss.next_set_with_char(3, 3), Some(4));

        // No later occurrences
        assert_eq!(ss.next_set_with_char(6, 0), None);
        assert_eq!(ss.next_set_with_char(6, 1), None);
        assert_eq!(ss.next_set_with_char(6, 2), None);
        assert_eq!(ss.next_set_with_char(6, 3), None);
    }
    #[test]
    fn test_next_set_with_char_across_word_boundary() {
        let mut sets = vec![vec![]; 130];

        sets[63] = vec![0];
        sets[64] = vec![1];
        sets[127] = vec![2];
        sets[128] = vec![3];

        let ss = SubsetCorrectionSets::new(sets, 4);

        assert_eq!(ss.next_set_with_char(60, 0), Some(63));
        assert_eq!(ss.next_set_with_char(60, 1), Some(64));
        assert_eq!(ss.next_set_with_char(60, 2), Some(127));
        assert_eq!(ss.next_set_with_char(60, 3), Some(128));

        assert_eq!(ss.next_set_with_char(64, 0), None);
        assert_eq!(ss.next_set_with_char(65, 1), None);
    }

    #[test]
    fn test_select() {
        let sets = vec![
            vec![0],
            vec![1],
            vec![0, 1],
            vec![],
            vec![2, 3],
            vec![0, 1, 2, 3],
        ];

        let ss = SubsetCorrectionSets::new(sets, 4);

        // A
        assert_eq!(ss.select(0, 0), Some(0));
        assert_eq!(ss.select(0, 1), Some(2));
        assert_eq!(ss.select(0, 2), Some(5));
        assert_eq!(ss.select(0, 3), None);

        // C
        assert_eq!(ss.select(1, 0), Some(1));
        assert_eq!(ss.select(1, 1), Some(2));
        assert_eq!(ss.select(1, 2), Some(5));
        assert_eq!(ss.select(1, 3), None);

        // G
        assert_eq!(ss.select(2, 0), Some(4));
        assert_eq!(ss.select(2, 1), Some(5));
        assert_eq!(ss.select(2, 2), None);

        // T
        assert_eq!(ss.select(3, 0), Some(4));
        assert_eq!(ss.select(3, 1), Some(5));
        assert_eq!(ss.select(3, 2), None);
    }
    #[test]
    fn test_select_sparse() {
        let mut sets = vec![vec![]; 600];

        sets[0] = vec![0];
        sets[63] = vec![0];
        sets[64] = vec![0];
        sets[255] = vec![0];
        sets[256] = vec![0];
        sets[511] = vec![0];
        sets[512] = vec![0];
        sets[599] = vec![0];

        let ss = SubsetCorrectionSets::new(sets, 4);

        assert_eq!(ss.select(0, 0), Some(0));
        assert_eq!(ss.select(0, 1), Some(63));
        assert_eq!(ss.select(0, 2), Some(64));
        assert_eq!(ss.select(0, 3), Some(255));
        assert_eq!(ss.select(0, 4), Some(256));
        assert_eq!(ss.select(0, 5), Some(511));
        assert_eq!(ss.select(0, 6), Some(512));
        assert_eq!(ss.select(0, 7), Some(599));
        assert_eq!(ss.select(0, 8), None);
    }
}
