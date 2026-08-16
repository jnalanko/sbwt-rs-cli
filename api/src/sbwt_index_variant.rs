use std::io::Write;

use crate::{LcsArray, PrefixLookupTable, SbwtIndex, StreamingIndex, SubsetMatrix, SubsetSeq, compact_int_vector::CompactIntVector, dbg::Dbg};


/// An enum listing SbwtIndex types built on different subset rank implementations provided in this crate. 
pub enum SbwtIndexVariant {
    SubsetMatrix(SbwtIndex<SubsetMatrix>),

    // This is work on progress
    //SubsetCorrectionSets(SbwtIndex<SubsetCorrectionSets>),
}

/// Calls `sbwt.$method$args` on the [SbwtIndex] wrapped inside whichever variant `$self` is.
/// Used to forward [SbwtIndex] methods to [SbwtIndexVariant] without repeating the match arms
/// for every variant at every method.
macro_rules! forward {
    ($self:expr, $method:ident $args:tt) => {
        match $self {
            SbwtIndexVariant::SubsetMatrix(sbwt) => sbwt.$method $args,
            //SbwtIndexVariant::SubsetCorrectionSets(sbwt) => sbwt.$method $args,
        }
    };
}

/// Forwarding methods to the [SbwtIndex] wrapped inside the variant, so that callers that don't
/// care which [SubsetSeq] implementation is in use don't have to match on the enum themselves.
impl SbwtIndexVariant {

    /// See [SbwtIndex::n_kmers].
    pub fn n_kmers(&self) -> usize { forward!(self, n_kmers()) }

    /// See [SbwtIndex::n_sets].
    pub fn n_sets(&self) -> usize { forward!(self, n_sets()) }

    /// See [SbwtIndex::k].
    pub fn k(&self) -> usize { forward!(self, k()) }

    /// See [SbwtIndex::char_idx].
    pub fn char_idx(&self, c: u8) -> usize { forward!(self, char_idx(c)) }

    /// See [SbwtIndex::alphabet].
    pub fn alphabet(&self) -> &[u8] { forward!(self, alphabet()) }

    /// See [SbwtIndex::interval_of_empty_string].
    pub fn interval_of_empty_string(&self) -> std::ops::Range<usize> { forward!(self, interval_of_empty_string()) }

    /// See [SbwtIndex::inlabel].
    pub fn inlabel(&self, i: usize) -> Option<u8> { forward!(self, inlabel(i)) }

    /// See [SbwtIndex::build_select].
    pub fn build_select(&mut self) { forward!(self, build_select()) }

    /// See [SbwtIndex::lf_step].
    pub fn lf_step(&self, i: usize, char_idx: usize) -> usize { forward!(self, lf_step(i, char_idx)) }

    /// See [SbwtIndex::inverse_lf_step].
    pub fn inverse_lf_step(&self, i: usize) -> Option<usize> { forward!(self, inverse_lf_step(i)) }

    /// See [SbwtIndex::push_kmer_to_vec].
    pub fn push_kmer_to_vec(&self, colex_rank: usize, buf: &mut Vec<u8>) { forward!(self, push_kmer_to_vec(colex_rank, buf)) }

    /// See [SbwtIndex::access_kmer].
    pub fn access_kmer(&self, colex_rank: usize) -> Vec<u8> { forward!(self, access_kmer(colex_rank)) }

    /// See [SbwtIndex::search].
    pub fn search(&self, pattern: &[u8]) -> Option<std::ops::Range<usize>> { forward!(self, search(pattern)) }

    /// See [SbwtIndex::search_from].
    pub fn search_from(&self, interval: std::ops::Range<usize>, pattern: &[u8]) -> Option<std::ops::Range<usize>> { forward!(self, search_from(interval, pattern)) }

    /// See [SbwtIndex::reconstruct_padded_spectrum].
    pub fn reconstruct_padded_spectrum(&self, n_threads: usize) -> Vec<u8> { forward!(self, reconstruct_padded_spectrum(n_threads)) }

    /// See [SbwtIndex::set_lookup_table].
    pub fn set_lookup_table(&mut self, prefix_lookup_table: PrefixLookupTable) { forward!(self, set_lookup_table(prefix_lookup_table)) }

    /// See [SbwtIndex::get_lookup_table].
    pub fn get_lookup_table(&self) -> &PrefixLookupTable { forward!(self, get_lookup_table()) }

    /// Returns the number of sets in the range `[0, i)` that have character `c`. See [SubsetSeq::rank].
    pub fn rank(&self, c: u8, i: usize) -> usize {
        match self {
            SbwtIndexVariant::SubsetMatrix(sbwt) => sbwt.sbwt().rank(c, i),
//            SbwtIndexVariant::SubsetCorrectionSets(sbwt) => sbwt.sbwt().rank(c, i),
        }
    }

    /// See [SbwtIndex::push_all_labels_forward].
    pub fn push_all_labels_forward(&self, labels_in: &[u8], labels_out: &mut [u8], n_threads: usize) { forward!(self, push_all_labels_forward(labels_in, labels_out, n_threads)) }

    /// See [SbwtIndex::push_all_labels_forward_compact].
    pub fn push_all_labels_forward_compact(&self, labels_in: &CompactIntVector<3>, labels_out: &mut CompactIntVector<3>, n_threads: usize) { forward!(self, push_all_labels_forward_compact(labels_in, labels_out, n_threads)) }

    /// See [SbwtIndex::build_last_column].
    pub fn build_last_column(&self) -> Vec<u8> { forward!(self, build_last_column()) }

    /// See [SbwtIndex::build_last_column_compact].
    pub fn build_last_column_compact(&self) -> CompactIntVector<3> { forward!(self, build_last_column_compact()) }

    /// See [SbwtIndex::compute_dummy_node_marks].
    pub fn compute_dummy_node_marks(&self) -> bitvec::vec::BitVec { forward!(self, compute_dummy_node_marks()) }

    /// Builds the LCS array for this index. See [LcsArray::from_sbwt].
    pub fn build_lcs(&self, n_threads: usize, optimize_peak_ram: bool) -> LcsArray {
        match self {
            SbwtIndexVariant::SubsetMatrix(sbwt) => LcsArray::from_sbwt(sbwt, n_threads, optimize_peak_ram),
//            SbwtIndexVariant::SubsetCorrectionSets(sbwt) => LcsArray::from_sbwt(sbwt, n_threads, optimize_peak_ram),
        }
    }

    pub fn build_lookup_table(&self, prefix_len: usize) -> PrefixLookupTable {
        match self {
            SbwtIndexVariant::SubsetMatrix(sbwt_index) => PrefixLookupTable::new(sbwt_index, prefix_len),
//            SbwtIndexVariant::SubsetCorrectionSets(sbwt_index) => PrefixLookupTable::new(sbwt_index, prefix_len),
        }
    }

    /// Loads an index that is wrapped in an enum describing the used subset rank structure type.
    /// The format includes a type identifier so the correct variant can later be loaded with [load_sbwt_index_variant].
    pub fn serialize(&self, out: &mut impl std::io::Write) -> std::io::Result<usize> {
        match self {
            SbwtIndexVariant::SubsetMatrix(sbwt) => {
                out.write_all(&(b"SubsetMatrix".len() as u64).to_le_bytes())?;
                out.write_all(b"SubsetMatrix")?;
                sbwt.serialize(out)
            },
            /*
            SbwtIndexVariant::SubsetCorrectionSets(sbwt) => {
                out.write_all(&(b"SubsetCorrectionSets".len() as u64).to_le_bytes())?;
                out.write_all(b"SubsetCorrectionSets")?;
                sbwt.serialize(out)
            }
            */
        }
    }

    /// Loads an index that was stored with [write_sbwt_index_variant]. This includes a type identifier
    /// to load the correct subset rank variant.
    pub fn load(input: &mut impl std::io::Read) -> Result<Self, Box<dyn std::error::Error>> {
        let type_id_len = byteorder::ReadBytesExt::read_u64::<byteorder::LittleEndian>(input).unwrap();
        let mut type_id = vec![0_u8; type_id_len as usize];
        input.read_exact(&mut type_id)?;

        if type_id == b"SubsetMatrix" {
            Ok(SbwtIndexVariant::SubsetMatrix(SbwtIndex::<SubsetMatrix>::load(input)?))
        } /*else if type_id == b"SubsetCorrectionSets" {
            Ok(SbwtIndexVariant::SubsetCorrectionSets(SbwtIndex::<SubsetCorrectionSets>::load(input)?))
        } */ else {
            Err("Unknown SBWT index type".into())
        }
    }

    pub fn get_streaming_index<'a>(&'a self, lcs: &'a LcsArray) -> StreamingIndex<'a, Self, LcsArray> {
        StreamingIndex{contract_left: lcs, extend_right: self, n: self.n_sets(), k: self.k()}
    }

    // &mut self because needs to build select support if it does not exist.
    pub fn parallel_export_unitigs<W: Write + Send + Sync>(&mut self, out: &mut W, lcs: Option<&LcsArray>, n_threads: usize) {
        match self {
            SbwtIndexVariant::SubsetMatrix(sbwt) => {
                sbwt.build_select();
                Dbg::new(sbwt, lcs, n_threads).parallel_export_unitigs(out, n_threads);
            }
            /*SbwtIndexVariant::SubsetCorrectionSets(sbwt) => {
                sbwt.build_select();
                Dbg::new(sbwt, lcs, n_threads).parallel_export_unitigs(out, n_threads);
            }*/
        }
    }
}

/// This implementation makes it possible to build a streaming index over the SbwtIndexVariant
/// without having to know the type of the underlying index. The downside is that there will be a
/// branch in every right extension call to resolve the type every time. It is a predictable branch,
/// but for best performance, use the inner index type directly as the right extension implementation.
impl crate::streaming_index::ExtendRight for SbwtIndexVariant {
    fn extend_right(&self, I: std::ops::Range<usize>, c: u8) -> std::ops::Range<usize> {
        forward!(&self, extend_right(I, c)) // The match statement happens here
    }
}

//pub fn set_lookup_table(&mut self, prefix_lookup_table: PrefixLookupTable) { forward!(self, set_lookup_table(prefix_lookup_table)) }