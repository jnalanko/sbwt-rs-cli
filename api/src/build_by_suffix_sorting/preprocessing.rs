use crate::SeqStream;
use super::input_structures::CHAR_TO_INDEX;
use super::input_structures::Lcp;
use super::input_structures::Bwt;

use std::path::PathBuf;
use simple_sds_sbwt::bit_vector::BitVector;
use simple_sds_sbwt::raw_vector::{RawVector, AccessRaw};
use simple_sds_sbwt::ops::{Rank, Select};
use byteorder::ReadBytesExt;

pub fn concatenate_sequences<SS, W>(reader: &mut SS, output: &mut W) -> std::io::Result<()>
where
    SS: SeqStream,
    W: std::io::Write
{
    write!(output, "#")?;
    while let Some(sequence) = reader.stream_next() {
        write!(output, "$")?;
        output.write_all(sequence)?;
    }
    write!(output, "$")?;
    Ok(())
}

pub fn truncate_lcp<R: std::io::Read, const LITTLE_ENDIAN: bool>(input: &mut R, length: usize, k: usize) -> std::io::Result<Lcp> {
    log::info!("[truncate_lcp] begin");
    let k_bit_width = usize::BITS - k.leading_zeros();
    let width = (k_bit_width.div_ceil(u8::BITS) as usize).next_power_of_two();
    let mut bytes = [0_u8; size_of::<u64>()];
    let mut data = Vec::<u8>::with_capacity(length * width);

    while input.read_exact(&mut bytes).is_ok() {
        let number = if LITTLE_ENDIAN {
            u64::from_le_bytes(bytes)
        } else {
            u64::from_be_bytes(bytes)
        };
        let truncated_number = number.min(k as u64);
        let result_bytes = &truncated_number.to_le_bytes()[..width];
        data.extend_from_slice(result_bytes);
    }

    let result = Lcp::new_with_width(data, width);
    log::info!("[truncate_lcp] done");
    Ok(result)
}

pub fn ascii_to_bwt<R: std::io::Read>(input: &mut R, length: usize) -> std::io::Result<Bwt> {
    log::info!("[ascii_to_bwt] begin");
    let mut raw_vectors = [
        RawVector::with_len(length, false), // $
        RawVector::with_len(length, false), // A
        RawVector::with_len(length, false), // C
        RawVector::with_len(length, false), // G
        RawVector::with_len(length, false), // T
    ];

    let mut index: usize = 0;
    while let Ok(byte) = input.read_u8() {
        let char_index = CHAR_TO_INDEX[byte as usize];
        if char_index < raw_vectors.len() {
            raw_vectors[char_index].set_bit(index, true);
        }
        index += 1;
    }

    let mut bit_vectors = raw_vectors.into_iter().map(BitVector::from).collect::<Vec<_>>();
    bit_vectors.iter_mut().for_each(|vector| {
        vector.enable_rank();
        vector.enable_select();
    });

    let result = Bwt::new(bit_vectors);
    log::info!("[ascii_to_bwt] done");
    Ok(result)
}

pub struct SeqReader<'a> {
    paths: &'a [PathBuf],
    next_idx: usize,
    current: Option<jseqio::reader::DynamicFastXReader>,
    local_buf: Vec<u8>,
}

impl<'a> SeqReader<'a> {
    pub fn new(paths: &'a [PathBuf]) -> Self {
        Self {
            paths,
            next_idx: 0,
            current: None,
            local_buf: vec![],
        }
    }
}

impl SeqStream for SeqReader<'_> {
    fn stream_next(&mut self) -> Option<&[u8]> {
        loop {
            if let Some(current) = &mut self.current {
                if let Some(rec) = current.read_next().unwrap() {
                    self.local_buf.clear();
                    self.local_buf.extend_from_slice(rec.seq);

                    self.local_buf.reverse();
                    sanitise(&mut self.local_buf);

                    return Some(&self.local_buf);
                } else {
                    self.current = None;
                }
            }

            if self.next_idx < self.paths.len() {
                let path = &self.paths[self.next_idx];
                self.next_idx += 1;
                self.current = Some(jseqio::reader::DynamicFastXReader::from_file(path).unwrap());
            } else {
                return None;
            }
        }
    }
}

impl<'a> Clone for SeqReader<'a> {
    fn clone(&self) -> Self {
        Self {
            paths: self.paths,
            next_idx: 0,
            current: None,
            local_buf: vec![],
        }
    }
}

pub(crate) fn sanitise(data: &mut [u8]) {
    for k in data {
        if CHAR_TO_INDEX[(*k) as usize] > 5 {
            *k = b'$';
        }
    }
}

