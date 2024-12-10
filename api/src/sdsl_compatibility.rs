use std::io::{BufReader, Cursor};

use byteorder::{LittleEndian, ReadBytesExt};
use rand::AsByteSliceMut;
use simple_sds_sbwt::serialize::Serialize;
use std::io::Read;

pub fn load_sdsl_bit_vector(input: &mut impl std::io::Read) -> std::io::Result<simple_sds_sbwt::raw_vector::RawVector> {
    // The simple_sds format is: [number of bits][number of words][data]
    // sdsl format is: [number of bits][data]

    // So we need to create a new Read implementation that injects the number of words into the byte stream.

    let n_bits = input.read_u64::<LittleEndian>()?;

    // The length of the serialized data is padded to a multiple of 64 bits.
    let n_bits_plus_pad = ((n_bits + 63) / 64) * 64;

    let n_bytes = n_bits_plus_pad / 8;
    let n_words = n_bytes / 8;

    let mut new_header = [0u64, 0u64]; // bits, words
    new_header[0] = n_bits; // Assumes little-endian byte order
    new_header[1] = n_words;

    let mut modified_input = new_header.as_byte_slice_mut().chain(input);

    simple_sds_sbwt::raw_vector::RawVector::load(&mut modified_input)

}

#[cfg(test)]
mod tests {

    use super::*;
    use hex_literal::hex;
    use simple_sds_sbwt::raw_vector::AccessRaw;

    #[test]
    fn test_load_sdsl_bit_vector(){
        // This is an sdsl vector of length 129 that is all zeroes except v[8] = 1 and v[9] = 1
        let data = hex!("81 00 00 00 00 00 00 00 00 03 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00");

        let v = load_sdsl_bit_vector(&mut std::io::Cursor::new(&data)).unwrap();
        assert!(v.len() == 129);

        for i in 0..129 {
            print!("{}", v.bit(i) as u8);
            if v.bit(i) {
                assert!(i == 8 || i == 9);
            } else {
                assert!(i != 8 && i != 9);
            }
        }
        println!();
    }
}