use byteorder::{LittleEndian, ReadBytesExt};
use rand::AsByteSliceMut;
use simple_sds_sbwt::serialize::Serialize;
use std::io::Read;

// Loads an sdsl::bit_vector
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

/// Loads an sdsl::int_vector<0> (width is determined at runtime)
pub fn load_runtime_len_sdsl_int_vector(input: &mut impl std::io::Read) -> std::io::Result<simple_sds_sbwt::int_vector::IntVector> {
    // The SDSL vector might or might include the bit width. If it does, then the format is:
    // [number of bits as u64][width as u8][data]
    // We assume that the format is this.

    // The simple_sds format is:
    // [number of elements as u64] [width as u64] [RawVector]
    // where RawVector is a bit vector in the format explained in load_sdsl_bit_vector.

    // So now, we need to read the header in sdsl format (9 bytes total), and replace it with
    // the header for simple_sds, including the header for RawVector.

    let n_bits = input.read_u64::<LittleEndian>()?; // Not including the header
    let width = input.read_u8()? as u64;

    assert!(n_bits % width == 0); 
    let n_elements = n_bits / width;

    // The length of the serialized data is padded to a multiple of 64 bits.
    let n_bits_plus_pad = ((n_bits + 63) / 64) * 64;
    let n_bytes = n_bits_plus_pad / 8;
    let n_words = n_bytes / 8;

    let mut new_header = [n_elements, width, n_bits, n_words]; // Assumes little endian byte order
    let mut modified_input = new_header.as_byte_slice_mut().chain(input);

    simple_sds_sbwt::int_vector::IntVector::load(&mut modified_input)
}

#[cfg(test)]
mod tests {

    use super::*;
    use hex_literal::hex;
    use simple_sds_sbwt::{ops::{Access, Vector}, raw_vector::AccessRaw};

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

    #[test]
    fn test_load_empty_sdsl_bit_vector(){

        let data = hex!("00 00 00 00 00 00 00 00");

        let v = load_sdsl_bit_vector(&mut std::io::Cursor::new(&data)).unwrap();
        assert!(v.is_empty());
    }

    #[test]
    fn test_load_sdsl_int_vector(){
        // This is an sdsl bit vector width 20 elements of width 5, such that v[3] = 7 and all other elements are zero
        let data = hex!("64 00 00 00 00 00 00 00 05 00 80 03 00 00 00 00 00 00 00 00 00 00 00 00 00");

        let v = load_runtime_len_sdsl_int_vector(&mut std::io::Cursor::new(&data)).unwrap();
        assert!(v.len() == 20);
        assert!(v.width() == 5);

        for i in 0..v.len() {
            print!("{} ", v.get(i));
            if i == 3 {
                assert!(v.get(i) == 7);
            } else {
                assert!(v.get(i) == 0);
            }
        }
        println!();
    }
}
