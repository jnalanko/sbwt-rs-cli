use read_exact::{self, ReadExactExt};

#[derive(Copy, Clone, PartialEq, Eq, Ord, PartialOrd, Hash, Debug)]
pub struct Kmer{
    data: u64
}

#[derive(Debug)]
pub enum KmerEncodingError{
    InvalidNucleotide(char), // contains the offending char
    TooLong(usize), // Contains the length of the k-mer which was too long
}

// B is the number of u64 in a kmer
// The k-mer will be a (32*B)-mer
// NOTE: k-mer comparison only works correctly for equal-length kmers
// because k-mers are padded with A's at the end to fill B*32 characters
#[derive(Copy, Clone, PartialEq, Eq, Ord, PartialOrd, Hash, Debug)]
pub struct LongKmer<const B: usize>{
    data: [u64; B] // Packed with 2 bits / nucleotide so that bitwise lex comparison is k-mer lex comparison
}

// TODO: always pass these by value since this type is Copy?
impl<const B: usize> LongKmer<B>{

    // If the length of the ASCII string is less than 32*B, the k-mer is padded with A's from the left
    pub fn from_ascii(ascii: &[u8]) -> Result<Self, KmerEncodingError>{
        if ascii.len() > B*32{
            return Err(KmerEncodingError::TooLong(ascii.len()));
        }
        let mut data = [0_u64; B];
        for (i, c) in ascii.iter().enumerate() {
            let bitpair: u64 = match *c{
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => {return Err(KmerEncodingError::InvalidNucleotide(*c as char))}
            };
            let block = i / 32;
            let off = 31 - i % 32;
            //eprintln!("Set {} {} {} {}", c, block, off, bitpair);
            data[block] |= bitpair << (2*off);
        }
        
        Ok(Self{data})
    }

    pub fn set_from_left(&self, i: usize, c: u8) -> Self {
        let pos = i;
        let block = pos / 32;
        let off = 31 - pos % 32;
        let mask = 3_u64 << (2*off);


        let mut data_copy = self.data;
        data_copy[block] = (data_copy[block] & !mask) | ((c as u64) << (2*off));

        Self{data: data_copy}
    }

    pub fn get_from_left(&self, i: usize) -> u8 {
        let pos = i;
        let block = pos / 32;
        let off = 31 - pos % 32;
        //eprintln!("Get {} {} {}", block, off, pos);
        ((self.data[block] >> (2*off)) & 3) as u8
    }

    // Extends with A's at the end
    pub fn right_shift(&self, chars: usize) -> Self{
        // TODO: this could be done without any branching
        let mut new_data = [0_u64; B];
        for block in 0..B{
            let b1 = block + chars / 32; // Which block the first char lands on
            let o1 = (chars % 32) * 2; // Which bit within block the first char lands on
            let b2 = block + (31 + chars) / 32; // Which block the last char lands on
            if b1 < B {
                new_data[b1] |= self.data[block] >> o1;
            }
            if b2 < B {
                let shift = 64 - o1; 

                // shift by 64 is panic
                if shift < 64 {
                    new_data[b2] |= self.data[block] << shift;
                }
            }
        }
        Self{data: new_data}
    }



    pub fn left_shift(&self, chars: usize) -> Self{
        // TODO: this could be done without any branching
        let chars = chars as isize;
        let mut new_data = [0_u64; B];
        for block in 0..(B as isize){
            let b1 = block - (chars + 31) / 32; // Which block the first char lands on
            let o1 = ((32 - (chars % 32)) * 2) % 64; // Which bit within block the first char lands on
            let b2 = block - chars / 32; // Which block the last char lands on
            if b1 >= 0 {
                new_data[b1 as usize] |= self.data[block as usize] >> o1;
            }
            if b2 >= 0 {
                let shift = 64 - o1; 

                // shift by 64 is panic
                if shift < 64 {
                    new_data[b2 as usize] |= self.data[block as usize] << shift;
                }
            }
        }
        Self{data: new_data}
    }


    #[allow(dead_code)]
    pub fn get_u64_data(&self) -> &[u64]{
        &self.data
    }

    pub fn from_u64_data(data: [u64; B]) -> Self{
        Self{data}
    }

    pub fn byte_size() -> usize {
        8*B
    }

    pub fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize>{
        // TODO: maybe use an unsafe single write to avoid the loop
        let mut written = 0;
        for block in self.data.iter(){
            let bytes = block.to_le_bytes();
            out.write_all(&bytes)?;
            written += bytes.len();
        }
        Ok(written)
    }

    // Returns Ok(None) if the stream gives an EOF
    pub fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Option<Self>>{
        // TODO: read with just 1 IO call
        // TODO: this should return an error if could not read 8*B bytes
        // These todos may seem easy but they are not because the const generic support is not good enough yet
        let mut data = [0_u64; B];
        let mut buf = [0_u8; 8];
        for block in data.iter_mut(){
            match input.read_exact_or_eof(&mut buf) {
                Ok(true) => {*block = u64::from_le_bytes(buf);},
                Ok(false) => return Ok(None), // EOF
                Err(e) => return Err(e),
            }
        }
        Ok(Some(Self::from_u64_data(data)))
    }

    pub fn lcp(a: &Self, b: &Self) -> usize{
        for i in 0..B{
            let xor = a.data[i] ^ b.data[i];
            if xor != 0{
                return 32*i + xor.leading_zeros() as usize / 2;
            }
        }
        B*32 // Full k-mer match: 32 nucleotids per block
    }

}

impl<const B: usize> std::fmt::Display for LongKmer<B>{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut s = String::new();
        for i in 0..B*32{
            s.push(match self.get_from_left(i){
                0 => 'A',
                1 => 'C',
                2 => 'G',
                3 => 'T',
                _ => panic!("Invalid character in DNA sequence"),
            });
        }
        write!(f, "{}", s)
    }
}

#[cfg(test)]
mod tests{
    use super::*;

    #[allow(clippy::ptr_arg)]
    fn left_shifted(s: &String) -> String {
        let mut s = s.clone();
        s.remove(0);
        s.push('A');
        s
    }

    #[allow(clippy::ptr_arg)]
    fn right_shifted(s: &String) -> String {
        let mut s = s.clone();
        s.pop();
        s.insert(0, 'A');
        s
    }

    #[test]
    fn test_lcp(){
        let ascii1 = b"ACGTACGTACGTACGTACGTACGTACGTACGTACATGCATTT";
        let ascii2 = b"ACGTACGTACGTACGTACGTACGTACGTACGTACATGCATAT";
        let x = LongKmer::<2>::from_ascii(ascii1).unwrap();        
        let y = LongKmer::<2>::from_ascii(ascii2).unwrap();        
        assert_eq!(LongKmer::<2>::lcp(&x, &y), 40);

        // Equal k-mers
        let ascii1 = b"ACGTACGTACGTACGTACGTACGTACGTACGTACATGCATTTCTAGCTAGCTGATCGATCGA";
        let ascii2 = b"ACGTACGTACGTACGTACGTACGTACGTACGTACATGCATTTCTAGCTAGCTGATCGATCGA";
        let x = LongKmer::<2>::from_ascii(ascii1).unwrap();        
        let y = LongKmer::<2>::from_ascii(ascii2).unwrap();        
        assert_eq!(LongKmer::<2>::lcp(&x, &y), 64);
    }

    #[test]
    #[allow(clippy::nonminimal_bool)]
    fn test_long_kmer(){
        let ascii = b"ACGTACGTACGTACGTACGTACGTACGTACGTACATGCATTT";
        let mut x = LongKmer::<2>::from_ascii(ascii).unwrap();        

        // Setting and getting

        x = x.set_from_left(0, 2);
        x = x.set_from_left(1, 3);
        x = x.set_from_left(62, 1);
        eprintln!("{}", x);

        let mut expected = String::from("GTGTACGTACGTACGTACGTACGTACGTACGTACATGCATTTAAAAAAAAAAAAAAAAAAAACA");
        let actual = format!("{}", x); // This currently uses get_from_left
        assert_eq!(expected, actual);

        // Shifts

        for i in 0..100{
            let our = x.left_shift(i);
            eprintln!("{}", our);
            assert_eq!(expected, format!("{}", our));
            expected = left_shifted(&expected);
        }

        expected = String::from("GTGTACGTACGTACGTACGTACGTACGTACGTACATGCATTTAAAAAAAAAAAAAAAAAAAACA");
        for i in 0..100{
            let our = x.right_shift(i);
            eprintln!("{}", our);
            assert_eq!(expected, format!("{}", our));
            expected = right_shifted(&expected);
        }

        // Comparison
        let x = LongKmer::<2>::from_ascii(b"AATCAGCTAGCTACTATCTACGTACTACGTACGGGCGTACGTAGCA").unwrap();
        let y = LongKmer::<2>::from_ascii(b"AATCAGCTAGCTACTATCTACGTACTACGTACGGGCGTACGTCAGC").unwrap();

        assert!(x < y);

        let x = LongKmer::<2>::from_ascii(b"GGGGAC").unwrap();
        let y = LongKmer::<2>::from_ascii(b"GGGGAC").unwrap();

        assert!(x == y);
        assert!(x <= y);
        assert!(!(x < y));
        assert!(!(x > y));

        // COMPARISON ONLY WORKS FOR EQUAL-LENGTH kmers
        /* 
        let x = LongKmer::<2>::from_ascii(b"GGGGAC").unwrap();
        let y = LongKmer::<2>::from_ascii(b"GGGGACAA").unwrap();

        assert!(x < y);

        let x = LongKmer::<2>::from_ascii(b"AATCAGCTAGCTACTATCTACGTACTACGTACGGGCGTACGTAGCAA").unwrap();
        let y = LongKmer::<2>::from_ascii(b"AATCAGCTAGCTACTATCTACGTACTACGTACGGGCGTACGTAGCA").unwrap();

        assert!(x > y);
        */
    }
}