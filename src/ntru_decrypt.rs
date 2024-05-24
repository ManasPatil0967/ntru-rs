use crate::ntru_util::{bits_to_str, string_to_array, Initializer, Polynomial, Field};
pub struct NtruDecrypt<T: Field> {
    pub p: i64,
    pub q: i64,
    pub N: i64,
    pub df: i64,
    pub dg: i64,
    pub dr: i64,
    pub f: Polynomial<T>,
    pub g: Polynomial<T>,
    pub h: Polynomial<T>,
    pub fp: Polynomial<T>,
    pub fq: Polynomial<T>,
    pub I: Polynomial<T>,
    pub input: String,
    pub output: String,
}

impl<T: Field> NtruDecrypt<T> {
    // Constructor
    pub fn new() -> NtruDecrypt<T> {
        let init = Initializer::new();
        NtruDecrypt {
            p: init.p,
            q: init.q,
            N: init.N,
            df: init.df,
            dg: init.dg,
            dr: init.dr,
            f: init.f,
            g: init.g,
            h: init.h,
            fp: init.fp,
            fq: init.fq,
            I: init.I,
            input: String::new(),
            output: String::new(),
        }
    }

    pub fn decrypt(&mut self, input: String) {
        self.input = input;
        println!("Ciphertext in Decrypt: {}", self.input.clone());
        let c = Polynomial::new(string_to_array(&self.input));
        println!("Ciphertext in Decrypt: {:?}", c.coeffs.clone());
        let cf = c.modulus(&self.fq);
        println!("Decrypted cf: {:?}", cf.coeffs.clone());
        let mut m = cf.multiply(&self.fp);
        println!("Decrypted m: {:?}", m.coeffs.clone());
        let mut mf = m.reduce_coeffs(self.p.into());
        println!("Decrypted mf: {:?}", mf.coeffs.clone());
        let mut binding = mf.modulus(&self.I);
        mf = binding.reduce_coeffs(self.p.into());
        println!("Decrypted mf: {:?}", mf.coeffs.clone());
        self.output = bits_to_str(&mf.coeffs.clone().into_iter().map(|x| x as u8).collect::<Vec<u8>>());
        println!("Decrypted: {}", self.output.clone());
    }


}
