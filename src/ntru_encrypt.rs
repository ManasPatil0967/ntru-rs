use crate::ntru_util::{gen_rand, pad_arr, string_to_array, Initializer, Polynomial};
pub struct NtruEncrypt {
    p: i64,
    q: i64,
    N: i64,
    df: i64,
    dg: i64,
    dr: i64,
    f: Polynomial,
    g: Polynomial,
    h: Polynomial,
    fp: Polynomial,
    fq: Polynomial,
    I: Polynomial,
    r: Polynomial,
    message: String,
    cipherpoly: Polynomial,
    ciphertext: String,
}

impl NtruEncrypt {
    // Constructor
    fn new() -> NtruEncrypt {
        let init = Initializer::new();
        NtruEncrypt {
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
            r: Polynomial::new(vec![0; init.N as usize]),
            message: String::new(),
            cipherpoly: Polynomial::new(vec![0; init.N as usize]),
            ciphertext: String::new(),
        }
    }

    pub fn gen_r(&mut self) {
        self.r = gen_rand(self.N, self.dr, self.dr);
    }

    pub fn encrypt(&mut self, message: String) {
        self.message = message;
        self.gen_r();
        let m_len = self.message.len();
        let mut bM = pad_arr(&string_to_array(&self.message), self.N);
        let mut m = Polynomial::new(bM);
        let mut binding = self.r.multiply(&self.h);
        let mut rh = binding.reduce_coeffs(self.q);
        let mut binding = rh.add(&m);
        let mut e = binding.reduce_coeffs(self.q);
        let mut binding = e.modulus(&self.I);
        e = binding.reduce_coeffs(self.q);
        e.pad(self.N);
        self.cipherpoly = e.clone();
    }

}