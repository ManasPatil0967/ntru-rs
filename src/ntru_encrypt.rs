use crate::ntru_util::{gen_rand, pad_arr, string_to_array, Initializer, Polynomial};
pub struct NtruEncrypt {
    pub p: i64,
    pub q: i64,
    pub N: i64,
    pub df: i64,
    pub dg: i64,
    pub dr: i64,
    pub f: Polynomial,
    pub g: Polynomial,
    pub h: Polynomial,
    pub fp: Polynomial,
    pub fq: Polynomial,
    pub I: Polynomial,
    pub r: Polynomial,
    pub message: String,
    pub cipherpoly: Polynomial,
    pub ciphertext: String,
}

impl NtruEncrypt {
    // Constructor
    pub fn new() -> NtruEncrypt {
        let mut init = Initializer::new();
        init.gen_keys();
        // println!("h in encrypt: {:?}", init.h.coeffs.clone());
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
        println!("r: {:?}", self.r.coeffs.clone());
    }

    pub fn encrypt(&mut self, message: String) {
        self.message = message;
        self.gen_r();
        // let m_len = self.message.len();
        println!("Message: {}", &self.message);
        let bM = pad_arr(&string_to_array(&self.message), self.N);
        println!("bM: {:?}", bM);
        let m = Polynomial::new(bM);
        let mut binding = self.r.multiply(&self.h);
        let rh = binding.reduce_coeffs(self.q);
        // println!("h: {:?}", &self.h.coeffs.clone());
        // println!("rh: {:?}", rh.coeffs.clone());
        let mut binding = rh.add(&m);
        let mut e = binding.reduce_coeffs(self.q);
        let mut binding = e.modulus(&self.I);
        e = binding.reduce_coeffs(self.q);
        e.pad(self.N);
        self.ciphertext = e.coeffs.iter().map(|x| x.to_string()).collect::<String>();
        // println!("Ciphertext: {}", self.ciphertext.clone());
        self.cipherpoly = e.clone();
    }

}