use crate::ntru_util::{pad_arr, string_to_array, Initializer, Polynomial, Field};
pub struct NtruEncrypt<T: Field> {
    pub p: T,
    pub q: T,
    pub N: T,
    pub df: T,
    pub dg: i64,
    pub dr: i64,
    pub f: Polynomial<T>,
    pub g: Polynomial<T>,
    pub h: Polynomial<T>,
    pub fp: Polynomial<T>,
    pub fq: Polynomial<T>,
    pub I: Polynomial<T>,
    pub r: Polynomial<T>,
    pub message: String,
    pub cipherpoly: Polynomial<T>,
    pub ciphertext: String,
}

impl<T: Field> NtruEncrypt<T> {
    // Constructor
    pub fn new(p: i64, q: i64, N: i64, df:i64, dg:i64, dr:i64, f:Polynomial<T>, g:Polynomial<T>, h:Polynomial<T>, fp:Polynomial<T>, 
               fq: Polynomial<T>, I:Polynomial<T>) -> NtruEncrypt<T> {
        // println!("h in encrypt: {:?}", init.h.coeffs.clone());
        NtruEncrypt {
            p,
            q,
            N,
            df,
            dg,
            dr,
            f,
            g,
            h,
            fp,
            fq,
            I,
            r: Polynomial::new(vec![0; N as usize]),
            message: String::new(),
            cipherpoly: Polynomial::new(vec![0; N as usize]),
            ciphertext: String::new(),
        }
    }

    pub fn gen_r(&mut self) {
        let init = Initializer::new();
        self.r = init.gen_rand(self.N, self.dr, self.dr);
        // println!("r: {:?}", self.r.coeffs.clone());
    }

    pub fn encrypt(&mut self, message: String) {
        self.message = message;
        self.gen_r();
        // let m_len = self.message.len();
        // println!("Message: {}", &self.message);
        let bM = pad_arr(&string_to_array(&self.message), self.N);
        // println!("bM: {:?}", bM);
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
        // Add a space between each element in the array
        self.ciphertext = e.coeffs.clone().into_iter().map(|x| x.to_string()).collect::<Vec<String>>().join(" ");
        // println!("Ciphertext: {}", self.ciphertext.clone());
        self.cipherpoly = e.clone();
    }

}
