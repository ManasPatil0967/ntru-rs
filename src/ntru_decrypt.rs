use crate::ntru_util::{Polynomial, gen_rand};
struct NtruDecrypt {
    // Initialize the NtruDecrypt struct
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
    input: String,
    output: String,
}

impl NtruDecrypt {
    // Constructor
    fn new(p: i64, q: i64, N: i64, df: i64, dg: i64, dr: i64) -> NtruDecrypt {
        NtruDecrypt {
            p,
            q,
            N,
            df,
            dg,
            dr,
            f: Polynomial::new(vec![0; N as usize]),
            g: Polynomial::new(vec![0; N as usize]),
            h: Polynomial::new(vec![0; N as usize]),
            fp: Polynomial::new(vec![0; N as usize]),
            fq: Polynomial::new(vec![0; N as usize]),
            I: Polynomial::new(vec![0; N as usize + 1]),
            input: String::new(),
            output: String::new(),
        }
    }

    pub fn gen_fg(&mut self) {
        let mut max_tries = 100;
        let fp_temp = Polynomial::new(vec![0; self.N as usize]);
        let fq_temp = Polynomial::new(vec![0; self.N as usize]);
        while max_tries > 0 {
            self.f = gen_rand(self.N, self.df, self.df - 1);
            self.g = gen_rand(self.N, self.dg, self.dg);
            self.fp = fp_temp.clone();
            self.fq = fq_temp.clone();
            self.fp = self.f.clone();
            self.fq = self.g.clone();
            if !self.I.is_zero() {
                break;
            }
            max_tries -= 1;
        }
    }


}