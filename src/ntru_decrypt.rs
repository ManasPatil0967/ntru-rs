use crate::ntru_util::{Polynomial, gen_rand, poly_euclid_inv, Initializer};
pub struct NtruDecrypt {
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
    fn new() -> NtruDecrypt {
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

    


}