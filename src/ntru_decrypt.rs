use crate::ntru_util::{bit2str, genRand10, poly_inv, checkPrime};
use std::io::{BufRead, BufReader, Write};
use std::fs::File;
use std::path::Path;
use num::{Integer, traits::Unsigned};
use symengine::Expr;

struct NTRUdecrypt {
    n: usize,
    p: u32,
    q: u32,
    df: usize,
    dg: usize,
    dr: usize,
    f: Vec<i32>,
    fp: Vec<i32>,
    fq: Vec<i32>,
    g: Vec<i32>,
    h: Vec<i32>,
    i: Vec<i32>,
    m: Option<String>,
}

impl NTRUdecrypt {
    fn new(n: usize, p: u32, q: u32, df: usize, dg: usize, d: usize) -> Self {
        let mut i = vec![0; n + 1];
        i[n] = -1;
        i[0] = 1;

        NTRUdecrypt {
            n,
            p,
            q,
            df,
            dg,
            dr: d,
            f: vec![0; n],
            fp: vec![0; n],
            fq: vec![0; n],
            g: vec![0; n],
            h: vec![0; n],
            i,
            m: None,
        }
    }

    fn set_npq(&mut self, n: Option<usize>, p: Option<u32>, q: Option<u32>, df: Option<usize>, dg: Option<usize>, d: Option<usize>) {
        if let Some(n) = n {
            if !checkPrime(n) {
                panic!("ERROR: Input value of N not prime");
            } else {
                if let Some(df) = df {
                    if 2 * df > n {
                        panic!("ERROR: Input N too small compared to default df {}", self.df);
                    }
                }
                if let Some(dg) = dg {
                    if 2 * dg > n {
                        panic!("ERROR: Input N too small compared to default dg {}", self.dg);
                    }
                }
                if let Some(d) = d {
                    if 2 * d > n {
                        panic!("ERROR: Input N too small compared to default dr {}", self.dr);
                    }
                }
                self.n = n;
                self.f = vec![0; n];
                self.fp = vec![0; n];
                self.fq = vec![0; n];
                self.g = vec![0; n];
                self.h = vec![0; n];
                self.i = vec![0; n + 1];
                self.i[n] = -1;
                self.i[0] = 1;
            }
        }

        if let (Some(p), Some(q)) = (p, q) {
            if 8 * p > q {
                panic!("ERROR: We require 8p <= q");
            } else if p.gcd(&q) != 1 {
                panic!("ERROR: Input p and q are not coprime");
            } else {
                self.p = p;
                self.q = q;
            }
        } else if p.is_some() || q.is_some() {
            panic!("Error: Can only set p and q together, not individually");
        }

        if let Some(df) = df {
            if 2 * df > self.n {
                panic!("ERROR: Input df such that 2*df>N");
            } else {
                self.df = df;
            }
        }

        if let Some(dg) = dg {
            if 2 * dg > self.n {
                panic!("ERROR: Input dg such that 2*dg>N");
            } else {
                self.dg = dg;
            }
        }

        if let Some(d) = d {
            if 2 * d > self.n {
                panic!("ERROR: Input dr such that 2*dr>N");
            } else {
                self.dr = d;
            }
        }
    }

    fn inv_f(&mut self) -> bool {
        let fp_tmp = poly_inv(&self.f, &self.i, self.p);
        let fq_tmp = poly_inv(&self.f, &self.i, self.q);
        if !fp_tmp.is_empty() && !fq_tmp.is_empty() {
            self.fp = fp_tmp;
            self.fq = fq_tmp;
            if self.fp.len() < self.n {
                self.fp = vec![0; self.n - self.fp.len()].into_iter().chain(self.fp.into_iter()).cloned().collect();
            }
            if self.fq.len() < self.n {
                self.fq = vec![0; self.n - self.fq.len()].into_iter().chain(self.fq.into_iter()).cloned().collect();
            }
            true
        } else {
            false
        }
    }

    fn gen_fg(&mut self) {
        let max_tries = 100;
        self.g = genRand10(self.n, self.dg, self.dg);
        for _ in 0..max_tries {
            self.f = genRand10(self.n, self.df, self.df - 1);
            if self.inv_f() {
                break;
            }
        }   // If we reach here, we have failed to generate f and g
        panic!("ERROR: Failed to generate f and g after {} tries", max_tries);
    }

    fn get_h(&mut self) {
        let x = Expr::new("x", &[]);
        let p_fq = Expr::from(self.p as i64) * Expr::from_poly(&self.fq);
        let g = Expr::from_poly(&self.g);
        let i = Expr::from_poly(&self.i);
        let h = (p_fq.trunc(self.q as i64) * g).trunc(self.q as i64) % i;
        self.h = h.into_poly().coeffs.into_iter().map(|c| c as i32).collect();
    }

    fn write_pub<P: AsRef<Path>>(&self, filename: P) {
        let mut file = File::create(filename).unwrap();
        let pub_head = format!("p ::: {}\nq ::: {}\nN ::: {}\nd ::: {}\nh :::", self.p, self.q, self.n, self.dr);
        file.write_all(pub_head.as_bytes()).unwrap();
        for h in &self.h {
            file.write_all(format!(" {}", h).as_bytes()).unwrap();
        }
    }

    fn read_pub<P: AsRef<Path>>(&mut self, filename: P) {
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let mut lines = reader.lines().map(|l| l.unwrap());
        self.p = lines.next().unwrap().split(" ").last().unwrap().parse().unwrap();
        self.q = lines.next().unwrap().split(" ").last().unwrap().parse().unwrap();
        self.n = lines.next().unwrap().split(" ").last().unwrap().parse().unwrap();
        self.dr = lines.next().unwrap().split(" ").last().unwrap().parse().unwrap();
        self.h = lines.next().unwrap().split(" ").skip(3).take_while(|s| *s != "").map(|s| s.parse().unwrap()).collect();
        self.i = vec![0; self.n + 1];
        self.i[self.n] = -1;
        self.i[0] = 1;
    }

    fn write_priv<P: AsRef<Path>>(&self, filename: P) {
        let mut file = File::create(filename).unwrap();
        let priv_head = format!("p ::: {}\nq ::: {}\nN ::: {}\ndf ::: {}\ndg ::: {}\nd ::: {}\nf/fp/fq/g :::", self.p, self.q, self.n, self.df, self.dg, self.dr);
        file.write_all(priv_head.as_bytes()).unwrap();
        file.write_all(b"\n").unwrap();
        for f in &self.f {
            file.write_all(format!(" {}", f).as_bytes()).unwrap();
        }
        file.write_all(b"\n").unwrap();
        for fp in &self.fp {
            file.write_all(format!(" {}", fp).as_bytes()).unwrap();
        }
        file.write_all(b"\n").unwrap();
        for fq in &self.fq {
            file.write_all(format!(" {}", fq).as_bytes()).unwrap();
        }
        file.write_all(b"\n").unwrap();
        for g in &self.g {
            file.write_all(format!(" {}", g).as_bytes()).unwrap();
        }
    }

    fn read_priv<P: AsRef<Path>>(&mut self, filename: P) {
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let mut lines = reader.lines().map(|l| l.unwrap());
        self.p = lines.next().unwrap().split(" ").last().unwrap().parse().unwrap();
        self.q = lines.next().unwrap().split(" ").last().unwrap().parse().unwrap();
        self.n = lines.next().unwrap().split(" ").last().unwrap().parse().unwrap();
        self.df = lines.next().unwrap().split(" ").last().unwrap().parse().unwrap();
        self.dg = lines.next().unwrap().split(" ").last().unwrap().parse().unwrap();
        self.dr = lines.next().unwrap().split(" ").last().unwrap().parse().unwrap();
        lines.next(); // Skip header line
        self.f = lines.next().unwrap().split(" ").map(|s| s.parse().unwrap()).collect();
        self.fp = lines.next().unwrap().split(" ").map(|s| s.parse().unwrap()).collect();
        self.fq = lines.next().unwrap().split(" ").map(|s| s.parse().unwrap()).collect();
        self.g = lines.next().unwrap().split(" ").map(|s| s.parse().unwrap()).collect();
        self.i = vec![0; self.n + 1];
        self.i[self.n] = -1;
        self.i[0] = 1;
    }

    fn gen_pub_priv<P: AsRef<Path>>(&mut self, keyfilename: P) {
        self.gen_fg();
        self.get_h();
        self.write_pub(keyfilename.as_ref().with_extension("pub"));
        self.write_priv(keyfilename.as_ref().with_extension("priv"));
    }

    fn decrypt(&self, e: &[i32]) -> Vec<i32> {
        if e.len() > self.n {
            panic!("Encrypted message has degree > N");
        }
        let x = Expr::new("x", &[]);
        let f = Expr::from_poly(&self.f);
        let e = Expr::from_poly(e);
        let i = Expr::from_poly(&self.i);
        let a = (f * e).trunc(self.q as i64) % i;
        let b = a.trunc(self.p as i64);
        let fp = Expr::from_poly(&self.fp);
        let c = (fp * b) % i;
        c.into_poly().coeffs.into_iter().map(|c| c as i32).collect()
    }

    fn decrypt_string(&mut self, e: &str) {
        let me: Vec<i32> = e.split(" ").map(|s| s.parse().unwrap()).collect();
        if me.len() % self.n != 0 {
            panic!("ERROR : Input decrypt string is not integer multiple of N");
        }
        let mut marr = Vec::new();
        for d in (0..me.len() / self.n).map(|d| d * self.n) {
            marr.extend(self.decrypt(&me[d..d + self.n]));
        }
        self.m = Some(bit2str(&marr));
    }
}