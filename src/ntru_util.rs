use std::{vec, fs, str::FromStr};
use num::ToPrimitive as _;
use rand::seq::SliceRandom;
use rand::thread_rng;

pub fn check_prime(n: i64) -> bool {
    if n <= 1 {
        return false;
    }
    if n <= 3 {
        return true;
    }
    if n % 2 == 0 || n % 3 == 0 {
        return false;
    }
    let mut i = 5;
    while i * i <= n {
        if n % i == 0 || n % (i + 2) == 0 {
            return false;
        }
        i += 6;
    }
    true
}

pub fn check_coprime(a: i64, b: i64) -> bool {
    let mut a = a;
    let mut b = b;
    while b != 0 {
        let temp = b;
        b = a % b;
        a = temp;
    }
    a == 1
}

pub fn array_to_string(arr: &Vec<i64>) -> String {
    let mut s = String::new();
    for i in 0..arr.len() {
        s.push_str(&arr[i].to_string());
        if i != arr.len() - 1 {
            s.push_str(" ");
        }
    }
    s
}

pub fn string_to_array(s: &str) -> Vec<i64> {
    s.split(" ").map(|x| i64::from_str(x).unwrap()).collect()
}

pub fn str_to_bits(st: &str) -> Vec<i64> {
    let bytes = st.as_bytes();
    
    let mut num = 0u128;
    for &byte in bytes {
        num = (num << 8) | byte as u128;
    }
    
    // Convert the integer to a binary string
    let bin_str = format!("{:b}", num);
    
    let res: Vec<u8> = bin_str.chars().map(|c| c.to_digit(2).unwrap() as u8).collect();
    
    res.into_iter().map(|x| x as i64).collect()
}

use std::str;

pub fn bits_to_str(bits: &[u8]) -> String {
    let mut bits_vec = bits.to_vec();
    
    while bits_vec.len() % 8!= 0 {
        bits_vec.insert(0, 0);
    }

    let bytes: Vec<u8> = bits_vec.chunks_exact(8).map(|chunk| {
        chunk.iter().rev().enumerate().fold(0, |acc, (i, &bit)| acc | ((bit & 1) << i))
    }).collect();

    bytes
    .iter().map(|&c| c as char).collect::<String>()
}

pub fn rem(vec: &mut Vec<i64>) {
    let mut index = 0;
    // Find the first non-zero element
    for (i, &value) in vec.iter().enumerate() {
        if value!= 0 {
            index = i;
            break;
        }
    }

    vec.drain(..index);
}

pub fn pad_arr(arr: &Vec<i64>, N: i64) -> Vec<i64> {
    let mut padded = vec![0; N as usize];
    for i in 0..arr.len() {
        padded[i] = arr[i];
    }
    padded
}

#[derive(PartialEq)]
pub struct Polynomial {
    pub coeffs: Vec<i64>,
}

impl Polynomial {
    pub fn new(coeffs: Vec<i64>) -> Self {
        Polynomial { coeffs }
    }

    pub fn zero() -> Self {
        Polynomial { coeffs: vec![0] }
    }

    pub fn one() -> Self {
        Polynomial { coeffs: vec![1] }
    }

    pub fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|&c| c == 0)
    }

    pub fn degree(&self) -> usize {
        self.coeffs.len() - 1
    }

    pub fn clone(&self) -> Self {
        Polynomial { coeffs: self.coeffs.clone() }
    }

    pub fn create_i() -> Self {
        let mut I = vec![0; 108 as usize];
        I[0] = 1;
        I[107] = -1;
        Polynomial::new(I)
    }

    pub fn add(&self, other: &Self) -> Self {
        let max_degree = std::cmp::max(self.degree(), other.degree());
        let mut coeffs = vec![0; max_degree + 1];
        for (i, coeff) in other.coeffs.iter().rev().enumerate() {
            coeffs[i] += coeff;
        }

        for (i, coeff) in self.coeffs.iter().rev().enumerate() {
            coeffs[i] += coeff;
        }
        // Reverse the coeffs
        coeffs.reverse();
    
        let mut z = Polynomial { coeffs };
        rem(&mut z.coeffs);
        z
    }

    pub fn subtract(&self, other: &Self) -> Self {
        let max_degree = std::cmp::max(self.degree(), other.degree());
        let mut coeffs = vec![0; max_degree + 1];
    
        for (i, coeff) in self.coeffs.iter().rev().enumerate() {
            coeffs[i] += coeff;
        }
    
        for (i, coeff) in other.coeffs.iter().rev().enumerate() {
            coeffs[i] -= coeff;
        }
        coeffs.reverse();
        let mut z = Polynomial { coeffs };
        rem(&mut z.coeffs);
        z
    }

    pub fn multiply(&self, other: &Self) -> Self {
        let total_degree = self.degree() + other.degree();
        let mut coeffs = vec![0; total_degree + 1];
    
        for (i, coeff_i) in self.coeffs.iter().enumerate() {
            for (j, coeff_j) in other.coeffs.iter().enumerate() {
                coeffs[i + j] += coeff_i * coeff_j;
            }
        }
    
        let mut z = Polynomial { coeffs };
        rem(&mut z.coeffs);
        z
    }

    pub fn divide<'a>(&'a self, other: &'a Self) -> (Self, Self) {
        let mut q_coeffs = vec![0; self.degree() - other.degree() + 1];
        let mut r_coeffs = self.coeffs.clone();

        for k in (0..=self.degree() - other.degree()).rev() {
            let qk = r_coeffs[k + other.degree()] / other.coeffs[other.degree()];
            q_coeffs[k] = qk;

            for j in 0..=other.degree() {
                r_coeffs[k + j] = r_coeffs[k + j] - qk * other.coeffs[j];
            }
        }

        (
            Polynomial::new(q_coeffs),
            Polynomial::new(r_coeffs[..other.degree()].to_vec()),
        )
    }

    pub fn reduce_coeffs(&mut self, n: i64) -> &mut Self{
        for coeff in &mut self.coeffs {
            *coeff %= n;
            if *coeff > n/2 {
                *coeff -= n;
            }
            if *coeff <= -n/2 {
                *coeff += n;
            }
        }
        let mut z = Polynomial { coeffs: self.coeffs.clone() };
        rem(&mut z.coeffs);
        self.coeffs = z.coeffs.clone();
        self
    }

    pub fn equal(&self, other: &Self) -> bool {
        self.coeffs == other.coeffs
    }

    pub fn modulus(&self, other: &Self) -> Self {
        let (_, mut rem) = self.divide(other);
        rem
    }

    pub fn multiply_scalar(&self, scalar: i64) -> Self {
        let coeffs = self.coeffs.iter().map(|&c| c * scalar).collect();
        let mut z = Polynomial { coeffs };
        rem(&mut z.coeffs);
        z
    }

    pub fn pad(&mut self, N: i64) -> &mut Self{
        while self.coeffs.len() < N as usize {
            self.coeffs.push(0);
        }
        self
    }
}

fn ntruprime_inv_poly(a: &Polynomial, modulus: u16) -> Option<Polynomial> {
    let n = a.degree() + 1;
    let mut b = Polynomial::zero(); // Corresponds to NtruIntPoly *b
    b.coeffs[0] = 1;
    let mut c = Polynomial::zero(); // Corresponds to NtruIntPoly *c

    let mut f = a.clone(); // Corresponds to NtruIntPoly *f
    f.pad(n as i64); // Ensure f has enough coefficients

    let mut g = Polynomial::create_i(); // Corresponds to NtruIntPoly *g
    g.coeffs[0] = modulus as i64 - 1;
    g.coeffs[1] = modulus as i64 - 1;

    let mut k = 0;

    loop {
        while f.is_zero() || f.coeffs[0] == 0 {
            f.coeffs.rotate_left(1); 
            c.coeffs.rotate_right(1); 
            k += 1;
            if f.is_zero() {
                return None;
            }
        }

        if f.degree() == 0 {
            let f0_inv_option = modular_inverse(f.coeffs[0], modulus as i64);
            match f0_inv_option {
                Some(f0_inv) => {
                    let b_times_f0_inv = b.multiply_scalar(f0_inv);
                    let mod_poly = Polynomial::new(vec![modulus as i64, 1]);
                    let mut inv = b_times_f0_inv.modulus(&mod_poly); 
                    for _ in 0..k {
                        inv.coeffs.rotate_right(1); 
                    }
                    return Some(inv);
                },
                None => return None,
            }
        }

        if f.degree() < g.degree() {
            std::mem::swap(&mut f, &mut g);
            std::mem::swap(&mut b, &mut c);
        }

        let g0_inv_option = modular_inverse(g.coeffs[0], modulus as i64);
        match g0_inv_option {
            Some(g0_inv) => {
                let u = ((f.coeffs[0] as u64 * g0_inv as u64) % (modulus as u64)) as i64;
                f.subtract(&g.multiply_scalar(u));
                b.subtract(&c.multiply_scalar(u));
            },
            None => return None,
        }
    }
}

fn modular_inverse(a: i64, m: i64) -> Option<i64> {
    let mut mn = (m, 0);
    let mut a = a;
    let mut m = m;

    let mut xy = (0, 1);

    while m > 0 {
        let q = a / m;
        let r = a % m;
        a = m;
        m = r;
        let t = mn.0 - q * mn.1;
        mn = (mn.1, t);
        let t = xy.0 - q * xy.1;
        xy = (xy.1, t);
    }

    if xy.0 < 0 {
        xy.0 += m;
    }

    if xy.0!= 1 {
        None
    } else {
        Some(xy.0)
    }
}

pub fn poly_euclid_inv(f: &Polynomial, g: &Polynomial, n: i64) -> Polynomial {
    let mut x0 = Polynomial::one(); // Represents 1
    let mut x1 = Polynomial::zero();
    let mut y0 = Polynomial::zero();
    let mut a = f.clone();
    let mut b = g.clone();

    while!b.is_zero() {
        let (q, r) = a.divide(&b); // Polynomial division
        a = b;
        b = r;

        let temp_x0 = x0.clone();
        x0 = y0.clone().subtract(&q.multiply(&x1));
        y0 = x1.clone();
        x1 = temp_x0.subtract(&q.multiply(&x1));

        // Reduce coefficients modulo N
        x0.reduce_coeffs(n);
        y0.reduce_coeffs(n);
    }

    if a.degree()!= 0 { 
        return Polynomial::zero();
    }

    x0
}

pub fn trunc(x: i64, q: i64) -> i64 {
    let mut x = x;
    x %= q;
    if x > q / 2 {
        x = q - x;
    }
    x
}

pub struct Initializer {
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
    pub ciphertext: String,
    pub outputpoly: Polynomial,
    pub output: String,
}


impl Initializer {
    pub fn new () -> Self {
        Initializer {
            p: 3,
            q: 64,
            N: 107,
            df: 15,
            dg: 12,
            dr: 5,
            f: Polynomial::new(vec![0; 107 as usize]),
            g: Polynomial::new(vec![0; 107 as usize]),
            h: Polynomial::new(vec![0; 107 as usize]),
            fp: Polynomial::new(vec![0; 107 as usize]),
            fq: Polynomial::new(vec![0; 107 as usize]),
            I: Polynomial::create_i(),
            r: Polynomial::new(vec![0; 107 as usize]),
            message: String::new(),
            ciphertext: String::new(),
            // cipherpoly: Polynomial::new(vec![0; 107 as usize]),
            outputpoly: Polynomial::new(vec![0; 107 as usize]),
            output: String::new(),
        }
    }

    pub fn gen_rand(&self, L: i64, P: i64, M: i64) -> Polynomial {
        let mut R: Vec<i64> = vec![0; L.to_usize().unwrap()];
        let mut rng = thread_rng();

        for i in 0..L {
            if i < P {
                R[i.to_usize().unwrap()] = 1;
            } else if i < P + M {
                R[i.to_usize().unwrap()] = -1;
            } else {
                break;
            }
        }

        let r_slice = R.as_mut_slice();
        r_slice.shuffle(&mut rng);
        Polynomial::new(R)
    }
    pub fn gen_fg(&mut self) {
        let mut max_tries = 100;
        self.f = self.gen_rand(self.N, self.df, self.df - 1);
        self.g = self.gen_rand(self.N, self.dg, self.dg);
        let zero = Polynomial::new(vec![0; self.N as usize]);
    
        while max_tries > 0 {
            // Assuming self.p and self.q are u16 values compatible with ntruprime_inv_poly
            let fp_temp_opt = ntruprime_inv_poly(&self.f, self.p as u16);
            let fq_temp_opt = ntruprime_inv_poly(&self.f, self.q as u16);
    
            if let (Some(fp_temp), Some(fq_temp)) = (fp_temp_opt, fq_temp_opt) {
                if fp_temp!= zero || fq_temp!= zero {
                    self.fp = fp_temp;
                    self.fq = fq_temp;
                    break;
                }
                else {
                    max_tries -= 1;
                }
            }
    
            if!self.I.is_zero() {
                break;
            }
            max_tries -= 1;
        }
    }

    pub fn gen_h(&mut self) {
        self.h = self.fq.clone();
        self.h = self.h.multiply_scalar(self.p);
        self.h = self.h.multiply(&self.g);
        self.h = self.h.modulus(&self.I);
    }

    pub fn write_pubkey(&self, filename: &str) -> std::io::Result<()> {
        let filename = format!("{}.pub", filename);
        let pub_head = format!(
            "p ::: {}\nq ::: {}\nN ::: {}\nd ::: {}\nh :::",
            self.p, self.q, self.N, self.dr
        );
        let h_val = self.h.coeffs.iter().map(|&x| x.to_string()).collect::<Vec<String>>().join(" ");
        let contents = format!("{}\n{}", pub_head, h_val);
        fs::write(filename, contents)
    }


    pub fn read_pubkey(&mut self, filename: &str) -> std::io::Result<()> {
        let filename = format!("{}.pub", filename);
        let mut contents = String::new();
        fs::read_to_string(filename)?.lines().for_each(|line| {
            contents.push_str(line.trim());
            contents.push('\n');
        });

        let lines: Vec<&str> = contents.lines().collect();

        self.p = i64::from_str(lines[0].split_whitespace().last().unwrap()).unwrap();
        self.q = i64::from_str(lines[1].split_whitespace().last().unwrap()).unwrap();
        self.N = i64::from_str(lines[2].split_whitespace().last().unwrap()).unwrap();
        self.dr = i64::from_str(lines[3].split_whitespace().last().unwrap()).unwrap();

        let h_coeffs: Vec<i64> = lines[4]
            .split_whitespace()
            .skip(2)
            .map(|x| i64::from_str(x).unwrap()).collect();

        self.h.coeffs = h_coeffs;

        Ok(())
    }


    pub fn write_privkey(&self, filename: &str) -> std::io::Result<()> {
        let filename = format!("{}.priv", filename);
        let priv_head = format!(
            "p ::: {}\nq ::: {}\nN ::: {}\ndf ::: {}\ndg ::: {}\nd ::: {}\nf/fp/fq/g :::",
            self.p, self.q, self.N, self.df, self.dg, self.dr
        );
        let f_val = self.f.coeffs.iter().map(|&x| x.to_string()).collect::<Vec<String>>().join(" ");
        let fp_val = self.fp.coeffs.iter().map(|&x| x.to_string()).collect::<Vec<String>>().join(" ");
        let fq_val = self.fq.coeffs.iter().map(|&x| x.to_string()).collect::<Vec<String>>().join(" ");
        let g_val = self.g.coeffs.iter().map(|&x| x.to_string()).collect::<Vec<String>>().join(" ");
        let contents = format!("{}\n{}\n{}\n{}\n{}", priv_head, f_val, fp_val, fq_val, g_val);
        fs::write(filename, contents)
    }

    pub fn read_privkey(&mut self, filename: &str) -> std::io::Result<()> {
        let filename = format!("{}.priv", filename);
        let contents = fs::read_to_string(filename)?;
        let lines: Vec<&str> = contents.lines().collect();
        self.p = i64::from_str(lines[0].split_whitespace().last().unwrap()).unwrap();
        self.q = i64::from_str(lines[1].split_whitespace().last().unwrap()).unwrap();
        self.N = i64::from_str(lines[2].split_whitespace().last().unwrap()).unwrap();
        self.df = i64::from_str(lines[3].split_whitespace().last().unwrap()).unwrap();
        self.dg = i64::from_str(lines[4].split_whitespace().last().unwrap()).unwrap();
        self.dr = i64::from_str(lines[5].split_whitespace().last().unwrap()).unwrap();
        self.f.coeffs = lines[7].split_whitespace().map(|x| i64::from_str(x).unwrap()).collect();
        self.fp.coeffs = lines[8].split_whitespace().map(|x| i64::from_str(x).unwrap()).collect();
        self.fq.coeffs = lines[9].split_whitespace().map(|x| i64::from_str(x).unwrap()).collect();
        self.g.coeffs = lines[10].split_whitespace().map(|x| i64::from_str(x).unwrap()).collect();
        Ok(())
    }

    pub fn gen_keys(&mut self) {
        if std::path::Path::new("pubkey.pub").exists() && std::path::Path::new("privkey.priv").exists() {
            self.read_pubkey("pubkey").unwrap();
            self.read_privkey("privkey").unwrap();
        }
        else {
            self.gen_fg();
            self.gen_h();
            println!("Writing keys to files");
            self.write_pubkey("pubkey").unwrap();
            self.write_privkey("privkey").unwrap();
        }
    }

    pub fn gen_r(&mut self) {
        self.r = self.gen_rand(107, self.dr, self.dr);
        let mut z = Polynomial { coeffs:self.r.coeffs.clone() };
        rem(&mut z.coeffs);
        self.r = z;
    }

    pub fn encrypt(&mut self, mut bM: Vec<i64>) {
        let mut ciphertexts: Vec<String> = Vec::new();

        for chunk in bM.chunks(107) {
            self.gen_r();
            let mut m = Polynomial::new(chunk.to_vec());
            while m.coeffs.len() < 107 {
                m.coeffs.insert(0, 0);
            }
            
            let mut binding = self.r.multiply(&self.h);
            let rh = binding.reduce_coeffs(self.q);
            
            let mut binding = rh.add(&m);
            
            
            let x = binding.degree() + 1;
            
            binding = binding.modulus(&self.I);
            
            let y = trunc(x.try_into().unwrap(), 107);
            
            binding.coeffs.rotate_right(y as usize);
            
            let e = binding.reduce_coeffs(self.q);

            let ciphertext = e.coeffs.clone().into_iter().map(|x| x.to_string()).collect::<Vec<String>>().join(" ");
            ciphertexts.push(ciphertext);
        }

        self.ciphertext = ciphertexts.join(" ");       
    }

    pub fn decrypt(&mut self, input: String) {
        self.ciphertext = input;
        let ciphervec: Vec<i64> = self.ciphertext.split_whitespace().map(|x| i64::from_str(x).unwrap()).collect();
        let c = Polynomial::new(ciphervec);
        let mut output_vecs: Vec<Vec<i64>> = Vec::new();

        for chunk in c.coeffs.chunks(107) {
            let fc = &self.f.multiply(&Polynomial::new(chunk.to_vec()));
            let mut fc_modI = fc.modulus(&self.I);
            fc_modI.coeffs.rotate_right(2);
            let a = fc_modI.reduce_coeffs(self.q);
            let b = a.reduce_coeffs(self.p);
            let binding = &self.fp.multiply(&b);
            let mut m = binding.modulus(&self.I);
            m.coeffs.rotate_right(1);
            let ans = m.reduce_coeffs(self.p);
            while ans.coeffs.len() < 107 {
                ans.coeffs.insert(0, 0);
            }
            output_vecs.push(ans.coeffs.clone());
        }

        let output_vec = output_vecs.into_iter().flatten().collect::<Vec<i64>>();
        self.output = bits_to_str(&output_vec.clone().into_iter().map(|x| x as u8).collect::<Vec<u8>>());
        self.output = self.output.replace("\0", "");
        self.outputpoly = Polynomial { coeffs: output_vec };
    }
}
