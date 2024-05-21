use std::{vec, fs, str::FromStr};

use rand::seq::SliceRandom;

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
            s.push_str(", ");
        }
    }
    s
}

pub fn string_to_array(s: &str) -> Vec<i64> {
    s.split(", ").map(|x| i64::from_str(x).unwrap()).collect()
}

pub fn str_to_bits(s: &str) -> Vec<i64> {
    let mut bits = Vec::new();
    for c in s.chars() {
        let mut b = c as u8;
        for _ in 0..8 {
            bits.push((b & 1) as i64);
            b >>= 1;
        }
    }
    bits
}

pub fn bits_to_str(bits: &Vec<u8>) -> String {
    let mut s = String::new();
    for i in 0..bits.len() / 8 {
        let mut b = 0;
        for j in 0..8 {
            b |= bits[i * 8 + j] << j;
        }
        s.push(b as char);
    }
    s
}

pub fn gen_rand(L: i64, P: i64, M: i64) -> Polynomial {
    let mut R: Vec<i64> = vec![0; L as usize];
    let mut rng = rand::thread_rng();
    for i in 0..L {
        if i < P {
            R[i as usize] = 1;
        } else if i < P + M {
            R[i as usize] = -1;
        }
        else {
            break;
        }
    }
    let mut r_slice = R.as_mut_slice();
    (&mut r_slice).shuffle(&mut rng);
    Polynomial::new(R)
}

pub fn pad_arr(arr: &Vec<i64>, N: i64) -> Vec<i64> {
    // Pad the array with 0s to make it of length N
    // This is left pad or right pad? I think it's right pad
    let mut padded = vec![0; N as usize];
    for i in 0..arr.len() {
        padded[i] = arr[i];
    }
    padded
}

fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if a == 0 {
        (b, 0, 1)
    } else {
        let (gcd, x, y) = extended_gcd(b % a, a);
        (gcd, y - (b / a) * x, x)
    }
}

#[derive(PartialEq)]
pub struct Polynomial {
    pub coeffs: Vec<i64>, // Coefficients of the polynomial
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

    pub fn create_I() -> Self {
        let mut I = vec![0; 762 as usize];
        I[0] = -1;
        I[761] = 1;
        Polynomial::new(I)
    }

    pub fn add(&self, other: &Self) -> Self {
        let max_degree = std::cmp::max(self.degree(), other.degree());
        let mut coeffs = vec![0; max_degree + 1];
    
        for (i, coeff) in self.coeffs.iter().enumerate() {
            coeffs[i] += coeff;
        }
    
        for (i, coeff) in other.coeffs.iter().enumerate() {
            coeffs[i] += coeff;
        }
    
        Polynomial { coeffs }
    }

    pub fn subtract(&self, other: &Self) -> Self {
        let max_degree = std::cmp::max(self.degree(), other.degree());
        let mut coeffs = vec![0; max_degree + 1];
    
        for (i, coeff) in self.coeffs.iter().enumerate() {
            coeffs[i] += coeff;
        }
    
        for (i, coeff) in other.coeffs.iter().enumerate() {
            coeffs[i] -= coeff;
        }
    
        Polynomial { coeffs }
    }

    pub fn multiply(&self, other: &Self) -> Self {
        let total_degree = self.degree() + other.degree();
        let mut coeffs = vec![0; total_degree + 1];
    
        for (i, coeff_i) in self.coeffs.iter().enumerate() {
            for (j, coeff_j) in other.coeffs.iter().enumerate() {
                coeffs[i + j] += coeff_i * coeff_j;
            }
        }
    
        Polynomial { coeffs }
    }

    pub fn divide(&self, other: &Self) -> (Self, Self) {
        let mut quo = Polynomial::zero();
        let mut rem = self.clone();
        let mut curr = rem.degree();

        while curr >= other.degree() {
            if rem.coeffs[curr] == 0 {
                curr -= 1;
                continue;
            }
            let lc_quo = rem.coeffs[curr];
            let lc_rem = rem.coeffs[curr - other.degree() + 1];
            let lead_term = lc_quo / lc_rem;

            quo.coeffs.push(lead_term);
            rem.coeffs = rem.coeffs[..curr - other.degree()].iter().map(|&c| c - lead_term * other.coeffs[0]).collect::<Vec<_>>();
            curr -= other.degree();
        }

        (quo, rem)
    }

    pub fn reduce_coeffs(&mut self, n: i64) -> &mut Self{
        for coeff in &mut self.coeffs {
            *coeff %= n;
        }
        self
    }

    pub fn equal(&self, other: &Self) -> bool {
        self.coeffs == other.coeffs
    }

    pub fn modulus(&self, other: &Self) -> Self {
        let (_, rem) = self.divide(other);
        rem
    }

    pub fn multiply_scalar(&self, scalar: i64) -> Self {
        let coeffs = self.coeffs.iter().map(|&c| c * scalar).collect();
        Polynomial { coeffs }
    }

    pub fn pad(&mut self, N: i64) -> &mut Self{
        while self.coeffs.len() < N as usize {
            self.coeffs.push(0);
        }
        self
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
}

impl Initializer {
    pub fn new () -> Self {
        Initializer {
            p: 3,
            q: 256,
            N: 761,
            df: 61,
            dg: 20,
            dr: 18,
            f: Polynomial::new(vec![0; 761 as usize]),
            g: Polynomial::new(vec![0; 761 as usize]),
            h: Polynomial::new(vec![0; 761 as usize]),
            fp: Polynomial::new(vec![0; 761 as usize]),
            fq: Polynomial::new(vec![0; 761 as usize]),
            I: Polynomial::create_I(),
        }
    }

    pub fn gen_fg(&mut self) {
        let mut max_tries = 100;
        self.f = gen_rand(self.N, self.df, self.df - 1);
        self.g = gen_rand(self.N, self.dg, self.dg);
        let zero = Polynomial::new(vec![0; self.N as usize]);
        while max_tries > 0 {
            let fp_temp = poly_euclid_inv(&self.f, &self.I, self.p);
            let fq_temp = poly_euclid_inv(&self.f, &self.I, self.q);
            if fp_temp != zero && fq_temp != zero {
                self.fp = fp_temp.clone();
                self.fq = fq_temp.clone();
            }
            
            if !self.I.is_zero() {
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
        let contents = fs::read_to_string(filename)?;
        let lines: Vec<&str> = contents.lines().collect();
        self.p = i64::from_str(lines[0].split_whitespace().last().unwrap()).unwrap();
        self.q = i64::from_str(lines[1].split_whitespace().last().unwrap()).unwrap();
        self.N = i64::from_str(lines[2].split_whitespace().last().unwrap()).unwrap();
        self.dr = i64::from_str(lines[3].split_whitespace().last().unwrap()).unwrap();
        self.h.coeffs = lines[4].split_whitespace().skip(3).map(|x| i64::from_str(x).unwrap()).collect();
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
        self.gen_fg();
        self.gen_h();
        self.write_pubkey("pubkey.pub").unwrap();
        self.write_privkey("privkey.priv").unwrap();
    }
}