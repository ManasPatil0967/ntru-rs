use std::{vec, fs, str::FromStr};
use std::ops::{Add, Mul, Div, Sub};
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
    s.split(" ").map(|x| i64::from_str(x).unwrap()).collect()
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


// Trait for implementing field operations
pub trait Field: Sized + Copy + Add<Output = Self> + Mul<Output = Self> + Div<Output = Self> + Sub<Output = Self> {
    fn zero() -> Self;
    fn one() -> Self;
}

// Implementation for f64 as an example
impl Field for i64 {
    fn zero() -> Self { 0 }
    fn one() -> Self { 1 }
}

#[derive(PartialEq)]
pub struct Polynomial<T: Field> {
    pub coeffs: Vec<T>, // Coefficients of the polynomial
}

impl<T: Field> Polynomial<T> {
    pub fn new(coeffs: Vec<T>) -> Self {
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
        let mut I = vec![0; 108 as usize];
        I[0] = -1;
        I[107] = 1;
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

    pub fn divide<'a>(&'a self, other: &'a Polynomial<T>) -> (Polynomial<T>, Polynomial<T>) {
        let mut q_coeffs = vec![T::zero(); self.degree() - other.degree() + 1];
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
        }
        self
    }

    pub fn equal(&self, other: &Self) -> bool {
        self.coeffs == other.coeffs
    }

    pub fn modulus(&self, other: &Self) -> Self {
        println!("Divide used in modulus");
        // println!("Self: {:?}", self.coeffs.clone());
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

fn ntruprime_inv_poly<T: Field>(a: &Polynomial<T>, modulus: u16) -> Option<Polynomial<T>> {
    let n = a.degree() + 1;
    let mut b = Polynomial::zero(); // Corresponds to NtruIntPoly *b
    b.coeffs[0] = 1;
    let mut c = Polynomial::zero(); // Corresponds to NtruIntPoly *c

    let mut f = a.clone(); // Corresponds to NtruIntPoly *f
    f.pad(n as i64); // Ensure f has enough coefficients

    let mut g = Polynomial::create_I(); // Corresponds to NtruIntPoly *g
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

// Helper function to compute the modular inverse
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

pub fn poly_euclid_inv<T: Field>(f: &Polynomial<T>, g: &Polynomial<T>, n: i64) -> Polynomial<T> {
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

pub struct Initializer<T: Field> {
    pub p: T,
    pub q: T,
    pub N: T,
    pub df: T,
    pub dg: T,
    pub dr: T,
    pub f: Polynomial<T>,
    pub g: Polynomial<T>,
    pub h: Polynomial<T>,
    pub fp: Polynomial<T>,
    pub fq: Polynomial<T>,
    pub I: Polynomial<T>,
    pub r: Polynomial<T>,
    pub message: String,
    pub ciphertext: String,
    pub cipherpoly: Polynomial<T>,
    pub output: String,
}


impl<T: Field> Initializer<T>{
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
            I: Polynomial::create_I(),
            r: Polynomial::new(vec![0; 107 as usize]),
            message: String::new(),
            ciphertext: String::new(),
            cipherpoly: Polynomial::new(vec![0; 107 as usize]),
            output: String::new(),
        }
    }

    pub fn gen_fg(&mut self) {
        let mut max_tries = 100;
        self.f = gen_rand(self.N, self.df, self.df - 1);
        self.g = gen_rand(self.N, self.dg, self.dg);
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

        // Assuming the first four lines contain the values p, q, N, and d
        self.p = i64::from_str(lines[0].split_whitespace().last().unwrap()).unwrap();
        self.q = i64::from_str(lines[1].split_whitespace().last().unwrap()).unwrap();
        self.N = i64::from_str(lines[2].split_whitespace().last().unwrap()).unwrap();
        self.dr = i64::from_str(lines[3].split_whitespace().last().unwrap()).unwrap();

        // Parsing the coefficients of h
        // let h_coeffs_line: &str = lines[4];
        let h_coeffs: Vec<i64> = lines[4]
            .split_whitespace()
            .skip(2)
            .map(|x| i64::from_str(x).unwrap()).collect();

        self.h.coeffs = h_coeffs;

        // println!("h: {:?}", self.h.coeffs);

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
        // println!("f: {:?}", self.f.coeffs);
        // println!("fp: {:?}", self.fp.coeffs);
        // println!("fq: {:?}", self.fq.coeffs);
        // println!("g: {:?}", self.g.coeffs);
        Ok(())
    }

    pub fn gen_keys(&mut self) {
        // If pubkey.pub and privley.priv exist, then use read_pubkey and read_privkey
        // Otherwise, use write_pubkey and write_privkey
        if std::path::Path::new("pubkey.pub").exists() && std::path::Path::new("privkey.priv").exists() {
            // println!("Reading keys from files");
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
        println!("Keys generated");
    }

    pub fn gen_r(&mut self) {
        self.r = gen_rand(self.N, self.dr, self.dr);
        // println!("r: {:?}", self.r.coeffs.clone());
    }

    pub fn encrypt(&mut self, message: String) {
        self.message = message;
        self.gen_r();
        // let m_len = self.message.len();
        // println!("Message: {}", &self.message);
        let bM = pad_arr(&string_to_array(&self.message), self.N);
        println!("bM done");
        // println!("bM: {:?}", bM);
        let m = Polynomial::new(bM);
        let mut binding = self.r.multiply(&self.h);
        println!("r*h done");
        let rh = binding.reduce_coeffs(self.q);
        println!("r*h reduced");
        // println!("h: {:?}", &self.h.coeffs.clone());
        println!("rh: {:?}", rh.coeffs.clone());
        let mut binding = rh.add(&m);
        println!("r*h + m done");
        let mut e = binding.reduce_coeffs(self.q);
        println!("r*h + m reduced");
        let mut binding = e.modulus(&self.I);
        println!("Modulus done");
        e = binding.reduce_coeffs(self.q);
        println!("Modulus reduced");
        e.pad(self.N);
        println!("Ciphertext: {:?}", e.coeffs.clone());
        // Add a space between each element in the array
        self.ciphertext = e.coeffs.clone().into_iter().map(|x| x.to_string()).collect::<Vec<String>>().join(" ");
        // println!("Ciphertext: {}", self.ciphertext.clone());
        self.cipherpoly = e.clone();
    }

    pub fn decrypt(&mut self, input: String) {
        // self.input = input;
        println!("Ciphertext in Decrypt: {}", input.clone());
        let c = Polynomial::new(string_to_array(&input.clone()));
        println!("Cipherpoly in Decrypt: {:?}", c.coeffs.clone());
        let cf = c.modulus(&self.fq);
        println!("Decrypt cf: {:?}", cf.coeffs.clone());
        let mut m = cf.multiply(&self.fp);
        println!("Decrypt m: {:?}", m.coeffs.clone());
        let mut mf = m.reduce_coeffs(self.p.into());
        println!("Decrypt mf: {:?}", mf.coeffs.clone());
        let mut binding = mf.modulus(&self.I);
        mf = binding.reduce_coeffs(self.p.into());
        println!("Decrypt mf: {:?}", mf.coeffs.clone());
        self.output = bits_to_str(&mf.coeffs.clone().into_iter().map(|x| x as u8).collect::<Vec<u8>>());
        println!("Decrypt: {}", self.output.clone());
    }
}
