use std::fs::File;
use std::io::{BufRead, BufReader};
use std::collections::VecDeque;
use num::integer::Integer;
use num_traits::Pow;
use num::Num;
use crate::ntru_util::{genRand10, padArr, str2bit, arr2str};

struct NTRUencrypt {
    N: usize,
    p: i32,
    q: i32,
    dr: usize,
    g: Vec<i32>,
    h: Vec<i32>,
    r: Vec<i32>,
    m: Vec<i32>,
    e: Vec<i32>,
    I: Vec<i32>,
    readKey: bool,
    Me: Option<String>,
}

impl NTRUencrypt {
    fn new(N: usize, p: i32, q: i32, d: usize) -> Self {
        let mut I = vec![0; N + 1];
        I[N] = -1;
        I[0] = 1;
        NTRUencrypt {
            N,
            p,
            q,
            dr: d,
            g: vec![0; N],
            h: vec![0; N],
            r: vec![0; N],
            m: vec![0; N],
            e: vec![0; N],
            I,
            readKey: false,
            Me: None,
        }
    }

    fn readPub(&mut self, filename: &str) {
        let file = File::open(filename).unwrap();
        let mut reader = BufReader::new(file);
        self.p = reader.lines().next().unwrap().unwrap().split(" ").last().unwrap().parse().unwrap();
        self.q = reader.lines().next().unwrap().unwrap().split(" ").last().unwrap().parse().unwrap();
        self.N = reader.lines().next().unwrap().unwrap().split(" ").last().unwrap().parse().unwrap();
        self.dr = reader.lines().next().unwrap().unwrap().split(" ").last().unwrap().parse().unwrap();
        self.h = reader.lines().next().unwrap().unwrap().split(" ").skip(3).take(self.N).map(|x| x.parse().unwrap()).collect();
        self.I[self.N] = -1;
        self.I[0] = 1;
        self.genr();
        self.readKey = true;
    }

    fn genr(&mut self) {
        self.r = genRand10(self.N, self.dr, self.dr);
    }

    fn setM(&mut self, M: &[i32]) {
        if !self.readKey {
            panic!("ERROR : Public key not read before setting message");
        }
        if M.len() > self.N {
            panic!("ERROR : Message length longer than degree of polynomial ring ideal");
        }
        for &m in M {
            if m < -self.p / 2 || m > self.p / 2 {
                panic!("ERROR : Elements of message must be in [-p/2,p/2]");
            }
        }
        self.m = padArr(M, self.N);
    }

    fn encrypt(&mut self, m: Option<&[i32]>) {
        if !self.readKey {
            panic!("Error : Not read the public key file, so cannot encrypt");
        }
        if let Some(m) = m {
            if m.len() > self.N {
                panic!("\n\nERROR: Polynomial message of degree >= N");
            }
            self.m = padArr(m, self.N);
        }
        let x = num::bigint::BigInt::from(1);
        let mut e = (((Poly::from(&self.r) * Poly::from(&self.h)).trunc(self.q as i32) + Poly::from(&self.m)) % Poly::from(&self.I)).trunc(self.q as i32).all_coeffs();
        e = padArr(&e, self.N);
        self.e = e;
    }

    fn encryptString(&mut self, M: &str) {
        if !self.readKey {
            panic!("Error : Not read the public key file, so cannot encrypt");
        }
        let bM = str2bit(M);
        let mut bM = padArr(&bM, bM.len() - bM.len() % self.N + self.N);
        self.Me = Some(String::new());

        for E in 0..bM.len() / self.N {
            self.genr();
            self.setM(&bM[E * self.N..(E + 1) * self.N]);
            self.encrypt(None);
            self.Me += &arr2str(&self.e) + " ";
        }
    }
}

struct Poly {
    coeffs: VecDeque<i32>,
}

impl Poly {
    fn from(arr: &[i32]) -> Self {
        Poly { coeffs: arr.to_vec().into() }
    }

    fn trunc(&self, q: i32) -> Self {
        Poly { coeffs: self.coeffs.iter().map(|x| x.rem_euclid(q)).collect() }
    }

    fn all_coeffs(&self) -> Vec<i32> {
        self.coeffs.iter().cloned().collect()
    }

    fn modulo(&self, other: &Poly) -> Self {
        let mut result = self.clone();
        while result.coeffs.len() >= other.coeffs.len() {
            let coeff = result.coeffs.front().unwrap() / other.coeffs.front().unwrap();
            let mut temp = other.clone();
            temp.coeffs.rotate_left(result.coeffs.len() - other.coeffs.len());
            temp.coeffs.iter_mut().for_each(|x| *x *= coeff);
            result.coeffs = (0..temp.coeffs.len()).map(|i| result.coeffs[i] - temp.coeffs[i]).collect();
            result.coeffs.drain(0..1);
        }
        result
    }
}

impl std::ops::Mul<Poly> for Poly {
    type Output = Poly;

    fn mul(self, rhs: Poly) -> Self::Output {
        let mut result = Poly { coeffs: VecDeque::with_capacity(self.coeffs.len() + rhs.coeffs.len() - 1) };
        for i in 0..self.coeffs.len() {
            for j in 0..rhs.coeffs.len() {
                result.coeffs.push_back(self.coeffs[i] * rhs.coeffs[j]);
            }
            result.coeffs.push_front(0);
        }
        result.coeffs.drain(0..1);
        result
    }
}

impl std::ops::Add<Poly> for Poly {
    type Output = Poly;

    fn add(self, rhs: Poly) -> Self::Output {
        let mut result = self;
        for i in 0..rhs.coeffs.len() {
            result.coeffs[i] += rhs.coeffs[i];
        }
        result
    }
}

impl std::ops::Rem<Poly> for Poly {
    type Output = Poly;

    fn rem(self, rhs: Poly) -> Self::Output {
        self.modulo(&rhs)
    }
}