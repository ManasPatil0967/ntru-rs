use std::collections::HashSet;
use std::convert::TryInto;
use std::ops::{Div, Mul, Rem};
use std::sync::Arc;

use num::{BigInt, Integer, One, Zero};
use num_traits::{Pow, ToPrimitive};
use rand::Rng;
use symengine::{Poly, Symbol, SymEngine};

pub fn check_prime(p: &BigInt) -> bool {
    if p <= &BigInt::one() {
        return false;
    } else if p == &BigInt::from(2) || p == &BigInt::from(3) {
        return true;
    } else {
        for i in BigInt::from(4)..=(p / BigInt::from(2)) {
            if p % i == BigInt::zero() {
                return false;
            }
        }
    }
    true
}

pub fn poly_inv(
    poly_in: &[i64],
    poly_i: &[i64],
    poly_mod: &BigInt,
) -> Option<Vec<i64>> {
    let x = Symbol::new("x");
    let p_poly_i = Poly::from_vec(&x, poly_i);
    let n_poly_i = p_poly_i.degree() + 1;

    if check_prime(poly_mod) {
        match p_poly_i.invert(&Poly::from_vec(&x, poly_in), poly_mod) {
            Ok(inv) => Some(inv.coeffs().to_vec()),
            Err(_) => None,
        }
    } else if poly_mod.is_power_of_two() {
        let mut inv = match p_poly_i.invert(&Poly::from_vec(&x, poly_in), &BigInt::from(2)) {
            Ok(inv) => inv,
            Err(_) => return None,
        };
        let exp = poly_mod.log2().to_usize().unwrap();
        for _ in 1..exp {
            inv = ((BigInt::from(2) * inv - Poly::from_vec(&x, poly_in) * inv.pow(2)) % p_poly_i)
                .trunc(poly_mod);
        }
        inv.set_domain(poly_mod.clone());
        Some(inv.coeffs().to_vec())
    } else {
        None
    }
}

pub fn pad_arr(a_in: &[i64], a_out_size: usize) -> Vec<i64> {
    let mut a = a_in.to_vec();
    a.resize(a_out_size, 0);
    a
}

pub fn gen_rand10(l: usize, p: usize, m: usize) -> Vec<i64> {
    if p + m > l {
        panic!("ERROR: Asking for P+M>L.");
    }
    let mut rng = rand::thread_rng();
    let mut r = vec![0; l];
    for i in 0..l {
        if i < p {
            r[i] = 1;
        } else if i < p + m {
            r[i] = -1;
        } else {
            break;
        }
    }
    r.shuffle(&mut rng);
    r
}

pub fn arr2str(ar: &[i64]) -> String {
    ar.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(" ")
}

fn str2bit(st: &str) -> Vec<i64> {
    let bytes = st.as_bytes();
    let mut bits = Vec::with_capacity(bytes.len() * 8);
    for byte in bytes {
        for i in (0..8).rev() {
            bits.push(((byte >> i) & 1) as i64);
        }
    }
    bits
}

pub fn bit2str(bi: &[i64]) -> String {
    let padded = pad_arr(bi, (bi.len() + 7) / 8 * 8);
    let mut chars = Vec::with_capacity(padded.len() / 8);
    for i in (0..padded.len()).step_by(8) {
        let byte = padded[i..i + 8]
            .iter()
            .rev()
            .fold(0, |acc, &bit| (acc << 1) + bit);
        chars.push(byte as u8);
    }
    String::from_utf8_lossy(&chars).to_string()
}

