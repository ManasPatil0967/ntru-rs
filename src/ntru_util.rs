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

pub struct Polynomial {
    coeffs: Vec<i64>, // Coefficients of the polynomial
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
        let mut current_degree = rem.degree(); // Store the current degree in a mutable variable
    
        while current_degree >= other.degree() {
            let lead_coeff_quo = rem.coeffs[current_degree];
            let lead_coeff_rem = rem.coeffs[current_degree - other.degree() + 1];
            let lead_term = lead_coeff_quo / lead_coeff_rem;
    
            quo.coeffs.push(lead_term);
            rem.coeffs = rem.coeffs[..current_degree - other.degree()].iter().map(|&c| c - lead_term * other.coeffs[0]).collect::<Vec<_>>();
            current_degree -= other.degree(); // Modify the mutable variable instead
        }
    
        (quo, rem)
    }

    pub fn reduce_coeffs(&mut self, n: i64) {
        for coeff in &mut self.coeffs {
            *coeff %= n;
        }
    }
    
    pub fn poly_euclid_inv(f: &Self, g: &Self, n: i64) -> Option<Polynomial> {
        let mut x0 = Polynomial::one(); // Represents 1
        let mut x1 = Polynomial::zero();
        let mut y0 = Polynomial::zero();
        let mut y1 = Polynomial::one();

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

        if a.degree()!= 0 { // Check if a is a unit polynomial (degree 0)
            return None;
        }

        Some(x0)
    }

}
