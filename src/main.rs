pub mod ntru_decrypt;
pub mod ntru_util;
pub mod ntru_encrypt;

fn main() {
    println!("Hello, world!");
    let mut init = ntru_util::Initializer::new();
    println!("N: {}", init.N);
    init.gen_keys();
    let message = "1 0 0 1 0 0 0 0 1 1 0 0 1 0 1 0 1 1 0 1 1 0 0 0 1 1 0 1 1 0 0 0 1 1 0 1 1 1 1".to_string();
    init.message = message.clone();
    println!("Message: {}", init.message.clone());
    init.encrypt(message.clone());
    println!("Ciphertext: {}", init.ciphertext.clone());
    // let ciphervec: Vec<i64> = vec! [-3, -31, 17, -4, 1, -9, 23, -7, -28, -20, -29, -17, 15, 23, -16, -13, -26, -2, 17, -25, -19, 14, 17, -27, -6, -17, 11, 32, 5, -16, -8, -24, 19, -26, -4, -4, 24, 15, -26, -11, 11, -25, 0, -5, 25, 28, 29, -16, -25, -3, 27, 27, 3, 20, -16, 7, 20, -14, -9, 4, -19, -21, -13, 8, -11, 14, -28, 11, -27, 3, -11, 4, -28, 28, 8, 19, -24, 15, -10, 9, 32, -9, 20, -21, -20, 4, 23, -4, 12, 20, -16, 26, -12, 29, 28, 11, 7, -14, -5, -18, -9, 7, -22, 30, -14, 27, 14];
    // let cipher = ciphervec.into_iter().map(|x| x.to_string()).collect::<Vec<String>>().join(" ");
    init.decrypt(init.ciphertext.clone());
    // println!("Decrypted: {:?}", init.outputpoly.coeffs.clone());
    let bM = ntru_util::string_to_array(&message);
    let m = ntru_util::Polynomial::new(bM);
    if m.coeffs == init.outputpoly.coeffs {
        println!("Decryption successful!");
    } else {
        println!("Decryption failed!");
    }
}
