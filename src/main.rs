pub mod ntru_decrypt;
pub mod ntru_util;
pub mod ntru_encrypt;

fn main() {
    println!("Hello, world!");
    let mut ntru_util = ntru_util::Initializer::new();
    ntru_util.gen_keys();
    println!("f in main: {:?}", ntru_util.f.coeffs.clone());
    println!("g in main: {:?}", ntru_util.g.coeffs.clone());
    println!("h in main: {:?}", ntru_util.h.coeffs.clone());
    println!("h in main: {:?}", ntru_util.h.coeffs.clone());
    let message = "0 1 0 0 0 1 0 0".to_string();
    ntru_util.message = message.clone();
    println!("Message: {}", ntru_util.message.clone());
    ntru_util.encrypt(message);
    // println!("Ciphertext: {}", ntru_encrypt.ciphertext.clone());
    // ntru_decrypt.decrypt(ntru_encrypt.ciphertext.clone());
    // println!("Decrypted: {}", ntru_decrypt.output.clone());
}
