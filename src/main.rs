#![feature(specialization)]
pub mod ntru_decrypt;
pub mod ntru_util;
pub mod ntru_encrypt;

fn main() {
    println!("Hello, world!");
    let mut ntru_util = ntru_util::Initializer::new();
    ntru_util.gen_keys();
    let message = "0 1 0 0 0 1 0 0".to_string();
    ntru_util.message = message.clone();
    println!("Message: {}", ntru_util.message.clone());
    ntru_util.encrypt(message);
    println!("Ciphertext: {}", ntru_util.ciphertext.clone());
    ntru_util.decrypt(ntru_util.ciphertext.clone());
    println!("Decrypted: {}", ntru_util.output.clone());
}
