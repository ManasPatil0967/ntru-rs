use std::fs;

pub mod ntru_decrypt;
pub mod ntru_util;
pub mod ntru_encrypt;

fn main() {
    let mut init = ntru_util::Initializer::new();
    init.gen_keys();
    let message = "Manas Patil.1".to_string();
    let msg_array = ntru_util::str_to_bits(&message);
    init.message = message.clone();
    init.encrypt(msg_array.clone());
    println!("Ciphertext: {}", init.ciphertext.clone());
    fs::write("ciphertext.txt", init.ciphertext.clone()).unwrap();
    let cipher = fs::read_to_string("ciphertext.txt").unwrap();
    init.decrypt(cipher.clone());
    println!("Decrypted: {:?}", init.output.clone());
}
