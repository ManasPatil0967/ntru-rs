pub mod ntru_decrypt;
pub mod ntru_util;
pub mod ntru_encrypt;

fn main() {
    let mut ntru_util = ntru_util::Initializer::new();
    ntru_util.gen_keys();
    let mut ntru_encrypt = ntru_encrypt::NtruEncrypt::new();
    println!("h in main: {:?}", ntru_encrypt.h.coeffs.clone());
    let mut ntru_decrypt = ntru_decrypt::NtruDecrypt::new();
    let message = "010001001".to_string();
    ntru_encrypt.encrypt(message);
    println!("Ciphertext: {}", ntru_encrypt.ciphertext.clone());
    ntru_decrypt.decrypt(ntru_encrypt.ciphertext.clone());
    println!("Decrypted: {}", ntru_decrypt.output.clone());
}