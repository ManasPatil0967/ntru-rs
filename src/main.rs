pub mod ntru_decrypt;
pub mod ntru_encrypt;
pub mod ntru_util;

use ntru_utils::{gen_rand10, poly_inv, check_prime, pad_arr, arr2str, str2bit, bit2str};
use ntru_encrypt::NTRUencrypt;
use ntru_decrypt::NTRUdecrypt;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

fn main() {
    // Example usage
    let n = 107; // Example N value, adjust as needed
    let p = 3; // Example p value, adjust as needed
    let q = 128; // Example q value, adjust as needed
    let df = 16; // Example df value, adjust as needed
    let dg = 16; // Example dg value, adjust as needed
    let dr = 16; // Example dr value, adjust as needed

    // Generate keypair
    let mut decryptor = NTRUdecrypt::new(n, p, q, df, dg, dr);
    decryptor.gen_pub_priv("keys");

    // Encrypt a message
    let message = "Hello, world!";
    let mut encryptor = NTRUencrypt::new(n, p, q, dr);
    encryptor.readPub("keys.pub");
    encryptor.encryptString(message);
    let encrypted_message = encryptor.Me.unwrap();

    // Decrypt the message
    decryptor.read_priv("keys.priv");
    decryptor.decrypt_string(&encrypted_message);
    let decrypted_message = decryptor.m.unwrap();

    println!("Original Message: {}", message);
    println!("Encrypted Message: {}", encrypted_message);
    println!("Decrypted Message: {}", decrypted_message);
}
