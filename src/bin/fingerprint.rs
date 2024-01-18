use std::fs::read_to_string;

use rsearch::rdkit::{fingerprint::tanimoto, ROMol};

fn main() {
    // TODO run DBScan on the distance matrix to find clusters
    let mut fps = Vec::new();
    for smiles in read_to_string("t18b.smiles").unwrap().lines() {
        let mol = ROMol::from_smiles(smiles);
        let fp = mol.morgan_fingerprint_bit_vec::<1024>(4);
        fps.push(fp);
    }

    for i in 0..fps.len() {
        for j in 0..fps.len() {
            print!("{:8.4}", tanimoto(&fps[i], &fps[j]));
        }
        println!();
    }
}
