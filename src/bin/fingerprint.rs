use std::fs::read_to_string;

use rsearch::rdkit::ROMol;

fn main() {
    for line in read_to_string("t18b.smiles").unwrap().lines() {
        let mol = ROMol::from_smiles(line);
        let fp = mol.morgan_fingerprint(4);
        dbg!(fp);
    }
}
