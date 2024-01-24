//! parse the input SDF file and write out the smiles

use std::fs::File;
use std::io::Write;

use rsearch::rdkit::SDMolSupplier;

fn main() {
    let args: Vec<_> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: parse <input.sdf> <output.smiles>");
        std::process::exit(1);
    }
    let mut output = File::create(args[2].as_str()).unwrap();
    for smile in SDMolSupplier::new(args[1].as_str()).map(|m| m.to_smiles()) {
        writeln!(output, "{smile}").unwrap();
    }
}
