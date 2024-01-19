use std::{collections::HashSet, path::Path};

use openff_qcsubmit::results::OptimizationResultCollection;
use rsearch::rdkit::ROMol;

fn main() {
    let base = Path::new(
        "/home/brent/omsf/projects/valence-fitting/02_curate-data/datasets",
    );
    let sage_tm_opt = base.join("combined-opt.json");
    assert!(sage_tm_opt.exists(), "{:?}", sage_tm_opt);
    let tm_opt = OptimizationResultCollection::parse_file(sage_tm_opt).unwrap();
    let inchis: HashSet<_> = tm_opt
        .entries()
        .into_values()
        .flatten()
        .map(|e| e.inchi_key())
        .collect();
    let my_inchis: HashSet<_> = tm_opt
        .entries()
        .into_values()
        .flatten()
        .map(|e| ROMol::from_smiles(&e.cmiles()).to_inchi_key())
        .collect();
    // my_inchis only contains 1663 instead of 1675
    assert_eq!(inchis.len(), my_inchis.len());
    assert_eq!(inchis, my_inchis);
    println!("hello, world");
}
