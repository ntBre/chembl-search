use std::collections::{HashMap, HashSet};

use openff_toolkit::ForceField;
use rsearch::find_matches;
use rsearch::rdkit::{AromaticityModel, SanitizeFlags};
use rsearch::rdkit::{ROMol, SDMolSupplier};

fn load_want(path: &str) -> HashSet<String> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|s| s.trim().to_owned())
        .collect()
}

fn main() {
    let path = "/home/brent/omsf/chembl/chembl_33.sdf";
    let m = SDMolSupplier::new(path);

    let forcefield = "openff-2.1.0.offxml";
    let ff = ForceField::load(forcefield).unwrap();
    let h = ff.get_parameter_handler("ProperTorsions").unwrap();
    let mut params = Vec::new();
    for p in h.parameters() {
        params.push((p.id(), ROMol::from_smarts(&p.smirks())));
    }

    let want = load_want("../want.params");

    let mut res: HashMap<String, Vec<String>> = HashMap::new();
    for mut mol in m.into_iter().take(10000) {
        mol.sanitize(
            SanitizeFlags::ALL
                ^ SanitizeFlags::ADJUSTHS
                ^ SanitizeFlags::SETAROMATICITY,
        );
        mol.set_aromaticity(AromaticityModel::MDL);
        mol.assign_stereochemistry();
        mol.add_hs();

        let matches: HashSet<_> =
            find_matches(&params, &mol).into_iter().collect();

        for pid in matches.intersection(&want) {
            res.entry(pid.to_string())
                .or_default()
                .push(mol.to_smiles());
        }
    }

    for (pid, moles) in res {
        println!("{pid}");
        for mol in moles {
            println!("\t{mol}");
        }
        println!();
    }
}
