use std::collections::{HashMap, HashSet};

use openff_toolkit::ForceField;
use rayon::iter::{ParallelBridge, ParallelIterator};
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

    let results: Vec<_> = m
        .into_iter()
        .par_bridge()
        .map(|mut mol| {
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

            let mut res: HashMap<String, Vec<String>> = HashMap::new();
            for pid in matches.intersection(&want) {
                res.entry(pid.to_string())
                    .or_default()
                    .push(mol.to_smiles());
            }
            res
        })
        .collect();

    let mut res: HashMap<String, Vec<String>> = HashMap::new();
    for r in results {
        for (pid, mols) in r {
            for mol in mols {
                res.entry(pid.to_string()).or_default().push(mol);
            }
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
