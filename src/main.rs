use std::collections::{HashMap, HashSet};
use std::sync::atomic::AtomicUsize;

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

    let want = load_want("want.params");

    let progress = AtomicUsize::new(0);

    let map_op = |mut mol: ROMol| -> Vec<(String, String)> {
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

        let mut res: Vec<(String, String)> = Vec::new();
        let mut smiles = None;
        for pid in matches.intersection(&want) {
            if smiles.is_none() {
                smiles = Some(mol.to_smiles());
            }
            res.push((pid.to_string(), smiles.clone().unwrap()));
        }
        let cur = progress.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        if cur % 24_000 == 0 {
            eprintln!("{}% complete", cur / 24_000);
        }
        res
    };
    let results: Vec<_> = m.into_iter().par_bridge().flat_map(map_op).collect();

    let mut res: HashMap<String, Vec<String>> = HashMap::new();
    for (pid, mol) in results {
        res.entry(pid.to_string()).or_default().push(mol);
    }

    for (pid, moles) in res {
        println!("{pid}");
        for mol in moles {
            println!("\t{mol}");
        }
        println!();
    }
}
