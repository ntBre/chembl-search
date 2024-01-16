use std::io::Write;

use openff_toolkit::ForceField;
use rsearch::find_matches;
use rsearch::rdkit::SDMolSupplier;
use rsearch::rdkit::{AromaticityModel, SanitizeFlags};

fn main() {
    let path = "/home/brent/omsf/chembl/chembl_33.sdf";
    let m = SDMolSupplier::new(path);

    let forcefield = "openff-2.1.0.offxml";
    let ff = ForceField::load(forcefield).unwrap();
    let h = ff.get_parameter_handler("ProperTorsions").unwrap();
    let mut params = Vec::new();
    for p in h.parameters() {
        params.push((p.id(), p.smirks()));
    }

    // let mut out = std::fs::File::create("out.smiles").unwrap();
    let mut out = std::io::stdout().lock();

    let want = String::from("t18a");

    for (i, mut mol) in m.into_iter().enumerate() {
        mol.sanitize(
            SanitizeFlags::ALL
                ^ SanitizeFlags::ADJUSTHS
                ^ SanitizeFlags::SETAROMATICITY,
        );
        mol.set_aromaticity(AromaticityModel::MDL);
        mol.assign_stereochemistry();
        mol.add_hs();

        let matches = find_matches(&params, &mol);
        if matches.contains(&want) {
            writeln!(out, "{i} {}", mol.to_smiles()).unwrap();
        }
    }
}
