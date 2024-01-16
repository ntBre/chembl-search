use openff_toolkit::ForceField;
use rsearch::find_matches;
use rsearch::rdkit::SDMolSupplier;
use rsearch::rdkit::{AromaticityModel, SanitizeFlags};

fn main() {
    let path = "/home/brent/omsf/chembl/chembl_33.sdf";
    let m = SDMolSupplier::new(path);
    // let mut out = File::create("out.smiles").unwrap();

    let forcefield = "openff-2.1.0.offxml";
    let ff = ForceField::load(forcefield).unwrap();
    let h = ff.get_parameter_handler("ProperTorsions").unwrap();
    let mut params = Vec::new();
    for p in h.parameters() {
        params.push((p.id(), p.smirks()));
    }

    for mut mol in m.into_iter().take(5) {
        mol.sanitize(
            SanitizeFlags::ALL
                ^ SanitizeFlags::ADJUSTHS
                ^ SanitizeFlags::SETAROMATICITY,
        );
        mol.set_aromaticity(AromaticityModel::MDL);
        mol.assign_stereochemistry();
        mol.add_hs();

        let matches = find_matches(&params, &mol);
        dbg!(matches);
    }
}
