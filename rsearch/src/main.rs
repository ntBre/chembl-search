use std::collections::HashMap;

use openff_toolkit::ForceField;
use rdkit::{find_smarts_matches, ROMol, SDMolSupplier};

/// TODO move this to its own crate, possibly in a workspace with rdkit-sys
pub mod rdkit {
    use std::ffi::{CStr, CString};

    use rdkit_sys::{
        RDKit_MolToSmiles, RDKit_ROMol, RDKit_ROMol_delete,
        RDKit_SDMolSupplier, RDKit_create_mol_supplier,
        RDKit_delete_mol_supplier, RDKit_mol_supplier_at_end,
        RDKit_mol_supplier_next,
    };

    pub struct SDMolSupplier(*mut RDKit_SDMolSupplier);

    impl SDMolSupplier {
        /// construct an [SDMolSupplier] from a filepath that can be converted
        /// to a CString. panics if this conversion fails
        pub fn new(path: impl Into<Vec<u8>>) -> Self {
            let cpath = CString::new(path).expect("failed to create CString");
            unsafe {
                let inner = RDKit_create_mol_supplier(cpath.as_ptr());
                Self(inner)
            }
        }

        /// reports whether or not `self` is at the end of the underlying file
        pub fn at_end(&self) -> bool {
            unsafe { RDKit_mol_supplier_at_end(self.0) }
        }

        pub fn next(&mut self) -> ROMol {
            unsafe { ROMol(RDKit_mol_supplier_next(self.0)) }
        }
    }

    impl Drop for SDMolSupplier {
        fn drop(&mut self) {
            unsafe {
                RDKit_delete_mol_supplier(self.0);
            }
        }
    }

    pub struct ROMol(*mut RDKit_ROMol);

    impl ROMol {
        pub fn to_smiles(&self) -> String {
            unsafe {
                let smiles = RDKit_MolToSmiles(self.0);
                let s = CStr::from_ptr(smiles);
                s.to_str().unwrap().to_owned()
            }
        }
    }

    impl Drop for ROMol {
        fn drop(&mut self) {
            unsafe {
                RDKit_ROMol_delete(self.0);
            }
        }
    }

    pub fn find_smarts_matches(mol: &ROMol, smarts: &str) -> Vec<Vec<usize>> {
        let mut len = 0;
        let mut match_size = 0;
        let smarts = CString::new(smarts).unwrap();
        unsafe {
            let matches = rdkit_sys::find_smarts_matches(
                mol.0,
                smarts.as_ptr(),
                &mut len,
                &mut match_size,
            );
            let matches = Vec::from_raw_parts(matches, len, len);

            let mut ret = Vec::new();
            for mat in matches.chunks(match_size) {
                ret.push(mat.into_iter().map(|&x| x as usize).collect());
            }
            ret
        }
    }
}

type TorsEnv = (usize, usize, usize, usize);

/// returns a vector of parameter ids matching `mol`. matching starts with the
/// first parameter and proceeds through the whole sequence of parameters, so
/// this should follow the SMIRNOFF typing rules
fn find_matches(params: &[(String, String)], mol: &ROMol) -> Vec<String> {
    let mut matches: HashMap<TorsEnv, String> = HashMap::new();
    for (id, smirks) in params {
        let env_matches = find_smarts_matches(mol, smirks);
        for mat in env_matches {
            assert_eq!(mat.len(), 4);
            let mat = (mat[0], mat[1], mat[2], mat[3]);
            matches.insert(mat, id.clone());
        }
    }
    let mut ret: Vec<_> = matches.into_values().collect();
    ret.sort();
    ret.dedup();
    ret
}

fn main() {
    let path = "/home/brent/omsf/chembl/chembl_33.sdf";
    let mut m = SDMolSupplier::new(path);
    // let mut out = File::create("out.smiles").unwrap();

    let forcefield = "openff-2.1.0.offxml";
    let ff = ForceField::load(forcefield).unwrap();
    let h = ff.get_parameter_handler("ProperTorsions").unwrap();
    let mut params = Vec::new();
    for p in h.parameters() {
        params.push((p.id(), p.smirks()));
    }

    let mut count = 0;
    while !m.at_end() && count < 1 {
        let mol = m.next();
        println!("{}", mol.to_smiles());
        let matches = find_matches(&params, &mol);
        dbg!(matches);
        count += 1;
    }
}

#[test]
fn first_molecule() {
    let path = "/home/brent/omsf/chembl/chembl_33.sdf";
    let mut m = SDMolSupplier::new(path);

    let forcefield = "openff-2.1.0.offxml";
    let ff = ForceField::load(forcefield).unwrap();
    let h = ff.get_parameter_handler("ProperTorsions").unwrap();
    let mut params = Vec::new();
    for p in h.parameters() {
        params.push((p.id(), p.smirks()));
    }

    let mol = m.next();

    let got_smiles = mol.to_smiles();
    let want_smiles = "Cc1cc(-c2csc(N=C(N)N)n2)cn1C";
    assert_eq!(got_smiles, want_smiles);

    let got = find_matches(&params, &mol);
    let want = vec![
        "t115", "t116", "t117", "t17", "t20", "t43", "t45", "t64", "t75",
        "t79", "t80", "t82", "t83", "t86",
    ];
    assert_eq!(got, want);
}
