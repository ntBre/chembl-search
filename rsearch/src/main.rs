use std::fs::File;
use std::io::Write;

use openff_toolkit::ForceField;
use rdkit::{find_smarts_matches, SDMolSupplier};

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
        pub fn to_smiles(self) -> String {
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

fn main() {
    let path = "/home/brent/omsf/chembl/chembl_33.sdf";
    let mut m = SDMolSupplier::new(path);
    let mut out = File::create("out.smiles").unwrap();
    let smarts = "[#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]";

    let forcefield = "openff-2.1.0.offxml";
    let ff = ForceField::load(forcefield).unwrap();
    let h = ff.get_parameter_handler("ProperTorsions").unwrap();
    let mut pairs = Vec::new();
    for p in h.parameters() {
        pairs.push((p.id(), p.smirks()));
    }

    let mut count = 0;
    while !m.at_end() && count < 50 {
        let mol = m.next();
        let matches = find_smarts_matches(&mol, smarts);
        for (i, mat) in matches.iter().enumerate() {
            println!("Match {i} : {} {} {} {}", mat[0], mat[1], mat[2], mat[3]);
        }
        let smiles = mol.to_smiles();
        writeln!(out, "{smiles}").unwrap();
        count += 1;
    }
}
