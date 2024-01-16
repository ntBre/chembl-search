use std::io::Write;
use std::{
    ffi::{CStr, CString},
    fs::File,
};

use openff_toolkit::ForceField;
use rdkit::SDMolSupplier;
use rdkit_sys::{RDKit_MolToSmiles, RDKit_ROMol_delete};

/// TODO move this to its own crate, possibly in a workspace with rdkit-sys
pub mod rdkit {
    use std::ffi::CString;

    use rdkit_sys::{
        RDKit_ROMol, RDKit_SDMolSupplier, RDKit_create_mol_supplier,
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

        pub fn next(&mut self) -> *mut RDKit_ROMol {
            unsafe { RDKit_mol_supplier_next(self.0) }
        }
    }

    impl Drop for SDMolSupplier {
        fn drop(&mut self) {
            unsafe {
                RDKit_delete_mol_supplier(self.0);
            }
        }
    }
}

fn main() {
    let path = "/home/brent/omsf/chembl/chembl_33.sdf";
    let mut m = SDMolSupplier::new(path);
    let mut out = File::create("out.smiles").unwrap();
    let smarts = CString::new("[#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]").unwrap();

    let forcefield = "openff-2.1.0.offxml";
    let ff = ForceField::load(forcefield).unwrap();
    let h = ff.get_parameter_handler("ProperTorsions").unwrap();
    let mut pairs = Vec::new();
    for p in h.parameters() {
        pairs.push((p.id(), p.smirks()));
    }

    unsafe {
        let mut count = 0;
        while !m.at_end() && count < 50 {
            let mol = m.next();
            let mut len = 0;
            let mut match_size = 0;
            let matches = rdkit_sys::find_smarts_matches(
                mol,
                smarts.as_ptr(),
                &mut len,
                &mut match_size,
            );
            let matches = Vec::from_raw_parts(matches, len, len);
            for (i, mat) in matches.chunks(match_size).enumerate() {
                println!(
                    "Match {i} : {} {} {} {}",
                    mat[0], mat[1], mat[2], mat[3]
                );
            }
            let smiles = RDKit_MolToSmiles(mol);
            let s = CStr::from_ptr(smiles);
            writeln!(out, "{}", s.to_str().unwrap()).unwrap();
            count += 1;
            RDKit_ROMol_delete(mol);
        }
    }
}
