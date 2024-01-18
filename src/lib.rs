use std::collections::{HashMap, HashSet};

use rdkit::{find_smarts_matches_mol, ROMol};

/// TODO move this to its own crate, possibly in a workspace with rdkit-sys
pub mod rdkit {
    use std::{
        collections::HashMap,
        ffi::{c_uint, CStr, CString},
    };

    use bitflags::bitflags;
    use rdkit_sys::{
        RDKit_MolToSmiles, RDKit_ROMol, RDKit_ROMol_delete,
        RDKit_SDMolSupplier, RDKit_SmartsToMol, RDKit_SmilesToMol,
        RDKit_create_mol_supplier, RDKit_delete_mol_supplier,
        RDKit_mol_supplier_at_end, RDKit_mol_supplier_next,
    };

    pub struct SDMolSupplier(*mut RDKit_SDMolSupplier);

    unsafe impl Send for SDMolSupplier {}

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
    }

    impl Drop for SDMolSupplier {
        fn drop(&mut self) {
            unsafe {
                RDKit_delete_mol_supplier(self.0);
            }
        }
    }

    impl Iterator for SDMolSupplier {
        type Item = ROMol;

        fn next(&mut self) -> Option<Self::Item> {
            if self.at_end() {
                return None;
            }
            Some(unsafe { ROMol(RDKit_mol_supplier_next(self.0)) })
        }
    }

    pub struct ROMol(*mut RDKit_ROMol);

    impl ROMol {
        pub fn from_smiles(smiles: &str) -> Self {
            let s = CString::new(smiles).expect("failed to create CString");
            unsafe { Self(RDKit_SmilesToMol(s.as_ptr())) }
        }

        pub fn from_smarts(smarts: &str) -> Self {
            let s = CString::new(smarts).expect("failed to create CString");
            unsafe { Self(RDKit_SmartsToMol(s.as_ptr())) }
        }

        pub fn to_smiles(&self) -> String {
            unsafe {
                let smiles = RDKit_MolToSmiles(self.0);
                let s = CStr::from_ptr(smiles);
                s.to_str().unwrap().to_owned()
            }
        }

        pub fn sanitize(&mut self, ops: SanitizeFlags) {
            unsafe {
                rdkit_sys::RDKit_SanitizeMol(self.0, ops.bits());
            }
        }

        pub fn set_aromaticity(&mut self, mdl: AromaticityModel) {
            unsafe {
                rdkit_sys::RDKit_SetAromaticity(self.0, mdl.bits());
            }
        }

        pub fn assign_stereochemistry(&mut self) {
            unsafe {
                rdkit_sys::RDKit_AssignStereochemistry(self.0);
            }
        }

        pub fn add_hs(&mut self) {
            unsafe {
                rdkit_sys::RDKit_AddHs(self.0);
            }
        }

        pub fn morgan_fingerprint(
            &self,
            radius: c_uint,
        ) -> HashMap<usize, usize> {
            unsafe {
                let mut len = 0;
                let ptr = rdkit_sys::RDKit_MorganFingerprint(
                    self.0, radius, &mut len,
                );
                let v1 = Vec::from_raw_parts(ptr, len, len);
                v1.into_iter()
                    .map(|v| (v.bit as usize, v.count as usize))
                    .collect()
            }
        }

        pub fn morgan_fingerprint_bit_vec<const N: usize>(
            &self,
            radius: c_uint,
        ) -> [usize; N] {
            unsafe {
                let mut ret = [0; N];
                rdkit_sys::RDKit_MorganFingerprintBitVector(
                    self.0,
                    radius,
                    N,
                    ret.as_mut_ptr(),
                );
                ret
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

    unsafe impl Send for ROMol {}
    unsafe impl Sync for ROMol {}

    bitflags! {
        pub struct SanitizeFlags: c_uint {
            const NONE =                    0x0;
            const CLEANUP =                 0x1;
            const PROPERTIES =              0x2;
            const SYMMRINGS =               0x4;
            const KEKULIZE =                0x8;
            const FINDRADICALS =            0x10;
            const SETAROMATICITY =          0x20;
            const SETCONJUGATION =          0x40;
            const SETHYBRIDIZATION =        0x80;
            const CLEANUPCHIRALITY =        0x100;
            const ADJUSTHS =                0x200;
            const CLEANUP_ORGANOMETALLICS = 0x400;
            const ALL =                     0xFFFFFFF;
        }
    }

    bitflags! {
        pub struct AromaticityModel: c_uint {
            const DEFAULT = 0x0;
            const RDKIT = 0x1;
            const SIMPLE = 0x2;
            const MDL = 0x4;
            const CUSTOM = 0xFFFFFFF;
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
                ret.push(mat.iter().map(|&x| x as usize).collect());
            }
            ret
        }
    }

    pub fn find_smarts_matches_mol(
        mol: &ROMol,
        smarts: &ROMol,
    ) -> Vec<Vec<usize>> {
        let mut len = 0;
        let mut match_size = 0;
        unsafe {
            let matches = rdkit_sys::find_smarts_matches_mol(
                mol.0,
                smarts.0,
                &mut len,
                &mut match_size,
            );
            let matches = Vec::from_raw_parts(matches, len, len);

            let mut ret = Vec::new();
            for mat in matches.chunks(match_size) {
                ret.push(mat.iter().map(|&x| x as usize).collect());
            }
            ret
        }
    }

    pub mod fingerprint {
        /// print `bv` to stdout in 16 groups of 4 per row
        pub fn print_bit_vec(bv: &[usize]) {
            for line in bv.chunks(16 * 4) {
                for chunk in line.chunks(4) {
                    print!(" ");
                    for elt in chunk {
                        print!("{elt}");
                    }
                }
                println!();
            }
        }

        /// return the number of bits in the intersection of a and b
        pub fn intersect(a: &[usize], b: &[usize]) -> usize {
            let mut ret = 0;
            for (a, b) in a.iter().zip(b) {
                ret += a & b;
            }
            ret
        }

        /// Computes the Tanimoto distance between bit vectors a and b
        ///
        /// T(a, b) = (a ∩ b) / (a + b - a ∩ b), at least according to
        /// featurebase.com/blog/tanimoto-and-chemical-similarity-in-featurebase
        pub fn tanimoto(a: &[usize], b: &[usize]) -> f64 {
            let num = intersect(a, b);
            let den: usize =
                a.iter().sum::<usize>() + b.iter().sum::<usize>() - num;
            num as f64 / den as f64
        }
    }
}

/// returns the set of parameter ids matching `mol`. matching starts with the
/// first parameter and proceeds through the whole sequence of parameters, so
/// this should follow the SMIRNOFF typing rules
pub fn find_matches(
    params: &[(String, ROMol)],
    mol: &ROMol,
) -> HashSet<String> {
    let mut matches = HashMap::new();
    for (id, smirks) in params {
        let env_matches = find_smarts_matches_mol(mol, smirks);
        for mat in env_matches {
            matches.insert(mat, id.clone());
        }
    }
    matches.into_values().collect()
}

#[cfg(test)]
mod tests {
    use openff_toolkit::ForceField;

    use crate::{
        find_matches,
        rdkit::{
            find_smarts_matches, find_smarts_matches_mol, AromaticityModel,
            ROMol, SDMolSupplier, SanitizeFlags,
        },
    };

    #[test]
    fn first_molecule() {
        let path = "/home/brent/omsf/chembl/chembl_33.sdf";
        let mut m = SDMolSupplier::new(path);

        let forcefield = "openff-2.1.0.offxml";
        let ff = ForceField::load(forcefield).unwrap();
        let h = ff.get_parameter_handler("ProperTorsions").unwrap();
        let mut params = Vec::new();
        for p in h.parameters() {
            params.push((p.id(), ROMol::from_smarts(&p.smirks())));
        }

        let mut mol = m.next().unwrap();

        mol.sanitize(
            SanitizeFlags::ALL
                ^ SanitizeFlags::ADJUSTHS
                ^ SanitizeFlags::SETAROMATICITY,
        );
        mol.set_aromaticity(AromaticityModel::MDL);
        mol.assign_stereochemistry();
        mol.add_hs();

        let got = find_matches(&params, &mol);
        let want = vec![
            "t115", "t116", "t117", "t17", "t20", "t43", "t45", "t64", "t75",
            "t79", "t80", "t82", "t83", "t86",
        ];
        assert_eq!(got, want.into_iter().map(|s| s.to_owned()).collect());
    }

    #[test]
    fn single_match() {
        let mut mol = ROMol::from_smiles("Cc1cc(-c2csc(N=C(N)N)n2)cn1C");
        mol.sanitize(
            SanitizeFlags::ALL
                ^ SanitizeFlags::ADJUSTHS
                ^ SanitizeFlags::SETAROMATICITY,
        );

        mol.set_aromaticity(AromaticityModel::MDL);
        mol.assign_stereochemistry();
        mol.add_hs();
        let smarts = "[*:1]-[#16X2,#16X3+1:2]-[#6:3]~[*:4]"; // t115
        let got = find_smarts_matches(&mol, smarts);
        let want = vec![
            vec![5, 6, 7, 8],
            vec![5, 6, 7, 12],
            vec![7, 6, 5, 4],
            vec![7, 6, 5, 20],
        ];
        assert_eq!(got, want);

        let smarts = ROMol::from_smarts("[*:1]-[#16X2,#16X3+1:2]-[#6:3]~[*:4]");
        let got = find_smarts_matches_mol(&mol, &smarts);
        assert_eq!(got, want);
    }
}
