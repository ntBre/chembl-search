use std::{
    collections::{HashMap, HashSet},
    sync::Mutex,
};

use matrix::Matrix;
use rayon::iter::{ParallelBridge, ParallelIterator};
use rdkit::{
    bitvector::BitVector, find_smarts_matches_mol, fingerprint::tanimoto, ROMol,
};

/// TODO move this to its own crate, possibly in a workspace with rdkit-sys
pub mod rdkit;

pub mod cluster;
pub mod matrix;

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

pub fn distance_matrix(nfps: usize, fps: Vec<BitVector>) -> Matrix<f64> {
    let db = Mutex::new(Matrix::zeros(nfps, nfps));

    eprintln!("allocated {} elements", nfps * nfps);

    // combinations of indices under the diagonal: 0..N, 0..i
    let combos = (0..nfps)
        .into_iter()
        .flat_map(|i| (0..i).into_iter().map(move |j| (i, j)));

    let combos = combos
        .par_bridge()
        .map(|(i, j)| (i, j, tanimoto(&fps[i], &fps[j])));

    // computing 1 - tanimoto here because dbscan groups items with _low_
    // distance, rather than high similarity
    combos.for_each(|(i, j, t)| {
        let mut db = db.lock().unwrap();
        db[(i, j)] = 1.0 - t;
        db[(j, i)] = 1.0 - t;
    });
    db.into_inner().unwrap()
}

#[cfg(test)]
mod tests {
    use openff_toolkit::ForceField;

    use crate::{
        find_matches,
        rdkit::{
            find_smarts_matches, find_smarts_matches_mol, ROMol, SDMolSupplier,
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
        mol.openff_clean();

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
        mol.openff_clean();
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
