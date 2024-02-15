use std::fs::read_to_string;
use std::path::Path;

use serde::Deserialize;

#[derive(Deserialize)]
pub(crate) struct Dbscan {
    /// The maximum acceptable distance between a core point of a cluster and
    /// one of its neighbors.
    pub(crate) epsilon: f64,

    /// The minimum number of points required to form a dense region
    pub(crate) min_pts: usize,
}

impl Default for Dbscan {
    /// The default parameters are taken from the
    /// 2020-03-05-OpenFF-Training-Data-Selection/select_TrainingDS.ipynb in the
    /// qca-dataset-submission repo commit 79ee3a3
    fn default() -> Self {
        Self {
            epsilon: 0.5,
            min_pts: 1,
        }
    }
}

// impl Default for Config {
//     fn default() -> Self {
//         Self {
//             smiles_file: String::new(),
//             max_atoms: 80,
//             forcefield: "input/tm.v2.offxml".to_owned(),
//             parameter: None,
//             parameter_type: "ProperTorsions".to_owned(),
//             dbscan: DBSCAN::default(),
//             radius: 4,
//             threads: 0,
//             fragment: false,
//         }
//     }
// }

#[derive(Deserialize)]
pub(crate) struct Config {
    /// The file of SMILES strings to read as input, one SMILES per line.
    pub(crate) smiles_file: String,

    /// The maximum number of atoms to consider
    pub(crate) max_atoms: usize,

    /// The force field to use for parameter labeling.
    pub(crate) forcefield: String,

    /// The parameter to use when highlighting atoms in the molecules.
    pub(crate) parameter: String,

    /// The `Parameter` type for which to extract parameters. Allowed options
    /// are valid arguments to `ForceField.get_parameter_handler`, such as
    /// Bonds, Angles, or ProperTorsions.
    pub(crate) parameter_type: String,

    /// [DBSCAN] parameters
    pub(crate) dbscan: Dbscan,

    /// Morgan fingerprinting radius
    pub(crate) radius: u32,

    /// The number of threads to use. Defaults to the number of logical CPUs as
    /// detected by rayon.
    pub(crate) threads: usize,

    /// Whether or not to fragment the molecules before the fingerprinting
    /// analysis.
    pub(crate) fragment: bool,
}

impl Config {
    pub(crate) fn load(path: impl AsRef<Path>) -> Self {
        toml::from_str(&read_to_string(path).unwrap()).unwrap()
    }
}
