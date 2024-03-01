use std::{collections::HashSet, fs::File, path::Path};

use openff_qcsubmit::results::{
    BaseResultCollection, OptimizationResultCollection,
    TorsionDriveResultCollection,
};
use rdkit_rs::ROMol;

fn load_inchis<T>(path: impl AsRef<Path>) -> HashSet<String>
where
    T: BaseResultCollection,
{
    T::parse_file(path)
        .unwrap()
        .entries()
        .into_values()
        .flatten()
        .map(|e| ROMol::from_smiles(&e.cmiles()).to_inchi_key())
        .collect()
}

fn main() {
    let base = Path::new("/home/brent/omsf/projects");
    let vf = base.join("valence-fitting/02_curate-data/datasets");

    let sage_tm_opt = vf.join("combined-opt.json");
    let sage_tm_td = vf.join("combined-td.json");
    let bench = base.join("benchmarking/datasets/industry.json");

    let mut total = load_inchis::<OptimizationResultCollection>(sage_tm_opt);
    total.extend(load_inchis::<TorsionDriveResultCollection>(sage_tm_td));
    total.extend(load_inchis::<OptimizationResultCollection>(bench));

    use std::io::Write;
    let mut out = File::create("inchis.dat").unwrap();
    for mol in total {
        writeln!(out, "{mol}").unwrap();
    }
}
