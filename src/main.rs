use std::collections::HashMap;
use std::path::Path;
use std::sync::atomic::AtomicUsize;

use clap::Parser;
use log::trace;
use openff_toolkit::ForceField;
use rayon::iter::{ParallelBridge, ParallelIterator};
use rdkit_rs::{RDError, ROMol, SDMolSupplier};
use rsearch::{find_matches, load_want, print_output, write_output};

#[derive(Parser)]
struct Cli {
    /// The OpenFF force field to load from the toolkit.
    #[arg(short, long, default_value = "openff-2.1.0.offxml")]
    forcefield: String,

    /// The `Parameter` type for which to extract parameters. Allowed options
    /// are valid arguments to `ForceField.get_parameter_handler`, such as
    /// Bonds, Angles, or ProperTorsions.
    #[arg(short, long, default_value = "ProperTorsions")]
    parameter_type: String,

    /// The path to the SDF file from which to read Molecules.
    #[arg(short, long, default_value = "chembl_33.sdf")]
    molecule_file: String,

    /// The path to the file listing the parameter identifiers to match against,
    /// one per line. These must correspond to parameters in the provided
    /// forcefield.
    #[arg(short, long, default_value = "want.params")]
    search_params: String,

    /// The number of threads to use. Defaults to the number of logical CPUs as
    /// detected by rayon.
    #[arg(short, long, default_value_t = 0)]
    threads: usize,

    /// Where to write output SMILES files, one for each input parameter. If
    /// false, print the output to stdout.
    #[arg(short, long)]
    output_dir: Option<String>,
}

fn main() {
    env_logger::init();

    let cli = Cli::parse();
    trace!("initializing mol supplier from {}", cli.molecule_file);
    let m = SDMolSupplier::new(cli.molecule_file).unwrap();
    trace!("initializing force field from {}", cli.forcefield);
    let ff = ForceField::load(&cli.forcefield).unwrap();
    let h = ff.get_parameter_handler(&cli.parameter_type).unwrap();
    let mut params = Vec::new();
    for p in h.parameters() {
        params.push((p.id(), ROMol::from_smarts(&p.smirks())));
    }

    let want = load_want(&cli.search_params);

    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .unwrap();

    let progress = AtomicUsize::new(0);

    let map_op = |mol: Result<ROMol, RDError>| -> Vec<(String, String)> {
        let cur = progress.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        if cur % 24_000 == 0 {
            eprintln!("{}% complete", cur / 24_000);
        }
        let Ok(mut mol) = mol else {
            return Vec::new();
        };
        trace!("calling clean");
        mol.openff_clean();

        trace!("calling find_matches");
        let matches = find_matches(&params, &mol);

        let mut res: Vec<(String, String)> = Vec::new();
        let mut smiles = None;
        for pid in matches.intersection(&want) {
            if smiles.is_none() {
                smiles = Some(mol.to_smiles());
            }
            res.push((pid.to_string(), smiles.clone().unwrap()));
        }
        res
    };
    let results: Vec<_> = m.into_iter().par_bridge().flat_map(map_op).collect();

    let mut res: HashMap<String, Vec<String>> = HashMap::new();
    for (pid, mol) in results {
        res.entry(pid.to_string()).or_default().push(mol);
    }

    if let Some(dir) = cli.output_dir {
        let dir = Path::new(&dir);
        if !dir.exists() && std::fs::create_dir_all(dir).is_err() {
            eprintln!("failed to create output dir {dir:?}");
            eprintln!("falling back to stdout");
            print_output(res);
            return;
        }
        write_output(dir, res);
    } else {
        print_output(res);
    }
}
