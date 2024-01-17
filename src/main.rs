use std::collections::{HashMap, HashSet};
use std::sync::atomic::AtomicUsize;

use clap::Parser;
use openff_toolkit::ForceField;
use rayon::iter::{ParallelBridge, ParallelIterator};
use rsearch::find_matches;
use rsearch::rdkit::{AromaticityModel, SanitizeFlags};
use rsearch::rdkit::{ROMol, SDMolSupplier};

fn load_want(path: &str) -> HashSet<String> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|s| s.trim().to_owned())
        .collect()
}

fn print_output(res: HashMap<String, Vec<String>>) {
    for (pid, moles) in res {
        println!("{pid}");
        for mol in moles {
            println!("\t{mol}");
        }
        println!();
    }
}

fn write_output(res: HashMap<String, Vec<String>>) {
    use std::io::Write;
    for (pid, moles) in res {
        let mut f = std::fs::File::create(format!("{pid}.smiles")).unwrap();
        for mol in moles {
            writeln!(f, "{mol}").unwrap();
        }
    }
}

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

    /// The path to the file listing the parameters to match against, one per
    /// line.
    #[arg(short, long, default_value = "want.params")]
    search_params: String,

    /// The number of threads to use. Defaults to the number of logical CPUs as
    /// detected by rayon.
    #[arg(short, long, default_value_t = 0)]
    threads: usize,

    /// Whether or not to create output files, one for each input parameter. If
    /// false, print the output to stdout.
    #[arg(short, long)]
    write_output: bool,
}

fn main() {
    let cli = Cli::parse();
    let m = SDMolSupplier::new(cli.molecule_file);
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

    let map_op = |mut mol: ROMol| -> Vec<(String, String)> {
        use SanitizeFlags as S;
        mol.sanitize(S::ALL ^ S::ADJUSTHS ^ S::SETAROMATICITY);
        mol.set_aromaticity(AromaticityModel::MDL);
        mol.assign_stereochemistry();
        mol.add_hs();

        let matches = find_matches(&params, &mol);

        let mut res: Vec<(String, String)> = Vec::new();
        let mut smiles = None;
        for pid in matches.intersection(&want) {
            if smiles.is_none() {
                smiles = Some(mol.to_smiles());
            }
            res.push((pid.to_string(), smiles.clone().unwrap()));
        }
        let cur = progress.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        if cur % 24_000 == 0 {
            eprintln!("{}% complete", cur / 24_000);
        }
        res
    };
    let results: Vec<_> = m.into_iter().par_bridge().flat_map(map_op).collect();

    let mut res: HashMap<String, Vec<String>> = HashMap::new();
    for (pid, mol) in results {
        res.entry(pid.to_string()).or_default().push(mol);
    }

    if cli.write_output {
        write_output(res);
    } else {
        print_output(res);
    }
}