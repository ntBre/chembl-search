//! Like the main rsearch program, but instead of searching by parameters in a
//! force field, simply read a sequence of SMARTS patterns and return molecules
//! matching those.

use std::{collections::HashMap, path::Path, sync::atomic::AtomicUsize};

use clap::Parser;
use rayon::iter::{ParallelBridge, ParallelIterator};
use rsearch::{
    load_want, print_output,
    rdkit::{find_smarts_matches_mol, RDError, ROMol, SDMolSupplier},
    write_output,
};

#[derive(Parser)]
struct Cli {
    /// The path to the SDF file from which to read Molecules.
    #[arg(short, long, default_value = "chembl_33.sdf")]
    molecule_file: String,

    /// The path to the file listing the parameter SMIRKS to match against, one
    /// per line, followed by a label.
    #[arg(short, long, default_value = "want.smirks")]
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
    let cli = Cli::parse();
    let m = SDMolSupplier::new(cli.molecule_file);
    let params: Vec<_> = load_want(&cli.search_params)
        .into_iter()
        .map(|s| {
            let sp: Vec<_> = s.split_ascii_whitespace().collect();
            assert_eq!(sp.len(), 2);
            let smirks = sp[0];
            let label = sp[1];
            (label.to_owned(), ROMol::from_smarts(&smirks))
        })
        .collect();

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
        mol.openff_clean();

        let mut smiles = None; /* cache for to_smiles */
        let mut res: Vec<(String, String)> = Vec::new();
        for (pid, smirks) in &params {
            if !find_smarts_matches_mol(&mol, smirks).is_empty() {
                if smiles.is_none() {
                    smiles = Some(mol.to_smiles());
                }
                res.push((pid.to_string(), smiles.clone().unwrap()));
            }
        }
        res
    };
    let results: Vec<_> = m.into_iter().par_bridge().flat_map(map_op).collect();

    let mut res: HashMap<String, Vec<String>> = HashMap::new();
    for (label, _) in params {
        res.insert(label, Vec::new());
    }
    for (pid, mol) in results {
        res.entry(pid.to_string()).and_modify(|mols| mols.push(mol));
    }

    println!("Summary:");
    for (pid, mols) in &res {
        println!("{pid} {}", mols.len());
    }
    println!("-----------------------------------");

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
