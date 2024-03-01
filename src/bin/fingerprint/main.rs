use std::{
    collections::{HashMap, HashSet},
    fs::read_to_string,
    io::{self},
    path::Path,
    process::exit,
};

use config::Parameter;
use log::info;
use openff_toolkit::ForceField;
use rdkit_rs::{fingerprint::tanimoto, ROMol};
use rsearch::{
    cluster::{dbscan, Label},
    utils::{load_mols, make_fps, Report},
};

type Pid = String;
type Smirks = String;

mod config;

fn single(
    ff: &str,
    max_atoms: usize,
    radius: u32,
    param: Parameter,
) -> Result<(), io::Error> {
    info!("starting {}", param.id);

    let s = read_to_string(&param.smiles)
        .unwrap_or_else(|e| panic!("failed to read {} for {e}", param.smiles));
    let smiles: Vec<_> = s.lines().collect();

    info!("{} initial smiles", smiles.len());

    let map: HashMap<Pid, Smirks> = ForceField::load(ff)
        .unwrap()
        .get_parameter_handler(&param.typ)
        .unwrap()
        .parameters()
        .into_iter()
        .map(|p| (p.id(), p.smirks()))
        .collect();

    let mol_map: Vec<(Pid, ROMol)> = ForceField::load(ff)
        .unwrap()
        .get_parameter_handler(&param.typ)
        .unwrap()
        .parameters()
        .into_iter()
        .map(|p| (p.id(), ROMol::from_smarts(&p.smirks())))
        .collect();

    let existing_inchis: HashSet<_> = read_to_string("data/inchis.dat")
        .unwrap()
        .split_ascii_whitespace()
        .map(|s| s.to_owned())
        .collect();

    let mols = load_mols(
        smiles,
        max_atoms,
        param.fragment,
        &param.id,
        existing_inchis,
        &mol_map,
    );

    info!("making fingerprints");

    let fps: Vec<_> = make_fps(&mols, radius);
    let nfps = fps.len();
    let distance_fn = |i, j| {
        if i == j {
            0.0
        } else {
            1.0 - tanimoto(&fps[i], &fps[j])
        }
    };

    info!("running DBSCAN");

    let labels = dbscan(
        nfps,
        nfps,
        distance_fn,
        param.dbscan.epsilon,
        param.dbscan.min_pts,
    );

    let max = match labels
        .iter()
        .filter_map(|l| match l {
            Label::Cluster(n) => Some(n),
            _ => None,
        })
        .max()
    {
        Some(n) => *n,
        None => {
            dbg!(labels);
            eprintln!("error: all noise points, exiting");
            exit(1);
        }
    };

    let mut clusters: Vec<Vec<usize>> = vec![vec![]; max + 1];
    let mut noise = 0;
    for (i, l) in labels.iter().enumerate() {
        match l {
            Label::Cluster(n) => clusters[*n].push(i),
            _ => noise += 1,
        }
    }

    info!("generating report");

    let output = Path::new(&param.smiles).with_extension("html");
    Report {
        args: std::env::args().collect::<Vec<_>>(),
        max,
        nfps,
        noise,
        clusters,
        mols,
        parameter: &param.id,
        map,
        mol_map,
    }
    .generate(output)?;
    Ok(())
}

fn main() -> io::Result<()> {
    env_logger::init();

    let args: Vec<_> = std::env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: fingerprint <config.toml>");
        std::process::exit(1);
    }

    let cli = config::Config::load(&args[1]);

    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .unwrap();

    for param in cli.parameters {
        single(&cli.forcefield, cli.max_atoms, cli.radius, param)?;
    }

    Ok(())
}
