use std::fs::read_to_string;

use clap::Parser;
use rsearch::{
    cluster::dbscan,
    matrix::Matrix,
    rdkit::{fingerprint::tanimoto, ROMol},
};

// default parameters from 2020-03-05-OpenFF-Training-Data-Selection

#[derive(Parser)]
struct Cli {
    /// The file of SMILES strings to read as input, one SMILES per line.
    smiles_file: String,

    /// DBSCAN parameter specifying the acceptable distance between a core point
    /// of a cluster and one of its neighbors.
    #[arg(short, long, default_value_t = 0.5)]
    epsilon: f64,

    /// DBSCAN parameter specifying the minimum number of points required to
    /// form a dense region.
    #[arg(short, long, default_value_t = 5)]
    min_pts: usize,

    /// Morgan fingerprinting radius
    #[arg(short, long, default_value_t = 4)]
    radius: u32,
}

fn main() {
    let cli = Cli::parse();

    let fps: Vec<_> = read_to_string(cli.smiles_file)
        .unwrap()
        .lines()
        .map(|smiles| {
            let mut mol = ROMol::from_smiles(smiles);
            mol.openff_clean();
            mol.morgan_fingerprint_bit_vec::<1024>(cli.radius)
        })
        .collect();

    let n = fps.len();
    let mut db = Matrix::zeros(n, n);
    // computing 1 - tanimoto here because dbscan groups items with _low_
    // distance, rather than high similarity
    for i in 0..n {
        db[(i, i)] = 0.0;
        for j in 0..i {
            let t = tanimoto(&fps[i], &fps[j]);
            db[(i, j)] = 1.0 - t;
            db[(j, i)] = 1.0 - t;
        }
    }

    let labels = dbscan(&db, cli.epsilon, cli.min_pts);

    // TODO some kind of useful output
    for group in labels.chunks(20) {
        for g in group {
            print!("{g:?}");
        }
        println!();
    }
}
