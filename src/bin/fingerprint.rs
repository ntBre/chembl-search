use std::fs::read_to_string;

use clap::Parser;
use rsearch::{
    cluster::{dbscan, Label},
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

fn main() -> std::io::Result<()> {
    let cli = Cli::parse();

    let s = read_to_string(cli.smiles_file).unwrap();
    let smiles: Vec<_> = s.lines().collect();
    let mols: Vec<_> = smiles
        .iter()
        .map(|smiles| {
            let mut mol = ROMol::from_smiles(smiles);
            mol.openff_clean();
            mol
        })
        .collect();
    let fps: Vec<_> = mols
        .iter()
        .map(|mol| mol.morgan_fingerprint_bit_vec::<1024>(cli.radius))
        .collect();

    let nfps = fps.len();
    let mut db = Matrix::zeros(nfps, nfps);
    // computing 1 - tanimoto here because dbscan groups items with _low_
    // distance, rather than high similarity
    for i in 0..nfps {
        db[(i, i)] = 0.0;
        for j in 0..i {
            let t = tanimoto(&fps[i], &fps[j]);
            db[(i, j)] = 1.0 - t;
            db[(j, i)] = 1.0 - t;
        }
    }

    let labels = dbscan(&db, cli.epsilon, cli.min_pts);

    let max = labels
        .iter()
        .filter_map(|l| match l {
            Label::Cluster(n) => Some(n),
            _ => None,
        })
        .max()
        .unwrap();

    // each entry contains a vec of molecule indices (smiles line numbers)
    // corresponding to that cluster. clusters[i] is the ith cluster with
    // members clusters[i][0..n], where n is however many members there are
    let mut clusters: Vec<Vec<usize>> = vec![vec![]; *max + 1];

    let mut noise = 0;
    for (i, l) in labels.iter().enumerate() {
        match l {
            Label::Cluster(n) => clusters[*n].push(i),
            _ => noise += 1,
        }
    }

    // largest molecule in current training and benchmark sets:
    // sage-tm opt: 76
    // sage-tm td: 67
    // bench: 154

    // TODO consider molecule size, not just centroid. some of these centroids
    // are way too big

    // TODO highlight involved atoms - this will require additional printing in
    // the original search, maybe tsv of smiles and involved atoms

    use std::io::Write;
    let mut out = std::fs::File::create("prints.html").unwrap();
    writeln!(out, "<html>")?;

    writeln!(
        out,
        "{nfps} molecules, {max} clusters, {noise} noise points"
    )?;

    for (i, c) in clusters.iter().enumerate() {
        let smile = smiles[c[0]];
        let mol = &mols[c[0]];
        let svg = mol.draw_svg(400, 300, "", &[]);
        writeln!(out, "<h1>Cluster {i}</h1>")?;
        writeln!(out, "<p>Molecule 0/{}</p>", c.len())?;
        writeln!(out, "<p>{} atoms</p>", mol.num_atoms())?;
        writeln!(out, "<p>SMILES: {smile}</p>")?;
        writeln!(out, "{svg}")?;
    }
    writeln!(out, "</html>")?;

    Ok(())
}
