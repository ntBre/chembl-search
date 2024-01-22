use std::{
    collections::{HashMap, HashSet},
    fs::read_to_string,
    path::Path,
};

use clap::Parser;
use openff_toolkit::ForceField;
use rsearch::{
    cluster::{dbscan, Label},
    matrix::Matrix,
    rdkit::{find_smarts_matches, fingerprint::tanimoto, ROMol},
};

struct Report<'a> {
    max: usize,
    nfps: usize,
    noise: usize,
    clusters: Vec<Vec<usize>>,
    smiles: Vec<&'a str>,
    mols: Vec<ROMol>,
    ff: String,
    param: Option<String>,
}

impl Report<'_> {
    fn generate(&self, path: impl AsRef<Path>) -> std::io::Result<()> {
        use std::io::Write;
        let mut out = std::fs::File::create(path).unwrap();
        writeln!(out, "<html>")?;
        writeln!(
            out,
            "{nfps} molecules, {max} clusters, {noise} noise points, \
        pruned {} empty clusters",
            self.max + 1 - self.clusters.len(),
            nfps = self.nfps,
            max = self.max,
            noise = self.noise
        )?;
        // TODO allow other handlers, another cli option...
        let map: HashMap<_, _> = ForceField::load(&self.ff)
            .unwrap()
            .get_parameter_handler("ProperTorsions")
            .unwrap()
            .parameters()
            .into_iter()
            .map(|p| (p.id(), p.smirks()))
            .collect();
        for (i, c) in self.clusters.iter().enumerate() {
            let smile = self.smiles[c[0]];
            let mol = &self.mols[c[0]];

            let mut hl_atoms = Vec::new();
            if let Some(pid) = &self.param {
                if let Some(smirks) = map.get(pid) {
                    let tmp = find_smarts_matches(mol, &smirks);
                    if !tmp.is_empty() {
                        hl_atoms = tmp[0].clone();
                    }
                }
            }

            let svg = mol.draw_svg(400, 300, "", &hl_atoms);
            writeln!(out, "<h1>Cluster {i}</h1>")?;
            writeln!(out, "<p>Molecule 1/{}</p>", c.len())?;
            writeln!(out, "<p>{} atoms</p>", mol.num_atoms())?;
            writeln!(out, "<p>SMILES: {smile}</p>")?;
            writeln!(out, "{svg}")?;
        }
        writeln!(out, "</html>")?;
        Ok(())
    }
}

/// The default DBSCAN parameters are taken from the
/// 2020-03-05-OpenFF-Training-Data-Selection/select_TrainingDS.ipynb in the
/// qca-dataset-submission repo commit 79ee3a3
#[derive(Parser)]
struct Cli {
    /// The file of SMILES strings to read as input, one SMILES per line.
    smiles_file: String,

    /// The force field to use for parameter labeling.
    #[arg(short, long, default_value = "openff-2.1.0.offxml")]
    forcefield: String,

    /// The parameter to use when highlighting atoms in the molecules.
    #[arg(short, long)]
    parameter: Option<String>,

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

    let max = *labels
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
    let mut clusters: Vec<Vec<usize>> = vec![vec![]; max + 1];

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
    const MAX_ATOMS: usize = 154;

    let existing_inchis: HashSet<_> = read_to_string("inchis.dat")
        .unwrap()
        .split_ascii_whitespace()
        .map(|s| s.to_owned())
        .collect();
    let new_inchis: Vec<_> = mols.iter().map(ROMol::to_inchi_key).collect();

    // filter out molecules larger than MAX_ATOMS or with inchi_keys already
    // covered by our existing data sets; then filter any empty clusters
    clusters.iter_mut().for_each(|cluster| {
        cluster.retain(|mol_idx| {
            mols[*mol_idx].num_atoms() <= MAX_ATOMS
                && !existing_inchis.contains(&new_inchis[*mol_idx])
        });
    });
    clusters.retain(|c| !c.is_empty());

    // TODO highlight involved atoms - this will require additional printing in
    // the original search, maybe tsv of smiles and involved atoms

    Report {
        max,
        nfps,
        noise,
        clusters,
        smiles,
        mols,
        ff: cli.forcefield,
        param: cli.parameter,
    }
    .generate("prints.html")?;

    Ok(())
}
