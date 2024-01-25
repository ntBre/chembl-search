use std::{
    collections::{HashMap, HashSet},
    fs::read_to_string,
    io::{self, Write},
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
    mols: Vec<ROMol>,
    cli: &'a Cli,
}

type Pid = String;
type Smirks = String;

impl Report<'_> {
    fn generate(&self, path: impl AsRef<Path>) -> io::Result<()> {
        let mut out = std::fs::File::create(path).unwrap();
        writeln!(out, "<html>")?;
        writeln!(
            out,
            "{nfps} molecules, {max} clusters, {noise} noise points, \
        pruned {} empty clusters",
            self.max + 1 - self.clusters.len(),
            nfps = self.nfps,
            max = self.max + 1,
            noise = self.noise
        )?;

        let map: HashMap<Pid, Smirks> = ForceField::load(&self.cli.forcefield)
            .unwrap()
            .get_parameter_handler(&self.cli.parameter_type)
            .unwrap()
            .parameters()
            .into_iter()
            .map(|p| (p.id(), p.smirks()))
            .collect();

        for (i, c) in self.clusters.iter().enumerate() {
            for mol_idx in c {
                self.add_svg(&mut out, "mol: ", &map, *mol_idx)?;
            }
            writeln!(out, "<h1>Cluster {}, {} molecules</h1>", i + 1, c.len())?;

            self.add_svg(&mut out, "Central Molecule", &map, c[0])?;

            // find the index of the smallest molecule
            let (idx, _) = c
                .iter()
                .enumerate()
                .map(|(i, m)| (i, self.mols[*m].num_atoms()))
                .min_by_key(|x| x.1)
                .unwrap();
            if idx != 0 {
                self.add_svg(&mut out, "Smallest Molecule", &map, c[idx])?;
            }
        }
        writeln!(out, "</html>")?;
        Ok(())
    }

    fn add_svg(
        &self,
        out: &mut impl Write,
        msg: &str,
        map: &HashMap<String, String>,
        idx: usize,
    ) -> io::Result<()> {
        let mol = &self.mols[idx];
        let smile = mol.to_smiles();
        let svg = self.make_svg(map, mol);
        writeln!(out, "<p>{msg}</p>")?;
        writeln!(out, "<p>{} atoms</p>", mol.num_atoms())?;
        writeln!(out, "<p>SMILES: {smile}</p>")?;
        writeln!(out, "{svg}")?;
        Ok(())
    }

    fn make_svg(&self, map: &HashMap<String, String>, mol: &ROMol) -> String {
        let mut hl_atoms = Vec::new();
        if let Some(pid) = &self.cli.parameter {
            if let Some(smirks) = map.get(pid) {
                eprintln!("{}", mol.to_smiles());
                let tmp = find_smarts_matches(mol, smirks);
                if !tmp.is_empty() {
                    hl_atoms = tmp[0].clone();
                } else {
                    // panic!("smirks doesn't match any more");
                }
            }
        }
        mol.draw_svg(400, 300, "", &hl_atoms)
    }
}

/// The default DBSCAN parameters are taken from the
/// 2020-03-05-OpenFF-Training-Data-Selection/select_TrainingDS.ipynb in the
/// qca-dataset-submission repo commit 79ee3a3
#[derive(Parser)]
struct Cli {
    /// The file of SMILES strings to read as input, one SMILES per line.
    smiles_file: String,

    /// The maximum number of atoms to consider
    #[arg(long, default_value_t = 80)]
    max_atoms: usize,

    /// The force field to use for parameter labeling.
    #[arg(short, long, default_value = "openff-2.1.0.offxml")]
    forcefield: String,

    /// The parameter to use when highlighting atoms in the molecules.
    #[arg(short, long)]
    parameter: Option<String>,

    /// The `Parameter` type for which to extract parameters. Allowed options
    /// are valid arguments to `ForceField.get_parameter_handler`, such as
    /// Bonds, Angles, or ProperTorsions.
    #[arg(short, long, default_value = "ProperTorsions")]
    parameter_type: String,

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

fn main() -> io::Result<()> {
    let cli = Cli::parse();

    let s = read_to_string(&cli.smiles_file).unwrap();
    let smiles: Vec<_> = s
        .lines()
        // .map(|smiles| smiles.split('.').max_by_key(|s| s.len()).unwrap())
        .collect();
    let mols: Vec<_> = smiles
        .iter()
        .flat_map(|smiles| {
            let mut mol = ROMol::from_smiles(smiles);
            eprintln!("{}", mol.to_smiles());
            mol.openff_clean();
            if mol.num_atoms() <= cli.max_atoms {
                Some(mol)
            } else {
                None
            }
        })
        .collect();
    eprintln!("......");
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

    let existing_inchis: HashSet<_> = read_to_string("inchis.dat")
        .unwrap()
        .split_ascii_whitespace()
        .map(|s| s.to_owned())
        .collect();
    let new_inchis: Vec<_> = mols.iter().map(ROMol::to_inchi_key).collect();

    // filter out molecules with inchi_keys already covered by our existing data
    // sets; then filter any empty clusters
    clusters.iter_mut().for_each(|cluster| {
        cluster
            .retain(|mol_idx| !existing_inchis.contains(&new_inchis[*mol_idx]));
    });
    clusters.retain(|c| !c.is_empty());

    let output = Path::new(&cli.smiles_file).with_extension("html");

    Report {
        max,
        nfps,
        noise,
        clusters,
        mols,
        cli: &cli,
    }
    .generate(output)?;

    Ok(())
}
