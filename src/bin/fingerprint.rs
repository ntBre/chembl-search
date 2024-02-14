use std::{
    collections::{HashMap, HashSet},
    fs::read_to_string,
    io::{self, Write},
    path::Path,
    process::exit,
    sync::atomic::{AtomicUsize, Ordering},
};

use clap::Parser;
use log::{debug, info};
use openff_toolkit::ForceField;
use rayon::iter::{
    IntoParallelIterator, IntoParallelRefIterator, ParallelIterator,
};
use rsearch::{
    cluster::{dbscan, Label},
    rdkit::{
        bitvector::BitVector, find_smarts_matches, fingerprint::tanimoto,
        fragment::recap_decompose, ROMol,
    },
};

struct Report<'a> {
    args: Vec<String>,
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
        let mut out = std::fs::File::create(path)?;
        writeln!(out, "<html>")?;
        writeln!(out, "<pre>args: {:?}</pre>", self.args)?;
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

        if let Some(pid) = &self.cli.parameter {
            if let Some(smirks) = map.get(pid) {
                writeln!(out, "<p>PID: {pid}, SMIRKS: {smirks}</p>")?;
            }
        }

        let mut clusters = self.clusters.clone();
        clusters.sort_by_key(|c| self.mols[c[0]].num_atoms());

        for (i, c) in clusters.iter().enumerate() {
            writeln!(out, "<h1>Cluster {}, {} molecules</h1>", i + 1, c.len())?;
            self.add_svg(&mut out, "Central Molecule", &map, c[0])?;
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
        println!("{smile}");
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
                let tmp = find_smarts_matches(mol, smirks);
                if !tmp.is_empty() {
                    hl_atoms = tmp[0].clone();
                } else {
                    panic!("smirks doesn't match any more");
                }
            }
        }
        mol.draw_svg(400, 300, "", &hl_atoms)
    }
}

fn make_fps(mols: &Vec<ROMol>, radius: u32) -> Vec<BitVector> {
    mols.par_iter()
        .map(|mol| mol.morgan_fingerprint_bit_vec::<1024>(radius))
        .collect()
}

fn load_mols(
    smiles: Vec<&str>,
    max_atoms: usize,
    do_fragment: bool,
    smirks: &str,
    inchis: HashSet<String>,
) -> Vec<ROMol> {
    // this gives each of the "fragments" from the original smiles
    let mut mols: Vec<_> = smiles
        .into_iter()
        .flat_map(|s| s.split('.').map(ROMol::from_smiles))
        .collect();

    info!("collected {} initial molecules", mols.len());

    if do_fragment {
        mols = fragment(mols);
    }

    let too_big = AtomicUsize::new(0);
    let no_match = AtomicUsize::new(0);
    let overlap = AtomicUsize::new(0);
    let fragments = AtomicUsize::new(0);

    let mut ret: Vec<_> = mols
        .into_par_iter()
        .flat_map(|mut mol| {
            mol.openff_clean();
            fragments.fetch_add(1, Ordering::Relaxed);
            if mol.num_atoms() > max_atoms {
                too_big.fetch_add(1, Ordering::Relaxed);
                return None;
            }
            if find_smarts_matches(&mol, smirks).is_empty() {
                no_match.fetch_add(1, Ordering::Relaxed);
                return None;
            }
            if inchis.contains(&mol.to_inchi_key()) {
                overlap.fetch_add(1, Ordering::Relaxed);
                return None;
            }
            Some((mol.to_smiles(), mol))
        })
        .collect();

    info!(
        "expanded initial smiles to {} fragments",
        fragments.into_inner()
    );

    let presort = ret.len();

    ret.sort_by_key(|(smiles, _mol)| smiles.clone());
    ret.dedup_by_key(|(smiles, _mol)| smiles.clone());

    info!(
        "filtered {} for size, {} for smirks, {} for inchi, {} for duplicates",
        too_big.into_inner(),
        no_match.into_inner(),
        overlap.into_inner(),
        presort - ret.len(),
    );

    let (_, ret): (Vec<String>, _) = ret.into_iter().unzip();

    ret
}

fn fragment(mols: Vec<ROMol>) -> Vec<ROMol> {
    info!("starting fragment");

    // from Lily's example
    let dummy_replacements = [
        // Handle the special case of -S(=O)(=O)[*] -> -S(=O)(-[O-])
        (
            ROMol::from_smiles("S(=O)(=O)*"),
            ROMol::from_smiles("S(=O)([O-])"),
        ),
        // Handle the general case
        (ROMol::from_smiles("*"), ROMol::from_smiles("[H]")),
    ];

    /// the maximum molecule size to try fragmenting. 100 seems like a good
    /// limit, but maybe only because we're filtering a particular molecule with
    /// 103 atoms
    const MAX_FRAG_ATOMS: usize = 100;

    // apply the replacements above to a molecule and then round-trip through
    // SMILES to prevent radical formation issue
    let replace_fn = |mut m: ROMol| {
        for (inp, out) in &dummy_replacements {
            m = m.replace_substructs(inp, out, true).remove(0);
        }
        let smiles = m.to_smiles();
        ROMol::from_smiles(&smiles)
    };

    let ret = mols
        .into_par_iter()
        .flat_map(|mol| {
            let natoms = mol.num_atoms();
            if natoms > MAX_FRAG_ATOMS {
                debug!("filtered a molecule with {natoms} atoms");
                return vec![mol];
            }
            let leaves = recap_decompose(&mol, None).get_leaves();
            if leaves.is_empty() {
                return vec![mol];
            }
            leaves.into_values().map(replace_fn).collect::<Vec<_>>()
        })
        .collect();

    info!("finished fragment");

    ret
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
    #[arg(short, long, default_value = "input/tm.v2.offxml")]
    forcefield: String,

    /// The parameter to use when highlighting atoms in the molecules.
    #[arg(short, long)]
    parameter: Option<String>,

    /// The `Parameter` type for which to extract parameters. Allowed options
    /// are valid arguments to `ForceField.get_parameter_handler`, such as
    /// Bonds, Angles, or ProperTorsions.
    #[arg(long, default_value = "ProperTorsions")]
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

    /// The number of threads to use. Defaults to the number of logical CPUs as
    /// detected by rayon.
    #[arg(short, long, default_value_t = 0)]
    threads: usize,

    /// Whether or not to fragment the molecules before the fingerprinting
    /// analysis.
    #[arg(short = 'x', long, default_value_t = false)]
    fragment: bool,
}

fn main() -> io::Result<()> {
    env_logger::init();

    let cli = Cli::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .unwrap();

    let s = read_to_string(&cli.smiles_file).unwrap_or_else(|e| {
        panic!("failed to read {} for {e}", cli.smiles_file)
    });

    let parameter = cli.parameter.clone().unwrap_or_else(|| {
        Path::new(&cli.smiles_file)
            .file_stem()
            .unwrap()
            .to_str()
            .unwrap()
            .to_string()
    });

    let smiles: Vec<_> = s.lines().collect();

    info!("{} initial smiles", smiles.len());

    let map: HashMap<Pid, Smirks> = ForceField::load(&cli.forcefield)
        .unwrap()
        .get_parameter_handler(&cli.parameter_type)
        .unwrap()
        .parameters()
        .into_iter()
        .map(|p| (p.id(), p.smirks()))
        .collect();

    let existing_inchis: HashSet<_> = read_to_string("data/inchis.dat")
        .unwrap()
        .split_ascii_whitespace()
        .map(|s| s.to_owned())
        .collect();

    let mols = load_mols(
        smiles,
        cli.max_atoms,
        cli.fragment,
        &map[&parameter],
        existing_inchis,
    );

    info!("making fingerprints");

    let fps: Vec<_> = make_fps(&mols, cli.radius);

    let nfps = fps.len();
    let distance_fn = |i, j| {
        if i == j {
            0.0
        } else {
            1.0 - tanimoto(&fps[i], &fps[j])
        }
    };

    info!("running DBSCAN");

    let labels = dbscan(nfps, nfps, distance_fn, cli.epsilon, cli.min_pts);

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

    let output = Path::new(&cli.smiles_file).with_extension("html");

    info!("generating report");

    Report {
        args: std::env::args().collect::<Vec<_>>(),
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
