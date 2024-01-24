use clap::Parser;
use openff_toolkit::ForceField;
use rsearch::find_matches;
use rsearch::rdkit::ROMol;
use std::collections::{HashMap, HashSet};
use std::fs::read_to_string;
use std::io;

fn load_want(path: &str) -> HashSet<String> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|s| s.trim().to_owned())
        .collect()
}

fn print_output(res: HashMap<String, Vec<String>>) {
    for (pid, moles) in res {
        for mol in moles {
            println!("{pid}\t{mol}");
        }
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

    /// The path to the SMILES file from which to read Molecules.
    #[arg(short, long, default_value = "chembl_33.smiles")]
    molecule_file: String,

    /// The path to the file listing the parameters to match against, one per
    /// line.
    #[arg(short, long, default_value = "want.params")]
    search_params: String,

    /// Whether or not to create output files, one for each input parameter. If
    /// false, print the output to stdout.
    #[arg(short, long)]
    write_output: bool,
}

fn main() -> io::Result<()> {
    let cli = Cli::parse();
    let ff = ForceField::load(&cli.forcefield).unwrap();
    let h = ff.get_parameter_handler(&cli.parameter_type).unwrap();
    let mut params = Vec::new();
    for p in h.parameters() {
        params.push((p.id(), ROMol::from_smarts(&p.smirks())));
    }

    let want = load_want(&cli.search_params);

    let map_op = |s: &str| -> Vec<(String, String)> {
        let mut mol = ROMol::from_smiles(s);
        mol.openff_clean();
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

    let results: Vec<_> = read_to_string(cli.molecule_file)?
        .split_ascii_whitespace()
        .flat_map(map_op)
        .collect();

    let mut res: HashMap<String, Vec<String>> = HashMap::new();
    for (pid, mol) in results {
        res.entry(pid.to_string()).or_default().push(mol);
    }

    if cli.write_output {
        write_output(res);
    } else {
        print_output(res);
    }

    Ok(())
}
