use std::{
    collections::{HashMap, HashSet},
    fs::read_to_string,
    sync::{Arc, Mutex},
};

use askama::Template;
use axum::{
    extract::{Path, Query, State},
    response::Html,
};
use log::debug;
use openff_toolkit::ForceField;
use rsearch::{
    cluster::{dbscan, Label},
    rdkit::{find_smarts_matches, fingerprint::tanimoto, ROMol},
    utils::{load_mols, make_fps, Report},
};

use crate::{
    config::Parameter,
    templates::{Body, Index, Param},
    AppState,
};

pub(crate) async fn index(
    State(state): State<Arc<Mutex<AppState>>>,
) -> Html<String> {
    Index {
        parameter_ids: state
            .lock()
            .unwrap()
            .cli
            .parameters
            .iter()
            .map(|p| p.id.clone())
            .collect(),
    }
    .render()
    .unwrap()
    .into()
}

/// load a sequence of SMILES from file, splitting on `.` to separate ionic
/// fragments, checking that `smarts` matches (NOTE that this is not a full
/// SMIRNOFF application check, just that the pattern matches at all), and
/// filtering out duplicates
fn load_smiles(
    filename: impl AsRef<std::path::Path>,
    smarts: &str,
) -> Vec<String> {
    let mut ret: Vec<String> = read_to_string(filename)
        .unwrap()
        .lines()
        .flat_map(|s| s.split('.').map(ROMol::from_smiles))
        .filter_map(|mut mol| {
            mol.openff_clean();
            if !find_smarts_matches(&mol, smarts).is_empty() {
                Some(mol.to_smiles())
            } else {
                None
            }
        })
        .collect();
    ret.sort();
    ret.dedup();
    ret
}

fn single(
    ff: &str,
    max_atoms: usize,
    radius: u32,
    param: Parameter,
) -> Result<String, std::io::Error> {
    debug!("preparing to read smiles file: {}", param.smiles);
    let s = read_to_string(&param.smiles)
        .unwrap_or_else(|e| panic!("failed to read {} for {e}", param.smiles));
    let smiles: Vec<_> = s.lines().collect();

    debug!("building map");

    // pid to smirks
    let map: HashMap<String, String> = ForceField::load(ff)
        .unwrap()
        .get_parameter_handler(&param.typ)
        .unwrap()
        .parameters()
        .into_iter()
        .map(|p| (p.id(), p.smirks()))
        .collect();

    // pid to mol
    let mol_map: Vec<(String, ROMol)> = ForceField::load(ff)
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

    debug!("loading molecules");

    let mols = load_mols(
        smiles,
        max_atoms,
        dbg!(param.fragment),
        &param.id,
        existing_inchis,
        &mol_map,
    );

    let fps: Vec<_> = make_fps(&mols, radius);
    let nfps = fps.len();
    let distance_fn = |i, j| {
        if i == j {
            0.0
        } else {
            1.0 - tanimoto(&fps[i], &fps[j])
        }
    };

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
            std::process::exit(1);
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

    let mut output = String::new();
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
    .write(&mut output)
    .unwrap();
    Ok(output)
}

pub(crate) async fn param(
    State(state): State<Arc<Mutex<AppState>>>,
    Path(pid): Path<String>,
    Query(params): Query<HashMap<String, String>>,
) -> Html<String> {
    let mut state = state.lock().unwrap();
    let smarts = state.pid_to_smarts[&pid].clone();
    let param: Parameter = state.param_by_id(&pid).unwrap().clone();
    let smiles_list = {
        let ps = state.param_states.entry(pid.clone()).or_default();
        if ps.smiles_list.is_none() {
            let collect = load_smiles(&param.smiles, &smarts);
            ps.smiles_list = Some(collect);
        }
        ps.smiles_list.clone().unwrap()
    };
    let tmpl = if params.get("eps").is_some() {
        {
            let param = state.param_by_id_mut(&pid).unwrap();
            param.dbscan.epsilon = params.get("eps").unwrap().parse().unwrap();
            param.dbscan.min_pts =
                params.get("min_pts").unwrap().parse().unwrap();
            param.fragment =
                dbg!(params.get("fragment")) == Some(&"true".to_owned());
        }
        let param = state.param_by_id_mut(&pid).unwrap().clone();
        let report = single(
            &state.cli.forcefield,
            state.cli.max_atoms,
            state.cli.radius,
            param.clone(),
        )
        .unwrap();
        Param {
            do_fragment: param.fragment,
            dbscan: param.dbscan,
            pid,
            body: Body::Report(report),
        }
    } else {
        Param {
            do_fragment: param.fragment,
            dbscan: param.dbscan,
            pid,
            body: Body::SmilesList(smiles_list),
        }
    };
    tmpl.render().unwrap().into()
}
