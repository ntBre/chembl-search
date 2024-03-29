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
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rdkit_rs::{find_smarts_matches, fingerprint::tanimoto, ROMol};
use rsearch::{
    cluster::{dbscan, Label},
    find_matches_full,
    utils::{load_mols, make_fps, Report},
};

use crate::{
    config::Parameter,
    templates::{Body, DrawMol, Index, Param},
    AppState,
};

pub(crate) async fn index(
    State(state): State<Arc<Mutex<AppState>>>,
) -> Html<String> {
    let state = state.lock().unwrap();
    let parameter_ids: Vec<_> =
        state.cli.parameters.iter().map(|p| p.id.clone()).collect();
    Index {
        cluster_counts: parameter_ids
            .iter()
            .map(|pid| {
                if let Some(ps) = state.param_states.get(pid) {
                    ps.nclusters
                } else {
                    0
                }
            })
            .collect(),
        parameter_ids,
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

/// returns the generated clustering report from [Report::write] as a String,
/// along with the number of clusters
fn single(
    ff: &str,
    max_atoms: usize,
    radius: u32,
    param: Parameter,
) -> Result<(String, usize), std::io::Error> {
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
    let mol_map = get_mol_map(ff, &param.typ);

    let existing_inchis: HashSet<_> = read_to_string("data/inchis.dat")
        .unwrap()
        .split_ascii_whitespace()
        .map(|s| s.to_owned())
        .collect();

    debug!("loading molecules");

    let mols = load_mols(
        smiles,
        max_atoms,
        param.fragment,
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
    Ok((output, max + 1))
}

fn get_mol_map(ff: &str, param_type: &str) -> Vec<(String, ROMol)> {
    ForceField::load(ff)
        .unwrap()
        .get_parameter_handler(param_type)
        .unwrap()
        .parameters()
        .into_iter()
        .map(|p| (p.id(), ROMol::from_smarts(&p.smirks())))
        .collect()
}

pub(crate) async fn param(
    State(state): State<Arc<Mutex<AppState>>>,
    Path(pid): Path<String>,
    Query(params): Query<HashMap<String, String>>,
) -> Html<String> {
    let mut state = state.lock().unwrap();
    let param: Parameter = state.param_by_id(&pid).unwrap().clone();
    let smarts = state.pid_to_smarts[&pid].clone();
    let smiles_list = {
        let ps = state.param_states.entry(pid.clone()).or_default();
        if ps.smiles_list.is_none() {
            let collect = load_smiles(&param.smiles, &smarts);
            ps.smiles_list = Some(collect);
        }
        ps.smiles_list.clone().unwrap()
    };
    // invalidate the cached page if max is provided. TODO save max_draw so that
    // we dont invalidate if it doesn't change
    let mut invalidate = false;
    const MAX_DRAW: usize = 50; /* the maximum number of mols to draw */
    let max_draw = if let Some(Ok(max_draw)) =
        params.get("max").map(|s| s.parse::<usize>())
    {
        invalidate = true;
        max_draw
    } else {
        MAX_DRAW
    };
    // this means Cluster button was pressed, so we overwrite whatever was there
    // before
    let tmpl = if params.contains_key("eps") {
        {
            let param = state.param_by_id_mut(&pid).unwrap();
            param.dbscan.epsilon = params.get("eps").unwrap().parse().unwrap();
            param.dbscan.min_pts =
                params.get("min_pts").unwrap().parse().unwrap();
            param.fragment = params.get("fragment") == Some(&"true".to_owned());
        }
        let param = state.param_by_id_mut(&pid).unwrap().clone();
        let (report, nclusters) = single(
            &state.cli.forcefield,
            state.cli.max_atoms,
            state.cli.radius,
            param.clone(),
        )
        .unwrap();
        state.param_states.get_mut(&pid).unwrap().nclusters = nclusters;
        Param {
            smarts,
            do_fragment: param.fragment,
            dbscan: param.dbscan,
            pid: pid.clone(),
            body: Body::Report(report),
        }
    } else if !invalidate
        && state
            .param_states
            .get(&pid)
            .is_some_and(|s| s.param_page.is_some())
    {
        state
            .param_states
            .get(&pid)
            .unwrap()
            .param_page
            .clone()
            .unwrap()
    } else {
        let mol_map = get_mol_map(&state.cli.forcefield, &param.typ);
        let mut mols: Vec<_> = smiles_list
            .into_par_iter()
            .map(|s| {
                let mut mol = ROMol::from_smiles(&s);
                // I don't really want to do this because it's expensive
                mol.openff_clean();
                let natoms = mol.num_atoms();
                (mol, s, natoms)
            })
            .collect();
        // sort first and only draw MAX_DRAW of them
        mols.sort_by_key(|(_mol, _smiles, natoms)| *natoms);
        let total_mols = mols.len();
        let mols = mols
            .into_iter()
            .take(max_draw)
            .map(|(mut mol, smiles, _natoms)| {
                mol.openff_clean();
                let (hl_atoms, _pid) = find_matches_full(&mol_map, &mol)
                    .into_iter()
                    .find(|(_atoms, param_id)| param_id == &pid)
                    .unwrap_or_default();
                let svg = mol.draw_svg(300, 300, "", &hl_atoms);
                let natoms = mol.num_atoms();
                DrawMol {
                    smiles,
                    natoms,
                    svg,
                }
            })
            .collect();
        Param {
            smarts,
            do_fragment: param.fragment,
            dbscan: param.dbscan,
            pid: pid.clone(),
            body: Body::SmilesList { mols, total_mols },
        }
    };
    let slot = state.param_states.get_mut(&pid).unwrap();
    slot.param_page = Some(tmpl.clone());
    tmpl.render().unwrap().into()
}
