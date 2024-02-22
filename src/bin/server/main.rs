use std::{
    collections::HashMap,
    sync::{Arc, Mutex},
};

use axum::{routing::get, Router};
use config::{Config, Parameter};
use openff_toolkit::ForceField;
use templates::Param;
use tokio::net::TcpListener;

mod config;
mod handlers;
mod templates;

#[derive(Default)]
struct ParamState {
    smiles_list: Option<Vec<String>>,
    param_page: Option<Param>,
    nclusters: usize,
}

struct AppState {
    cli: Config,
    pid_to_smarts: HashMap<String, String>,
    param_states: HashMap<String, ParamState>,
}

impl AppState {
    #[inline]
    fn param_by_id(&self, id: &str) -> Option<&Parameter> {
        self.cli.parameters.iter().find(|p| p.id == id)
    }

    fn param_by_id_mut(&mut self, id: &str) -> Option<&mut Parameter> {
        self.cli.parameters.iter_mut().find(|p| p.id == id)
    }
}

#[tokio::main]
async fn main() {
    env_logger::init();

    let args: Vec<_> = std::env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: server <config.toml>");
        std::process::exit(1);
    }

    let cli = config::Config::load(&args[1]);

    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .unwrap();

    let pid_to_smarts: HashMap<String, String> =
        ForceField::load(&cli.forcefield)
            .unwrap()
            .get_parameter_handler("ProperTorsions")
            .unwrap()
            .parameters()
            .into_iter()
            .map(|p| (p.id(), p.smirks()))
            .collect();

    let state = Arc::new(Mutex::new(AppState {
        cli,
        param_states: HashMap::new(),
        pid_to_smarts,
    }));

    let app = Router::new()
        .route("/", get(handlers::index))
        .route("/param/:pid", get(handlers::param))
        .with_state(state);
    let listener = TcpListener::bind("0.0.0.0:3000").await.unwrap();
    axum::serve(listener, app).await.unwrap();
}