use askama::Template;

use crate::config::Dbscan;

#[derive(Template)]
#[template(path = "index.html")]
pub(crate) struct Index {
    pub(crate) parameter_ids: Vec<String>,
    pub(crate) cluster_counts: Vec<usize>,
}

#[derive(Clone)]
pub struct DrawMol {
    pub smiles: String,
    pub natoms: usize,
    pub svg: String,
}

#[derive(Clone)]
pub enum Body {
    SmilesList {
        total_mols: usize,
        mols: Vec<DrawMol>,
    },
    Report(String),
}

#[derive(Clone, Template)]
#[template(path = "param.html")]
pub(crate) struct Param {
    pub smarts: String,
    pub dbscan: Dbscan,
    pub do_fragment: bool,
    pub pid: String,
    pub body: Body,
}
