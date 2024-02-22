use askama::Template;

use crate::config::Dbscan;

#[derive(Template)]
#[template(path = "index.html")]
pub(crate) struct Index {
    pub(crate) parameter_ids: Vec<String>,
}

pub enum Body {
    SmilesList(Vec<String>),
    Report(String),
}

#[derive(Template)]
#[template(path = "param.html")]
pub(crate) struct Param {
    pub dbscan: Dbscan,
    pub do_fragment: bool,
    pub pid: String,
    pub body: Body,
}
