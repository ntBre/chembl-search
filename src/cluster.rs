use std::collections::HashSet;

use crate::matrix::Matrix;

#[derive(Clone, Copy, Default, PartialEq)]
pub enum Label {
    Cluster(usize),
    Noise,
    #[default]
    None,
}

impl std::fmt::Debug for Label {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let w = f.width().unwrap_or(4);
        match self {
            Label::Cluster(n) => write!(f, "{n:w$}"),
            Label::Noise => write!(f, "{:>w$}", "-1"),
            Label::None => write!(f, "{:>w$}", "N"),
        }
    }
}

fn range_query(db: &Matrix<f64>, p: usize, eps: f64) -> HashSet<usize> {
    let mut n = HashSet::new();
    let (_, c) = db.shape();
    for q in 0..c {
        if db[(p, q)] <= eps {
            n.insert(q);
        }
    }
    n
}

/// takes a distance matrix and DBSCAN parameters `eps` and `min_pts`. Using the
/// "original query-based algorithm" from Wikipedia for now
pub fn dbscan(db: &Matrix<f64>, eps: f64, min_pts: usize) -> Vec<Label> {
    let (r, c) = db.shape();
    assert_eq!(r, c, "db must be a square matrix");

    let mut label = vec![Label::None; r];
    let mut c = 0; // cluster counter

    // each row corresponds to one entry's distance from every other molecule
    for p in 0..r {
        if label[p] != Label::None {
            continue;
        }
        // here we have each of those distances, count up the neighbors of the
        // ith molecule. this is basically range_query
        let n = range_query(db, p, eps);
        if n.len() < min_pts {
            label[p] = Label::Noise;
            continue;
        }
        label[p] = Label::Cluster(c);
        let mut s = n; // S := N \ {P} since p is counted as its own neighbor
        s.remove(&p);

        let mut s: Vec<_> = s.into_iter().collect();

        let mut qi = 0;
        while qi < s.len() {
            let q = s[qi];
            if label[q] == Label::Noise {
                label[q] = Label::Cluster(c);
            }
            if label[q] != Label::None {
                qi += 1;
                continue;
            }
            label[q] = Label::Cluster(c);
            let neighbors = range_query(db, q, eps);
            if neighbors.len() >= min_pts {
                for n in neighbors {
                    if !s.contains(&n) {
                        s.push(n);
                    }
                }
            }
            qi += 1;
        }
        c += 1;
    }
    label
}

/// the input matrix and output results are taken from this page:
/// scikit-learn.org/stable/auto_examples/cluster/plot_dbscan.html, with the
/// distances precomputed using this code:
/// ``` python
/// import numpy as np
/// from sklearn import metrics
/// from sklearn.cluster import DBSCAN
/// from sklearn.datasets import make_blobs
/// from sklearn.preprocessing import StandardScaler
///
/// centers = [[1, 1], [-1, -1], [1, -1]]
/// X, labels_true = make_blobs(
///     n_samples=750, centers=centers, cluster_std=0.4, random_state=0
/// )
///
/// X = StandardScaler().fit_transform(X)
///
/// dist = np.zeros((len(X), len(X)))
/// for i in range(len(X)):
///     for j in range(len(X)):
///         dist[(i, j)] = np.linalg.norm(X[i] - X[j])
///
/// db = DBSCAN(eps=0.3, min_samples=10, metric="precomputed").fit(dist)
/// labels = db.labels_
/// ```
#[test]
fn test_dbscan() {
    let s = std::fs::read_to_string("testfiles/dist.mat").unwrap();
    let mut v: Vec<Vec<f64>> = Vec::new();
    for line in s.lines() {
        v.push(
            line.split_ascii_whitespace()
                .map(|s| s.parse().unwrap())
                .collect(),
        );
    }

    let m = Matrix::new(v);
    let got = dbscan(&m, 0.3, 10);
    let want: Vec<Label> = std::fs::read_to_string("testfiles/dist.want")
        .unwrap()
        .split_ascii_whitespace()
        .map(|s| {
            let n = s.parse::<isize>().unwrap();
            if n < 0 {
                Label::Noise
            } else {
                Label::Cluster(n as usize)
            }
        })
        .collect();
    assert_eq!(got, want);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tanimoto() {
        let mut fps = Vec::new();
        for smiles in std::fs::read_to_string("testfiles/t18b.smiles")
            .unwrap()
            .lines()
        {
            let mol = crate::rdkit::ROMol::from_smiles(smiles);
            let fp = mol.morgan_fingerprint_bit_vec::<1024>(4);
            fps.push(fp);
        }

        let n = fps.len();
        let mut got = Matrix::zeros(n, n);
        for i in 0..n {
            got[(i, i)] = 1.0;
            for j in 0..i {
                let t = crate::rdkit::fingerprint::tanimoto(&fps[i], &fps[j]);
                got[(i, j)] = t;
                got[(j, i)] = t;
            }
        }

        let want = [
            vec![
                1.0000, 0.0652, 0.0526, 0.0693, 0.9167, 0.0500, 0.0490, 0.0490,
            ],
            vec![
                0.0652, 1.0000, 0.0610, 0.8958, 0.0659, 0.0529, 0.0647, 0.0710,
            ],
            vec![
                0.0526, 0.0610, 1.0000, 0.0888, 0.0532, 0.5210, 0.5250, 0.6789,
            ],
            vec![
                0.0693, 0.8958, 0.0888, 1.0000, 0.0700, 0.0678, 0.0791, 0.0914,
            ],
            vec![
                0.9167, 0.0659, 0.0532, 0.0700, 1.0000, 0.0505, 0.0495, 0.0495,
            ],
            vec![
                0.0500, 0.0529, 0.5210, 0.0678, 0.0505, 1.0000, 0.5537, 0.4803,
            ],
            vec![
                0.0490, 0.0647, 0.5250, 0.0791, 0.0495, 0.5537, 1.0000, 0.6102,
            ],
            vec![
                0.0490, 0.0710, 0.6789, 0.0914, 0.0495, 0.4803, 0.6102, 1.0000,
            ],
        ];

        let g: Vec<_> = got.data().iter().flatten().collect();
        let w: Vec<_> = want.iter().flatten().collect();
        approx::assert_abs_diff_eq!(g.as_slice(), w.as_slice(), epsilon = 1e-4);
    }
}
