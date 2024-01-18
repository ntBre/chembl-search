use std::{
    collections::HashSet,
    fmt::Display,
    fs::read_to_string,
    io::Write,
    ops::{Index, IndexMut},
};

use rsearch::rdkit::{fingerprint::tanimoto, ROMol};

struct Matrix<T>(Vec<Vec<T>>);

impl<T: Default + Clone> Matrix<T> {
    fn zeros(rows: usize, cols: usize) -> Self {
        Self(vec![vec![T::default(); cols]; rows])
    }
}

impl<T> Matrix<T> {
    /// (rows, cols)
    fn shape(&self) -> (usize, usize) {
        (self.0.len(), self.0.get(0).map(|v| v.len()).unwrap_or(0))
    }
}

impl<T> Index<(usize, usize)> for Matrix<T> {
    type Output = T;

    fn index(&self, (x, y): (usize, usize)) -> &Self::Output {
        &self.0[x][y]
    }
}

impl<T> IndexMut<(usize, usize)> for Matrix<T> {
    fn index_mut(&mut self, (x, y): (usize, usize)) -> &mut Self::Output {
        &mut self.0[x][y]
    }
}

impl<T: Display> Display for Matrix<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let width = f.width().unwrap_or(8);
        let prec = f.precision().unwrap_or(4);
        for row in &self.0 {
            for col in row {
                write!(f, "{col:width$.prec$}")?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
enum Label {
    Cluster(usize),
    Noise,
    #[default]
    None,
}

/// takes a distance matrix and DBSCAN parameters `eps` and `min_pts`. Using the
/// "original query-based algorithm" from Wikipedia for now
fn dbscan(db: &Matrix<f64>, eps: f64, min_pts: usize) -> Vec<Label> {
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

fn range_query(db: &Matrix<f64>, p: usize, eps: f64) -> HashSet<usize> {
    let mut n = HashSet::new();
    for q in 0..db.0[p].len() {
        if db[(p, q)] <= eps {
            n.insert(q);
        }
    }
    n
}

#[test]
fn test_dbscan() {
    let s = read_to_string("testfiles/dist.mat").unwrap();
    let mut v: Vec<Vec<f64>> = Vec::new();
    for line in s.lines() {
        v.push(
            line.split_ascii_whitespace()
                .map(|s| s.parse().unwrap())
                .collect(),
        );
    }

    let m = Matrix(v);
    let got = dbscan(&m, 0.3, 10);
    let want: Vec<Label> = read_to_string("testfiles/dist.want")
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

fn main() {
    let mut fps = Vec::new();
    // line counts for possible input smiles
    // 17562 t123a.smiles
    //    96 t138a.smiles
    //  4245 t18a.smiles
    //     8 t18b.smiles
    //    14 t87a.smiles
    for smiles in read_to_string("t138a.smiles").unwrap().lines() {
        let mol = ROMol::from_smiles(smiles);
        let fp = mol.morgan_fingerprint_bit_vec::<1024>(4);
        fps.push(fp);
    }

    let n = fps.len();
    let mut db = Matrix::zeros(n, n);
    for i in 0..n {
        db[(i, i)] = 1.0;
        for j in 0..i {
            let t = tanimoto(&fps[i], &fps[j]);
            db[(i, j)] = t;
            db[(j, i)] = t;
        }
        print!("finished mol {i}\r");
        std::io::stdout().flush().unwrap();
    }
    println!();

    // parameters taken from 2020-03-05-OpenFF-Training-Data-Selection notebook
    let labels = dbscan(&db, 0.5, 5);

    dbg!(labels);
}
