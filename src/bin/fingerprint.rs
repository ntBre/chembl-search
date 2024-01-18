use std::{fs::read_to_string, io::Write};

use rsearch::{
    cluster::dbscan,
    matrix::Matrix,
    rdkit::{fingerprint::tanimoto, ROMol},
};

fn main() {
    let mut fps = Vec::new();
    // line counts for possible input smiles
    // 17562 t123a.smiles
    //    96 t138a.smiles
    //  4245 t18a.smiles
    //     8 t18b.smiles
    //    14 t87a.smiles
    // let inp = "t138a.smiles";
    let inp = "t18a.smiles";
    // let inp = "t87a.smiles";
    for smiles in read_to_string(inp).unwrap().lines() {
        let mut mol = ROMol::from_smiles(smiles);
        mol.openff_clean();
        let fp = mol.morgan_fingerprint_bit_vec::<1024>(4);
        fps.push(fp);
    }

    let start = std::time::Instant::now();

    let n = fps.len();
    let mut db = Matrix::zeros(n, n);
    // computing 1 - tanimoto here because dbscan groups items with _low_
    // distance, rather than high similarity
    for i in 0..n {
        db[(i, i)] = 0.0;
        for j in 0..i {
            let t = tanimoto(&fps[i], &fps[j]);
            db[(i, j)] = 1.0 - t;
            db[(j, i)] = 1.0 - t;
        }
        print!("finished mol {i}\r");
        std::io::stdout().flush().unwrap();
    }

    println!(
        "finished tanimoto after {:.1} s",
        start.elapsed().as_millis() as f64 / 1000.0
    );

    // default parameters from 2020-03-05-OpenFF-Training-Data-Selection
    // notebook are eps=0.5, min_pts=5
    let labels = dbscan(&db, 0.5, 5);

    for group in labels.chunks(20) {
        for g in group {
            print!("{g:?}");
        }
        println!();
    }
}
