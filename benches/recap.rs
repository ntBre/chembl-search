use rsearch::rdkit::{
    fragment::{recap_decompose, RecapResult},
    ROMol,
};

fn main() {
    divan::main();
}

// 2.111 ms initial mean

#[divan::bench(args = [ROMol::from_smiles("[H]/C(=N/C([H])([H])C(=O)OC([H])([H])C([H])([H])[H])\
C(SSC(/C([H])=N/C([H])([H])C(=O)OC([H])([H])C([H])([H])[H])(C([H])([H])C([H])\
([H])[H])C([H])([H])C([H])([H])[H])(C([H])([H])C([H])([H])[H])C([H])([H])C([H])\
([H])[H]")])]
fn recap(mol: &ROMol) -> RecapResult {
    recap_decompose(&mol, None)
}
