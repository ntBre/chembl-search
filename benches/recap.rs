use rsearch::rdkit::{
    fragment::{recap_decompose, RecapResult},
    ROMol,
};

fn main() {
    divan::main();
}

// 5.9 ms initial mean

#[divan::bench(args = [ROMol::from_smiles("            [H]C1=C(C([H])([H])[C@@]([H])(C(=O)N([H])[C@@]2([H])C(=O)\
N([H])[C@@]([H])(C([H])([H])c3c([H])nc([H])c([H])c3[H])C(=O)N([H])[C@]([H])\
(C([H])([H])C3=C([H])N([H])c4c([H])c([H])c([H])c([H])c43)C(=O)N([H])[C@@]([H])\
(C([H])([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H])C(=O)N([H])[C@@]([H])\
(C([H])([H])[C@]([H])(C([H])([H])[H])C([H])([H])S3=C([H])C([H])([H])C([H])([H])\
C3([H])[H])C(=O)N([H])[C@]([H])(C(=O)N([H])[C@]([H])(C(=O)N([H])[H])C([H])([H])\
c3c([H])c([H])c4c([H])c([H])c([H])c([H])c4c3[H])C([H])([H])SSC2([H])[H])N([H])\
[H])c2c([H])c([H])c([H])c([H])c2N1[H]")])]
fn recap(mol: &ROMol) -> RecapResult {
    recap_decompose(&mol, None)
}
