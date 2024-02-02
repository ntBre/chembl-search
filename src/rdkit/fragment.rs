#![allow(unused)]

use std::{
    cell::RefCell,
    collections::{HashMap, HashSet},
    ffi::CString,
    ptr::null_mut,
    rc::Rc,
};

use super::ROMol;

const REACTION_DEFS: [&str; 12] = [
    "[#7;+0;D2,D3:1]!@C(!@=O)!@[#7;+0;D2,D3:2]>>*[#7:1].[#7:2]*", // urea
    "[C;!$(C([#7])[#7]):1](=!@[O:2])!@[#7;+0;!D1:3]>>*[C:1]=[O:2].*[#7:3]", // amide
    "[C:1](=!@[O:2])!@[O;+0:3]>>*[C:1]=[O:2].[O:3]*", // ester
    "[N;!D1;+0;!$(N-C=[#7,#8,#15,#16])](-!@[*:1])-!@[*:2]>>*[*:1].[*:2]*", // amines
    // "[N;!D1](!@[*:1])!@[*:2]>>*[*:1].[*:2]*", // amines
    // again: what about aromatics?
    "[#7;R;D3;+0:1]-!@[*:2]>>*[#7:1].[*:2]*", // cyclic amines
    "[#6:1]-!@[O;+0]-!@[#6:2]>>[#6:1]*.*[#6:2]", // ether
    "[C:1]=!@[C:2]>>[C:1]*.*[C:2]",           // olefin
    "[n;+0:1]-!@[C:2]>>[n:1]*.[C:2]*", // aromatic nitrogen - aliphatic carbon
    "[O:3]=[C:4]-@[N;+0:1]-!@[C:2]>>[O:3]=[C:4]-[N:1]*.[C:2]*", // lactam nitrogen - aliphatic carbon
    "[c:1]-!@[c:2]>>[c:1]*.*[c:2]", // aromatic carbon - aromatic carbon
    // aromatic nitrogen - aromatic carbon *NOTE* this is not part of the standard recap set.
    "[n;+0:1]-!@[c:2]>>[n:1]*.*[c:2]",
    "[#7;+0;D2,D3:1]-!@[S:2](=[O:3])=[O:4]>>[#7:1]*.*[S:2](=[O:3])=[O:4]", // sulphonamide
];

pub type Node = Rc<RefCell<RecapHierarchyNode>>;

pub fn recap_decompose(mol: &ROMol, all_nodes: Option<bool>) -> Node {
    // TODO make this lazy static
    let reactions: Vec<_> = REACTION_DEFS
        .iter()
        .map(|x| ChemicalReaction::from_smarts(x))
        .collect();

    // TODO passes 1 as arg, I think that's a bool asking for isomeric
    // smiles
    let msmi = mol.to_smiles();

    let mut all_nodes: HashMap<String, Node> = if all_nodes.is_none() {
        HashMap::new()
    } else {
        todo!("what is all_nodes supposed to be");
    };

    // TODO default arg
    let min_fragment_size = 0;
    let only_use_reactions: Option<HashSet<usize>> = None;

    if all_nodes.contains_key(&msmi) {
        return all_nodes.remove(&msmi).unwrap();
    }

    let mut res = RecapHierarchyNode::new(Rc::new(mol.clone()));
    res.smiles = Some(msmi.clone());
    let res = Rc::new(RefCell::new(res));
    let mut active_pool = HashMap::new();
    active_pool.insert(msmi.clone(), res.clone());
    all_nodes.insert(msmi, res.clone());

    while !active_pool.is_empty() {
        // this used to be while let Some((nsmi, node)) =
        // active_pool.next(), but that causes a disaster. instead, we need
        // to pop the entry out of the map like Python does
        let (nsmi, node) = {
            let nsmi = active_pool.keys().next().unwrap().clone();
            let node = active_pool.remove(&nsmi).unwrap();
            (nsmi, node)
        };
        if node.as_ref().borrow().mol.is_none() {
            continue;
        }

        for (rxn_idx, reaction) in reactions.iter().enumerate() {
            // TODO another default None arg
            if let Some(only_use_reactions) = &only_use_reactions {
                if !only_use_reactions.contains(&rxn_idx) {
                    continue;
                }
            }

            // TODO see GraphMol/ChemReactions/Reaction.h
            let ps =
                reaction.run_reactants(node.borrow().mol.as_ref().unwrap());

            // TODO maybe not is empty?
            if !ps.is_empty() {
                for prod_seq in ps {
                    let mut seq_ok = true;
                    // disqualify small fragments by sorting by size and
                    // look for "forbidden fragments"
                    let mut tseq: Vec<_> = prod_seq
                        .iter()
                        .enumerate()
                        // NOTE num_atoms might need onlyExplicit=True here
                        .map(|(idx, prod)| (prod.mol.num_atoms(), idx))
                        .collect();
                    tseq.sort();
                    // we can't move out of prod_seq itself, so turn it into
                    // a map of idx -> mol and the .remove on that
                    let mut prod_seq: HashMap<_, _> =
                        prod_seq.into_iter().enumerate().collect();
                    let ts: Vec<_> = tseq
                        .iter()
                        .map(|(x, y)| (x, prod_seq.remove(y).unwrap()))
                        .collect();
                    let mut prod_seq = ts;
                    for (nats, prod) in prod_seq.iter_mut() {
                        // default again
                        prod.mol.sanitize(super::SanitizeFlags::ALL);
                        // NOTE passing 1 again
                        let psmi = prod.mol.to_smiles();

                        // TODO default arg I think
                        if min_fragment_size > 0 {
                            let ndummies =
                                psmi.chars().filter(|&c| c == '*').count();
                            if *nats - ndummies < min_fragment_size {
                                seq_ok = false;
                                break;
                            }
                        } else if ["", "C", "CC", "CCC"].contains(
                            &psmi.replace('*', "").replace("()", "").as_str(),
                        ) {
                            // remove empty branches after replacing dummy
                            // atoms
                            seq_ok = false;
                            break;
                        }
                        prod.psmi = psmi;
                    }
                    if seq_ok {
                        for (_nats, prod) in prod_seq {
                            let psmi = prod.psmi;
                            if !all_nodes.contains_key(&psmi) {
                                let mut pnode =
                                    RecapHierarchyNode::new(Rc::new(prod.mol));
                                pnode.smiles = Some(psmi.clone());
                                // cycle, horrifying
                                pnode
                                    .parents
                                    .insert(nsmi.clone(), node.clone());
                                let pnode = Rc::new(RefCell::new(pnode));
                                node.borrow_mut()
                                    .children
                                    .insert(psmi.clone(), pnode.clone());
                                active_pool.insert(psmi.clone(), pnode.clone());
                                all_nodes.insert(psmi, pnode);
                            } else {
                                let pnode = &all_nodes[&psmi];
                                pnode
                                    .borrow_mut()
                                    .parents
                                    .insert(nsmi.clone(), node.clone());
                                node.borrow_mut()
                                    .children
                                    .insert(psmi, pnode.clone());
                            }
                        }
                    }
                }
            }
        }
    }
    res
}

struct ChemicalReaction(*mut rdkit_sys::RDKit_ChemicalReaction);

struct Product {
    mol: ROMol,
    psmi: String,
}

impl ChemicalReaction {
    fn run_reactants(&self, mol: &ROMol) -> Vec<Vec<Product>> {
        unsafe {
            let mut len = 0;
            let mut inner_lens: *mut usize = null_mut();
            let mut inner_lens_len = 0;
            let products = rdkit_sys::RDKit_RunReactants(
                self.0,
                mol.0,
                &mut len,
                &mut inner_lens,
                &mut inner_lens_len,
            );
            // Vec::from_raw_parts
            let inner_lens =
                Vec::from_raw_parts(inner_lens, inner_lens_len, inner_lens_len);
            let mut total = 0;
            for inner in &inner_lens {
                total += inner;
            }
            // returned as flat, so we need to rewrap it
            let mols = Vec::from_raw_parts(products, total, total);
            let mut idx = 0;
            let mut ret = Vec::new();
            for inner in inner_lens {
                let mut row = Vec::new();
                for _ in 0..inner {
                    let mol = ROMol(mols[idx]);
                    row.push(Product {
                        psmi: mol.to_smiles(),
                        mol,
                    });
                    idx += 1;
                }
                ret.push(row);
            }
            ret
        }
    }

    fn from_smarts(smarts: &str) -> ChemicalReaction {
        unsafe {
            let s = CString::new(smarts).expect("failed to create CString");
            let inner =
                rdkit_sys::RDKit_RxnSmartsToChemicalReaction(s.as_ptr());
            Self(inner)
        }
    }
}

pub struct RecapHierarchyNode {
    pub mol: Option<Rc<ROMol>>,
    pub children: HashMap<String, Node>,
    pub parents: HashMap<String, Node>,
    pub smiles: Option<String>,
}

impl RecapHierarchyNode {
    fn new(prod: Rc<ROMol>) -> Self {
        Self {
            mol: Some(prod),
            children: HashMap::new(),
            parents: HashMap::new(),
            smiles: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_recap() {
        let smiles = "[H]/C(=N/C([H])([H])C(=O)OC([H])([H])C([H])([H])[H])C(SSC(/C([H])=N/C([H])([H])C(=O)OC([H])([H])C([H])([H])[H])(C([H])([H])C([H])([H])[H])C([H])([H])C([H])([H])[H])(C([H])([H])C([H])([H])[H])C([H])([H])C([H])([H])[H]";
        let mol = ROMol::from_smiles(smiles);
        let node = recap_decompose(&mol, None);
        for (smi, leaf) in node.borrow().children.iter() {
            dbg!(smi);
        }
    }
}
