from rdkit import Chem


def canonical_smiles(rd_molecule: Chem.Mol) -> str:
    return Chem.MolToSmiles(
        Chem.AddHs(rd_molecule), isomericSmiles=False, canonical=True
    )
