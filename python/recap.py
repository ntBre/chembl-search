from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Recap

tot = 0
with open("output/t103.smiles") as inp:
    for line in inp:
        print(line.strip())
        print(line.strip().split("."))
        for s in line.strip().split("."):
            print(s)
            mol = Chem.MolFromSmiles(s)
            tot += len(Recap.RecapDecompose(mol).GetLeaves().values())

print(tot)
