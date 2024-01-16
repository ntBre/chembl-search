from openff.toolkit import ForceField, Molecule
from openff.toolkit.utils.toolkits import RDKitToolkitWrapper
from rdkit import Chem

smiles = "Cc1cc(-c2csc(N=C(N)N)n2)cn1C"
smarts = "[*:1]-[#16X2,#16X3+1:2]-[#6:3]~[*:4]"

mol = Molecule.from_smiles(smiles)
matches = RDKitToolkitWrapper()._find_smarts_matches(mol.to_rdkit(), smarts)
print("toolkit matches: ", matches)

rdmol = Chem.MolFromSmiles(smiles)
rdmol = Chem.AddHs(rdmol)
Chem.Kekulize(rdmol)
pat = Chem.MolFromSmarts(smarts)
print("substruct matches: ", rdmol.GetSubstructMatches(pat))

ff = ForceField("openff-2.1.0.offxml")
labels = ff.label_molecules(mol.to_topology())[0]["ProperTorsions"]
print("label matches: ")
for mat in matches:
    print(mat, labels[mat].id)

# these are from the rdkit example
# mol = Molecule.from_smiles("c1ccccc1O")
# smarts = "ccO"
# print("rdkit example: ", find_smarts_matches(wrapper, mol, smarts))

# TODO zip up fb-fit directory and include in Sage 2.2.0
