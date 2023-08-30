import logging

from openff.toolkit import Molecule
from rdkit.Chem.Draw import rdDepictor, rdMolDraw2D
from rdkit.Chem.rdmolops import RemoveHs
from tqdm import tqdm

logging.getLogger("openff").setLevel(logging.ERROR)


# copy pasta from known-issues
def draw_rdkit(mol: Molecule, filename, show_all_hydrogens=True):
    rdmol = mol.to_rdkit()
    if not show_all_hydrogens:
        rdmol = RemoveHs(rdmol, updateExplicitCount=True)
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(rdmol)
    rdmol = rdMolDraw2D.PrepareMolForDrawing(rdmol)

    drawer = rdMolDraw2D.MolDraw2DCairo(300, 300)
    drawer.DrawMolecule(rdmol)
    drawer.FinishDrawing()
    drawer.WriteDrawingText(filename)


with open("smiles.dat") as inp:
    win = 0
    failed = 0
    for i, line in tqdm(enumerate(inp), desc="Drawing molecules"):
        mol = Molecule.from_smiles(line.strip(), allow_undefined_stereo=True)
        draw_rdkit(mol, f"output/mol{i:03d}.png")

# now make images and look through them in output/*.png
