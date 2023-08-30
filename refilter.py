# filtering smiles.dat further to isolate cases where the P is not bound to 4
# other things. hopefully this is easier than inspecting each of the images
# myself

import logging

from openff.toolkit import Molecule
from rdkit.Chem.Draw import rdDepictor, rdMolDraw2D
from rdkit.Chem.rdmolops import RemoveHs
from tqdm import tqdm

logging.getLogger("openff").setLevel(logging.ERROR)

with open("smiles.dat") as inp:
    for line in tqdm(inp, desc="Searching molecules"):
        mol = Molecule.from_smiles(line.strip(), allow_undefined_stereo=True)
