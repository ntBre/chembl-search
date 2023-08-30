# filtering smiles.dat further to isolate cases where the P is not bound to 4
# other things. hopefully this is easier than inspecting each of the images
# myself

import logging

from openff.toolkit import Molecule
from tqdm import tqdm

from visualize import draw_rdkit

logging.getLogger("openff").setLevel(logging.ERROR)

with open("smiles.dat") as inp:
    for i, line in tqdm(enumerate(inp), desc="Searching molecules"):
        mol = Molecule.from_smiles(line.strip(), allow_undefined_stereo=True)
        if mol.chemical_environment_matches("[*:1]~[#7:2]=[#15&!X4:3]~[*:4]"):
            draw_rdkit(mol, f"ref_output/mol{i:03d}.png")
            print(mol.to_smiles())
