import logging

from openff.toolkit import ForceField, Molecule
from rdkit.Chem.Draw import MolsToGridImage, rdDepictor, rdMolDraw2D

logging.getLogger("openff").setLevel(logging.ERROR)


def draw_rdkit(rdmol, matches=None):
    if matches is None:
        matches = []
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(rdmol)
    rdmol = rdMolDraw2D.PrepareMolForDrawing(rdmol)
    return MolsToGridImage(
        [rdmol],
        useSVG=True,
        highlightAtomLists=[matches],
        subImgSize=(300, 300),
        molsPerRow=1,
    )


def all_matches(ff, mol):
    labels = ff.label_molecules(mol.to_topology())[0]["ProperTorsions"]
    return {p.id for p in labels.values()}


def get_matches(ff, mol, param):
    labels = ff.label_molecules(mol.to_topology())[0]["ProperTorsions"]
    ret = []
    for env, p in labels.items():
        if p.id == param:
            ret.append(env)

    return ret


ff = ForceField("input/tm.v2.offxml")
h = ff.get_parameter_handler("ProperTorsions")

with open("td.dat") as inp:
    last = None
    for line in inp:
        [param, smiles] = line.split()
        if param != last:
            print(f"<h1>{param}</h1>")
            smirks = h.get_parameter(dict(id=param))[0].smirks
            print(f"<p>{smirks}</p>")
        last = param
        mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
        matches = get_matches(ff, mol, param)
        svg = draw_rdkit(mol.to_rdkit(), matches[0])
        print(svg)
