import logging
import sys

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


def debug_matches(ff, mol):
    labels = ff.label_molecules(mol.to_topology())[0]["ProperTorsions"]
    for env, p in labels.items():
        print(env, p.id)


def get_matches(ff, mol, param):
    labels = ff.label_molecules(mol.to_topology())[0]["ProperTorsions"]
    ret = []
    for env, p in labels.items():
        if p.id == param:
            ret.append(env)
        else:
            print("matched", p.id)

    return ret


ff = ForceField("input/tm.v2.offxml")
h = ff.get_parameter_handler("ProperTorsions")

smiles = (
    "[H]C(=O)C([H])([H])C([H])([H])[N+]12N3[C@@]([H])(C([H])([H])C([H])([H])"
    "[C@]3([H])C1([H])[H])C2([H])[H]"
)
smiles = "[H]N(C([H])([H])C([H])([H])[H])[N+](C([H])([H])[H])(C([H])([H])[H])C([H])([H])[H]"
mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
print(smiles, end=" ")
for p in sorted(all_matches(ff, mol)):
    print(f"{p}", end=" ")
print()

if len(sys.argv) > 1:
    debug_matches(ff, mol)

# with open("td.dat") as inp:
#     last = None
#     for line in inp:
#         [param, smiles] = line.split()
#         if param != last:
#             print(f"<h1>{param}</h1>")
#             smirks = h.get_parameter(dict(id=param))[0].smirks
#             print(f"<p>{smirks}</p>")
#         last = param
#         mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
#         matches = get_matches(ff, mol, param)
#         print(smiles)
#         svg = draw_rdkit(mol.to_rdkit(), matches[0])
#         # print(svg)