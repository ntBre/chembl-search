import logging

import click
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


@click.command()
@click.option("--input-file", "-i", default="td.dat")
@click.option("--output-file", "-o", default="td.html")
def main(input_file, output_file):
    ff = ForceField("input/tm.v2.offxml")
    h = ff.get_parameter_handler("ProperTorsions")

    with open(input_file) as inp, open(output_file, "w") as out:
        last = None
        for line in inp:
            [param, smiles] = line.split()
            if param != last:
                out.write(f"<h1>{param}</h1>\n")
                smirks = h.get_parameter(dict(id=param))[0].smirks
                out.write(f"<p>{smirks}</p>\n")
            last = param
            mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
            matches = get_matches(ff, mol, param)
            svg = draw_rdkit(mol.to_rdkit(), matches[0])
            out.write(svg + "\n")


if __name__ == "__main__":
    main()
