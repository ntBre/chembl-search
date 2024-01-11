# search for structures matching ForceField parameters in the ChEMBL database

import logging

import click
from openff.toolkit import ForceField, Molecule
from openff.toolkit.utils.exceptions import RadicalsNotSupportedError
from tqdm import tqdm

logging.getLogger("openff").setLevel(logging.ERROR)


forcefield = "../projects/benchmarking/forcefields/tm-tm.offxml"
targets = ["t18b"]

ff = ForceField(forcefield, allow_cosmetic_attributes=True)


def label_molecules(self, molecule):
    from openff.toolkit import Topology

    top_mol = Topology.from_molecules([molecule])
    current_molecule_labels = dict()
    tag = "ProperTorsions"
    parameter_handler = self._parameter_handlers[tag]

    matches = parameter_handler.find_matches(top_mol)

    parameter_matches = matches.__class__()

    for match in matches:
        parameter_matches[match] = matches[match].parameter_type

    current_molecule_labels[tag] = parameter_matches

    return current_molecule_labels


@click.command()
@click.option("--output-file", "-o")
@click.option("--input-file", "-i")
def main(output_file, input_file):
    with open(input_file) as inp, open(output_file, "w") as out:
        for smiles in tqdm(inp, total=2372674):  # from wc -l
            try:
                mol = Molecule.from_smiles(
                    smiles.strip(), allow_undefined_stereo=True
                )
            except RadicalsNotSupportedError:
                continue
            labels = label_molecules(ff, mol)
            torsions = labels["ProperTorsions"]
            ids = {t.id for t in torsions.values() if t.id in targets}
            if len(ids) > 0:
                natoms = mol.n_atoms
                smiles = mol.to_smiles()
                out.write("{natoms} {ids} {smiles}\n")
                print(natoms, ids, smiles)


if __name__ == "__main__":
    main()
