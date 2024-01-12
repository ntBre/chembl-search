# search for structures matching ForceField parameters in the ChEMBL database

# using ib-dev-new environment, just random choice

import logging

import click
from openff.toolkit import ForceField, Molecule
from openff.toolkit.topology import ValenceDict
from openff.toolkit.utils.exceptions import RadicalsNotSupportedError
from openff.toolkit.utils.toolkits import RDKitToolkitWrapper
from tqdm import tqdm

logging.getLogger("openff").setLevel(logging.ERROR)


forcefield = "../projects/benchmarking/forcefields/tm-tm.offxml"
targets = {"t18b"}

ff = ForceField(forcefield, allow_cosmetic_attributes=True)
wrapper = RDKitToolkitWrapper()


def find_smarts_matches(
    self,
    molecule: "Molecule",
    smarts: str,
    aromaticity_model: str = "OEAroModel_MDL",
    unique: bool = False,
) -> list[tuple[int, ...]]:
    rdmol = self._connection_table_to_rdkit(
        molecule, aromaticity_model=aromaticity_model
    )
    return self._find_smarts_matches(
        rdmol,
        smarts,
        aromaticity_model="OEAroModel_MDL",
        unique=unique,
    )


def _find_matches(self, molecule):
    matches = ValenceDict()
    for parameter in self._parameters:
        env_matches = find_smarts_matches(
            wrapper,
            molecule,
            parameter.smirks,
            unique=False,
        )
        for environment_match in env_matches:
            matches[environment_match] = parameter
    return matches


def label_molecules(self, molecule) -> set[str]:
    "returns a set of ProperTorsion ids matched"
    tag = "ProperTorsions"
    parameter_handler = self._parameter_handlers[tag]
    matches = _find_matches(parameter_handler, molecule)
    return {m.id for m in matches.values()}


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
            matching = labels & targets
            if len(matching) > 0:
                natoms = mol.n_atoms
                smiles = mol.to_smiles()
                out.write("{natoms} {matching} {smiles}\n")
                print(natoms, labels, smiles)


if __name__ == "__main__":
    main()
