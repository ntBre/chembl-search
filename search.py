# search for structures matching ForceField parameters in the ChEMBL database

# using ib-dev-new environment, just random choice

import logging

import click
from openff.toolkit import ForceField, Molecule
from openff.toolkit.topology import ValenceDict
from openff.toolkit.utils.exceptions import RadicalsNotSupportedError
from openff.toolkit.utils.toolkits import RDKitToolkitWrapper
from rdkit import Chem
from tqdm import tqdm

logging.getLogger("openff").setLevel(logging.ERROR)


forcefield = "../projects/benchmarking/forcefields/tm-tm.offxml"
targets = {"t18b"}

ff = ForceField(forcefield, allow_cosmetic_attributes=True)
wrapper = RDKitToolkitWrapper()

# the big issue here is that we need to know *all* of the parameters that
# "match" a given molecule. it's not enough to know *that* a parameter matches
# because it could be superseded by a later parameter. similarly, it's not
# enough to say that something matches as soon as a single parameter matches
# because each molecule has multiple "chemical environments" that will match
# different parameters. so each molecule needs to map to a set of parameters
# that all "match" it somewhere, where match is defined by the force field to
# mean the last listed parameter that matches the environment.
#
# the comment in toolkit/parameters.py sounds easy to do, just "loop in reverse
# order, and [break] early once all environments have been matched." However,
# it's not clear (at least to me) how to get the list of all environments. if
# we had that, we could simply iterate as follows:
#
# for each environment:
#   for each parameter.rev():
#     if matches(parameter, environment):
#       save parameter
#       break 'environment
#
# it just occurred to me that this might be hard to do in general, but because
# I'm only looking at torsions, I may be able to do this exactly with
# Molecule.propers


def _find_smarts_matches(
    rdmol,
    smarts: str,
    unique: bool = False,
) -> list[tuple[int, ...]]:
    qmol = Chem.MolFromSmarts(smarts)

    assert qmol is not None

    idx_map = dict()
    for atom in qmol.GetAtoms():
        smirks_index = atom.GetAtomMapNum()
        if smirks_index != 0:
            idx_map[smirks_index - 1] = atom.GetIdx()
    map_list = [idx_map[x] for x in sorted(idx_map)]

    # UINT_MAX from limits.h
    max_matches = 4294967295
    match_kwargs = dict(
        uniquify=unique, maxMatches=max_matches, useChirality=True
    )
    full_matches = rdmol.GetSubstructMatches(qmol, **match_kwargs)

    matches = [tuple(mat[x] for x in map_list) for mat in full_matches]

    return matches


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
    # this is where I can stop with Python. if I load an RDMol straight from
    # SDF in C++, I can call the RDKit functions in _find_smarts_matches there
    return _find_smarts_matches(
        rdmol,
        smarts,
        unique=unique,
    )


def label_molecules(ff, molecule) -> set[str]:
    "returns a set of ProperTorsion ids matched"
    tag = "ProperTorsions"
    envs = ValenceDict()
    for a, b, c, d in molecule.propers:
        key = (
            a.molecule_atom_index,
            b.molecule_atom_index,
            c.molecule_atom_index,
            d.molecule_atom_index,
        )
        envs[key] = None
    remaining = len(envs)
    parameter_handler = ff._parameter_handlers[tag]
    for parameter in reversed(parameter_handler._parameters):
        env_matches = find_smarts_matches(
            wrapper,
            molecule,
            parameter.smirks,
            unique=False,
        )
        for environment_match in env_matches:
            if environment_match in envs and envs[environment_match] is None:
                envs[environment_match] = parameter.id
                remaining -= 1
            if remaining == 0:
                return {v for v in envs.values()}
    assert False


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
