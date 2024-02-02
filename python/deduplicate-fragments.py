import typing
from collections import defaultdict

import click
import tqdm
from rdkit import Chem

from common import canonical_smiles


@click.command()
@click.option(
    "--input",
    "input_file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    required=True,
)
@click.option(
    "--reference",
    "reference_files",
    multiple=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    required=True,
)
@click.option(
    "--output",
    "output_file",
    type=click.Path(exists=False, file_okay=True, dir_okay=False),
    required=True,
)
def deduplicate(
    input_file: str,
    reference_files: typing.List[str],
    output_file: str,
):
    from openff.nagl.toolkits.openff import capture_toolkit_warnings
    from openff.toolkit.topology.molecule import Molecule

    reference_smiles = defaultdict(set)
    for reference_file in reference_files:
        with open(reference_file, "r") as f:
            all_smiles = [x.strip() for x in f.readlines()]
        for smi in tqdm.tqdm(all_smiles, desc=f"Sorting {reference_file}"):
            rdmol = Chem.MolFromSmiles(smi)
            elements = []
            for atom in rdmol.GetAtoms():
                atom.SetAtomMapNum(0)
                elements.append(atom.GetSymbol())
            elements = tuple(sorted(elements))
            reference_smiles[elements].add(canonical_smiles(rdmol))

    target_smiles = defaultdict(set)
    with open(input_file, "r") as f:
        smiles = [x.strip() for x in f.readlines()]

    print(f"Loaded {len(smiles)} target SMILES")
    for smi in tqdm.tqdm(smiles, desc="Sorting target"):
        rdmol = Chem.MolFromSmiles(smi)
        elements = []
        for atom in rdmol.GetAtoms():
            atom.SetAtomMapNum(0)
            elements.append(atom.GetSymbol())
        elements = tuple(sorted(elements))
        canonical = canonical_smiles(rdmol)
        if canonical in reference_smiles[elements]:
            continue
        else:
            target_smiles[elements].add(canonical)

    smiles_to_keep = []
    with capture_toolkit_warnings():
        for key, targets in tqdm.tqdm(
            target_smiles.items(), total=len(target_smiles)
        ):
            references = reference_smiles[key]
            reference_molecules = [
                Molecule.from_smiles(ref, allow_undefined_stereo=True)
                for ref in references
            ]
            for target in tqdm.tqdm(targets, desc=f"{key}"):
                offmol = Molecule.from_smiles(
                    target, allow_undefined_stereo=True
                )
                if any(
                    offmol.is_isomorphic_with(other)
                    for other in reference_molecules
                ):
                    continue
                else:
                    smiles_to_keep.append(target)

    print(f"Keeping {len(smiles_to_keep)} SMILES")
    with open(output_file, "w") as f:
        f.write("\n".join(smiles_to_keep))


if __name__ == "__main__":
    deduplicate()
