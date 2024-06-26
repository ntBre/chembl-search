import logging
import pathlib
import typing
from collections import defaultdict

import click
import tqdm
from click_option_group import optgroup
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Recap

from common import canonical_smiles


def fragment_single(parent_smiles: str):
    allowed_elements = {"Br", "C", "Cl", "F", "H", "I", "N", "O", "P", "S"}
    rd_dummy_replacements = [
        # Handle the special case of -S(=O)(=O)[*] -> -S(=O)(-[O-])
        (Chem.MolFromSmiles("S(=O)(=O)*"), Chem.MolFromSmiles("S(=O)([O-])")),
        # Handle the general case
        (Chem.MolFromSmiles("*"), Chem.MolFromSmiles("[H]")),
    ]
    unique_fragments_by_n_heavy = defaultdict(set)

    rd_parent: Chem.Mol = Chem.MolFromSmiles(parent_smiles)
    if rd_parent is None:
        return unique_fragments_by_n_heavy

    if any(
        rd_atom.GetNumRadicalElectrons() != 0
        or rd_atom.GetIsotope() != 0
        or rd_atom.GetSymbol() not in allowed_elements
        for rd_atom in rd_parent.GetAtoms()
    ):
        return unique_fragments_by_n_heavy

    fragment_nodes = Recap.RecapDecompose(rd_parent).GetLeaves()
    for fragment_node in fragment_nodes.values():
        rd_fragment = fragment_node.mol

        # clean fragment
        for rd_dummy, rd_replacement in rd_dummy_replacements:
            rd_fragment = AllChem.ReplaceSubstructs(
                rd_fragment, rd_dummy, rd_replacement, True
            )[0]
            # Do a SMILES round-trip to avoid weird issues with radical
            # formation...
            rd_fragment = Chem.MolFromSmiles(Chem.MolToSmiles(rd_fragment))

        if Descriptors.NumRadicalElectrons(rd_fragment) > 0:
            logging.warning(
                f"A fragment of {parent_smiles} has a radical electron"
            )
            continue

        fragment_smiles = canonical_smiles(rd_fragment)

        if "." in fragment_smiles:
            # Skip dimers, trimers, etc.
            continue

        n_heavy = rd_fragment.GetNumHeavyAtoms()
        unique_fragments_by_n_heavy[n_heavy].add(fragment_smiles)
    return unique_fragments_by_n_heavy


def batch_fragment(smiles_list: list):
    unique_fragments_by_n_heavy = defaultdict(set)
    for parent_smiles in tqdm.tqdm(smiles_list):
        parent_fragments_by_n_heavy = fragment_single(parent_smiles)
        for n_heavy, fragment_smiles in parent_fragments_by_n_heavy.items():
            unique_fragments_by_n_heavy[n_heavy] |= fragment_smiles
    return unique_fragments_by_n_heavy


@click.option(
    "--input",
    "input_paths",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    required=True,
    multiple=True,
)
@click.option(
    "--output",
    "output_directory",
    type=click.Path(exists=False, dir_okay=True, file_okay=False),
    default="fragments",
    help="The directory to the save the generated fragments to.",
    show_default=True,
)
@optgroup.group("Parallelization configuration")
@optgroup.option(
    "--n-workers",
    help="The number of workers to distribute the labelling across. "
    "Use -1 to request one worker per batch.",
    type=int,
    default=1,
    show_default=True,
)
@optgroup.option(
    "--worker-type",
    help="The type of worker to distribute the labelling across.",
    type=click.Choice(["lsf", "local", "slurm"]),
    default="local",
    show_default=True,
)
@optgroup.option(
    "--batch-size",
    help="The number of molecules to processes at once on a single worker.",
    type=int,
    default=500,
    show_default=True,
)
@optgroup.group("LSF configuration", help="Options to configure LSF workers.")
@optgroup.option(
    "--memory",
    help="The amount of memory (GB) to request per LSF queue worker.",
    type=int,
    default=3,
    show_default=True,
)
@optgroup.option(
    "--walltime",
    help="The maximum wall-clock hours to request per LSF queue worker.",
    type=int,
    default=2,
    show_default=True,
)
@optgroup.option(
    "--queue",
    help="The LSF queue to submit workers to.",
    type=str,
    default="cpuqueue",
    show_default=True,
)
@optgroup.option(
    "--conda-environment",
    help="The conda environment that LSF workers should run using.",
    type=str,
)
@click.command()
def main(
    input_paths,
    output_directory,
    worker_type: typing.Literal["lsf", "local"] = "local",
    queue: str = "cpuqueue",
    conda_environment: str = "openff-nagl",
    memory: int = 4,  # GB
    walltime: int = 32,  # hours
    batch_size: int = 300,
    n_workers: int = -1,
):
    from dask import distributed
    from openff.nagl.utils._parallelization import batch_distributed

    # Load in the molecules to fragment
    all_parent_smiles = set()

    for input_path in input_paths:
        with open(input_path, "r") as file:
            smiles = [x.strip() for x in file.readlines()]
            smiles = [x for x in smiles if "." not in x]
            all_parent_smiles |= set(smiles)

    all_parent_smiles = sorted(all_parent_smiles, key=len, reverse=True)
    print(f"Found {len(all_parent_smiles)} parent molecules")

    uniq_frags_by_n_heavy = defaultdict(set)

    with batch_distributed(
        all_parent_smiles,
        batch_size=batch_size,
        worker_type=worker_type,
        queue=queue,
        conda_environment=conda_environment,
        memory=memory,
        walltime=walltime,
        n_workers=n_workers,
    ) as batcher:
        futures = list(batcher(batch_fragment))
        for future in tqdm.tqdm(
            distributed.as_completed(futures, raise_errors=False),
            total=len(futures),
            desc="Saving",
        ):
            fragments_by_n_heavy = future.result()
            for n_heavy, fragment_smiles in fragments_by_n_heavy.items():
                uniq_frags_by_n_heavy[n_heavy] |= fragment_smiles

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    full_path = output_directory / "fragments-full.smi"
    small_path = output_directory / "fragments-small.smi"
    # Save the fragments
    with full_path.open("w") as file:
        file.write(
            "\n".join(
                fragment_pattern
                for n_heavy, fragment_smiles in uniq_frags_by_n_heavy.items()
                for fragment_pattern in fragment_smiles
            )
        )

    small_fragments = [
        fragment_pattern
        for n_heavy in range(3, 13)
        for fragment_pattern in uniq_frags_by_n_heavy[n_heavy]
    ]

    with small_path.open("w") as file:
        file.write("\n".join(small_fragments))


if __name__ == "__main__":
    main()
