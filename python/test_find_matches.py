import logging
import sys
from collections.abc import MutableMapping

from openff.toolkit import ForceField, Molecule, Topology
from rdkit.Chem.Draw import MolsToGridImage, rdDepictor, rdMolDraw2D

logging.getLogger("openff").setLevel(logging.ERROR)


use_rdkit = False  # just to allay fears about toolkits

if use_rdkit:
    from openff.toolkit.utils import (
        GLOBAL_TOOLKIT_REGISTRY,
        OpenEyeToolkitWrapper,
    )

    GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(OpenEyeToolkitWrapper)


class _TransformedDict(MutableMapping):
    def __init__(self, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __getitem__(self, key):
        return self.store[self.__keytransform__(key)]

    def __setitem__(self, key, value):
        self.store[self.__keytransform__(key)] = value

    def __delitem__(self, key):
        del self.store[self.__keytransform__(key)]

    def __iter__(self):
        return iter(sorted(self.store, key=self.__sortfunc__))

    def __len__(self):
        return len(self.store)

    def __keytransform__(self, key):
        return key

    @staticmethod
    def __sortfunc__(key):
        return key

    @classmethod
    def _return_possible_index_of(cls, key, possible, permutations):
        raise NotImplementedError


class ValenceDict(_TransformedDict):
    @staticmethod
    def key_transform(key):
        """Reverse tuple if first element is larger than last element."""
        key = tuple(key)
        if key[0] > key[-1]:
            key = tuple(reversed(key))
        return key

    @classmethod
    def index_of(cls, key, possible=None):
        raise NotImplementedError

    def __keytransform__(self, key):
        return self.key_transform(key)


def label_molecules(self, topology):
    mols = [m for m in topology.molecules]  # stupid generators...
    assert len(mols) == 1

    top_mol = Topology.from_molecules(mols[0])
    parameter_handler = self.get_parameter_handler("ProperTorsions")
    matches = parameter_handler.find_matches(top_mol)

    parameter_matches = ValenceDict()

    for m in matches:
        # print(matches[m].parameter_type.id)
        parameter_matches[m] = matches[m].parameter_type

    return parameter_matches


ForceField.label_molecules = label_molecules


def all_matches(ff, mol):
    labels = ff.label_molecules(mol.to_topology())
    return {p.id for p in labels.values()}


def debug_matches(ff, mol):
    labels = ff.label_molecules(mol.to_topology())
    for env, p in labels.items():
        print(env, p.id)


ff = ForceField("input/tm.v2.offxml")

smiles = (
    "[H]C(=O)C([H])([H])C([H])([H])[N+]12N3[C@@]([H])(C([H])([H])C([H])([H])"
    "[C@]3([H])C1([H])[H])C2([H])[H]"
)
mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
print(smiles, end=" ")
for p in sorted(all_matches(ff, mol)):
    print(f"{p}", end=" ")
print()

if len(sys.argv) > 1:
    debug_matches(ff, mol)
