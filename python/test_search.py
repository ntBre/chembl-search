from openff.toolkit import ForceField, Molecule, Topology

from search import label_molecules

with open("../scratch/openff-toolkit/bench/smiles.100") as inp:
    moles = [
        Molecule.from_smiles(smiles.strip(), allow_undefined_stereo=True)
        for smiles in inp
    ]

forcefield = "../projects/benchmarking/forcefields/tm-tm.offxml"
ff = ForceField(forcefield, allow_cosmetic_attributes=True)


def default(self, molecule):
    top_mol = Topology.from_molecules([molecule])
    tag = "ProperTorsions"
    parameter_handler = self._parameter_handlers[tag]
    matches = parameter_handler.find_matches(top_mol)
    return {matches[m].parameter_type.id for m in matches}


got = []
want = []

for mol in moles:
    got.append(label_molecules(ff, mol))
    want.append(default(ff, mol))

assert len(got) == len(want)

for g, w in zip(got, want):
    if g != w:
        print("not ok")
        print("\t", sorted(g))
        print("\t", sorted(w))

assert got == want
