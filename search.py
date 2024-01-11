# search for a SMIRKS pattern in the ChEMBL database
import gzip
import logging
import sys
import time

import numpy as np
from openff.toolkit import Molecule
from openff.toolkit.utils import RDKitToolkitWrapper
from openff.toolkit.utils.exceptions import UndefinedStereochemistryError
from openff.units import Quantity, unit
from rdkit import Chem
from tqdm import tqdm

logging.getLogger("openff").setLevel(logging.ERROR)


def from_rdkit(
    self,
    rdmol,
    allow_undefined_stereo: bool = False,
    hydrogens_are_explicit: bool = False,
    _cls=None,
):
    if _cls is None:
        _cls = Molecule

    # Make a copy of the RDKit Mol as we'll need to change it (e.g. assign
    # stereo).
    rdmol = Chem.Mol(rdmol)

    if not hydrogens_are_explicit:
        rdmol = Chem.AddHs(rdmol, addCoords=True)

    # Sanitizing the molecule. We handle aromaticity and chirality manually.
    # This SanitizeMol(...) calls cleanUp, updatePropertyCache, symmetrizeSSSR,
    # assignRadicals, setConjugation, and setHybridization.
    Chem.SanitizeMol(
        rdmol,
        (
            Chem.SANITIZE_ALL
            ^ Chem.SANITIZE_SETAROMATICITY
            ^ Chem.SANITIZE_ADJUSTHS
            ^ Chem.SANITIZE_CLEANUPCHIRALITY
        ),
    )
    Chem.SetAromaticity(rdmol, Chem.AromaticityModel.AROMATICITY_MDL)
    # SetAromaticity set aromatic bonds to 1.5, but Molecule.bond_order is an
    # integer (contrarily to fractional_bond_order) so we need the Kekule order
    Chem.Kekulize(rdmol)

    # Make sure the bond stereo tags are set before checking for
    # undefined stereo. RDKit can figure out bond stereo from other
    # information in the Mol object like bond direction properties.
    # Do not overwrite eventual chiral tags provided by the user.
    Chem.AssignStereochemistry(rdmol, cleanIt=False)

    # Check for undefined stereochemistry.
    self._detect_undefined_stereo(
        rdmol,
        raise_warning=allow_undefined_stereo,
        err_msg_prefix="Unable to make OFFMol from RDMol: ",
    )

    # Create a new OpenFF Molecule
    offmol = _cls()

    # If RDMol has a title, use it
    if rdmol.HasProp("_Name"):
        offmol.name = rdmol.GetProp("_Name")

    # Store all properties
    # TODO: Should there be an API point for storing properties?
    properties = rdmol.GetPropsAsDict()
    offmol._properties = properties

    # setting chirality in openeye requires using neighbor atoms
    # therefore we can't do it until after the atoms and bonds are all added
    map_atoms = {}
    map_bonds = {}
    # if we are loading from a mapped smiles extract the mapping
    atom_mapping = {}
    # We need the elements of the lanthanides, actinides, and transition
    # metals as we don't want to exclude radicals in these blocks.
    d_and_f_block_elements = {
        *range(21, 31),
        *range(39, 49),
        *range(57, 81),
        *range(89, 113),
    }
    for rda in rdmol.GetAtoms():
        # See issues #1075 for some discussion on radicals
        if (
            rda.GetAtomicNum() not in d_and_f_block_elements
            and rda.GetNumRadicalElectrons() != 0
        ):
            return None

        rd_idx = rda.GetIdx()
        # if the molecule was made from a mapped smiles this has been hidden
        # so that it does not affect the sterochemistry tags
        try:
            map_id = int(rda.GetProp("_map_idx"))
        except KeyError:
            map_id = rda.GetAtomMapNum()

        # create a new atom
        # atomic_number = oemol.NewAtom(rda.GetAtomicNum())
        atomic_number = rda.GetAtomicNum()
        # implicit units of elementary charge
        formal_charge = rda.GetFormalCharge()
        is_aromatic = rda.GetIsAromatic()
        if rda.HasProp("_Name"):
            name = rda.GetProp("_Name")
        else:
            # check for PDB names
            try:
                name = rda.GetMonomerInfo().GetName().strip()
            except AttributeError:
                name = ""

        # If chiral, store the chirality to be set later
        stereochemistry = None
        # tag = rda.GetChiralTag()
        if rda.HasProp("_CIPCode"):
            stereo_code = rda.GetProp("_CIPCode")
            # if tag == Chem.CHI_TETRAHEDRAL_CCW:
            if stereo_code == "R":
                stereochemistry = "R"
            # if tag == Chem.CHI_TETRAHEDRAL_CW:
            elif stereo_code == "S":
                stereochemistry = "S"
            else:
                raise UndefinedStereochemistryError(
                    "In from_rdkit: Expected atom stereochemistry of R or S. "
                    "Got {} instead.".format(stereo_code)
                )

        res = rda.GetPDBResidueInfo()
        metadata = dict()
        if res is not None:
            metadata["residue_name"] = res.GetResidueName()
            metadata["residue_number"] = res.GetResidueNumber()
            metadata["insertion_code"] = res.GetInsertionCode()
            metadata["chain_id"] = res.GetChainId()

        atom_index = offmol._add_atom(
            atomic_number,
            formal_charge,
            is_aromatic,
            name=name,
            stereochemistry=stereochemistry,
            metadata=metadata,
            invalidate_cache=False,
        )
        map_atoms[rd_idx] = atom_index
        atom_mapping[atom_index] = map_id

    offmol._invalidate_cached_properties()

    # If we have a full / partial atom map add it to the molecule. Zeroes 0
    # indicates no mapping
    if {*atom_mapping.values()} != {0}:
        offmol._properties["atom_map"] = {
            idx: map_idx
            for idx, map_idx in atom_mapping.items()
            if map_idx != 0
        }

    # Similar to chirality, stereochemistry of bonds in OE is set relative to
    # their neighbors
    for rdb in rdmol.GetBonds():
        rdb_idx = rdb.GetIdx()
        a1 = rdb.GetBeginAtomIdx()
        a2 = rdb.GetEndAtomIdx()

        # Determine bond aromaticity and Kekulized bond order
        is_aromatic = rdb.GetIsAromatic()
        order = rdb.GetBondTypeAsDouble()
        # Convert floating-point bond order to integral bond order
        order = int(order)

        # create a new bond
        bond_index = offmol._add_bond(
            map_atoms[a1],
            map_atoms[a2],
            order,
            is_aromatic=is_aromatic,
            invalidate_cache=False,
        )
        map_bonds[rdb_idx] = bond_index

    offmol._invalidate_cached_properties()

    # Now fill in the cached (structure-dependent) properties. We have to have
    # the 2D structure of the molecule in place first, because each call to
    # add_atom and add_bond invalidates all cached properties
    for rdb in rdmol.GetBonds():
        rdb_idx = rdb.GetIdx()
        offb_idx = map_bonds[rdb_idx]
        offb = offmol.bonds[offb_idx]
        # determine if stereochemistry is needed. Note that RDKit has 6
        # possible values of bond stereo: CIS, TRANS, E, Z, ANY, or NONE The
        # logic below assumes that "ANY" and "NONE" mean the same thing.
        stereochemistry = None
        tag = rdb.GetStereo()
        if tag == Chem.BondStereo.STEREOZ:
            stereochemistry = "Z"
        elif tag == Chem.BondStereo.STEREOE:
            stereochemistry = "E"
        elif (
            tag == Chem.BondStereo.STEREOTRANS
            or tag == Chem.BondStereo.STEREOCIS
        ):
            raise ValueError(
                f"Expected RDKit stereochemistry of E or Z, got {tag} instead"
            )
        offb._stereochemistry = stereochemistry
        fractional_bond_order = None
        if rdb.HasProp("fractional_bond_order"):
            fractional_bond_order = rdb.GetDoubleProp("fractional_bond_order")
        offb.fractional_bond_order = fractional_bond_order

    # TODO: Save conformer(s), if present
    # If the rdmol has a conformer, store its coordinates
    if len(rdmol.GetConformers()) != 0:
        for conf in rdmol.GetConformers():
            n_atoms = offmol.n_atoms
            # Here we assume this always be angstrom
            positions = np.zeros((n_atoms, 3))
            for rd_idx, off_idx in map_atoms.items():
                atom_coords = conf.GetPositions()[rd_idx, :]
                positions[off_idx, :] = atom_coords
            offmol._add_conformer(Quantity(positions, unit.angstrom))

    partial_charges = np.zeros(shape=offmol.n_atoms, dtype=np.float64)

    any_atom_has_partial_charge = False
    for rd_idx, rd_atom in enumerate(rdmol.GetAtoms()):
        off_idx = map_atoms[rd_idx]
        if rd_atom.HasProp("PartialCharge"):
            charge = rd_atom.GetDoubleProp("PartialCharge")
            partial_charges[off_idx] = charge
            any_atom_has_partial_charge = True
        else:
            # If some other atoms had partial charges but this one doesn't,
            # raise an Exception
            if any_atom_has_partial_charge:
                raise ValueError(
                    "Only some atoms in rdmol have partial charges"
                )
    if any_atom_has_partial_charge:
        offmol.partial_charges = Quantity(
            partial_charges, unit.elementary_charge
        )
    else:
        offmol.partial_charges = None
    return offmol


# monkey patching to avoid crash on radical error. my version returns None
# instead.
RDKitToolkitWrapper.from_rdkit = from_rdkit

target = "[*:1]~[#7:2]=[#15:3]~[*:4]"

start = time.time()
with gzip.open("chembl_33.sdf.gz") as infile:
    supplier = Chem.ForwardSDMolSupplier(
        infile, removeHs=False, sanitize=False, strictParsing=True
    )

    molecules = []
    for mol in tqdm(
        RDKitToolkitWrapper()._process_sdf_supplier(
            supplier, allow_undefined_stereo=True, _cls=None
        ),
        total=2399743,  # from chembl website, might be correct?
    ):
        if mol is not None:
            try:
                res = mol.chemical_environment_matches(target)
            except Exception as e:
                print(f"warning: {e}", file=sys.stderr)
                res = False

            if res:
                print(mol.n_atoms, mol.to_smiles())
