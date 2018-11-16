__all__ = ['optimize_ligand', 'find_substructure', 'find_substructure_split', 'rotate_ligand',
           'combine_qd', 'qd_int', 'adf_connectivity', 'fix_h']

import itertools
import numpy as np

from scm.plams import Atom, Molecule
from scm.plams.core.functions import add_to_class
from scm.plams.core.errors import MoleculeError
import scm.plams.interfaces.molecule.rdkit as molkit
from rdkit import Chem
from rdkit.Chem import AllChem, Bond, rdMolTransforms

from .qd_database import compare_database
from .qd_import_export import export_mol


@add_to_class(Molecule)
def update_coords(self, plams_mol):
    """
    Update the atomic coordinates of self with coordinates from plams_mol.

    plams_mol <plams.Molecule>: A PLAMS molecule.
    """
    for at, at_self in zip(plams_mol, self):
        at_self.coords = at.coords


@add_to_class(Atom)
def get_index(self):
    """
    Return the index of an |Atom| (numbering starts with 1).
    """
    return self.mol.atoms.index(self) + 1


@add_to_class(Bond)
def get_index(self):
    """
    Return a tuple of two atomic indices defining a |Bond| (numbering starts with 1).
    """
    return self.atom1.get_index(), self.atom2.get_index()


@add_to_class(Atom)
def in_ring(self):
    """
    Check if this |Atom| is part of a ring. Returns a boolean.
    """
    mol = self.mol.copy()
    atom = mol[self.get_index()]
    before = len(mol.separate())
    neighbors = len(atom.bonds)
    bonds = [bond for bond in atom.bonds]
    for bond in bonds:
        mol.delete_bond(bond)
    after = len(mol.separate())
    if before != after - neighbors:
        return True
    return False


@add_to_class(Molecule)
def neighbors_mod(self, atom, exclude_atnum=[]):
    """Return a list of neighbors of *atom* within the molecule.

    *atom* has to belong to the molecule. Returned list follows the same order as
        the ``bonds`` attribute of *atom*.

    *exclude_atnum* is an optional list of atomic numbers (int).
    If the atomic number/symbol of an |Atom| matches with any entry in *exclude*,
        it will not be returned.
    """
    if atom.mol != self:
        raise MoleculeError('neighbors: passed atom should belong to the molecule')
    return [b.other_end(atom) for b in atom.bonds if b.other_end(atom).atnum not in exclude_atnum]


@add_to_class(Molecule)
def split_bond(self, bond, element='H', length=1.1):
    """
    Delete a bond and cap the resulting fragments.
    A link to the two atoms previously defining the bond & the two capping atoms is stored under
        self.properties.mark in a list of 4-tuples.

    self <plams.Molecule>: A PLAMS molecule.
    bond <plams.Bond>: A PLAMS bond.
    element <str> or <int>: The atomic symbol or number of the two to be created capping atoms.
    resize <float>: The length of the two new bonds in angstrom.
    """
    if isinstance(element, int):
        element = Atom(atnum=element).symbol
    at1, at2 = bond.atom1, bond.atom2
    at3, at4 = Atom(symbol=element, coords=at1.coords), Atom(symbol=element, coords=at2.coords)
    self.add_atom(at3, adjacent=[at2])
    self.add_atom(at4, adjacent=[at1])
    self.bonds[-1].resize(at1, length)
    self.bonds[-2].resize(at2, length)
    if self.properties.mark:
        self.properties.mark.append((at1, at4, at2, at3))
    else:
        self.properties.mark = [(at1, at4, at2, at3)]
    self.delete_bond(bond)


def split_mol(plams_mol, exclude_atnum=[1]):
    """
    Split a molecule into multiple smaller fragments for every branch within the molecule.


    plams_mol <plams.Molecule>: The input molecule.
    exclude_atom <list>[<plams.Atom>]: A list of one or more atoms that will not be returend by
        the neighbors_mod() function.
    exclude_element <list>[<int> or <str>]: A list of atomic numbers or symbols of atoms that will
        not be returned by the neighbors_mod() function.

    return <list>[<plams.Molecule>] A list of one or more plams molecules.
    """
    bond_list = [bond for bond in plams_mol.bonds if
                 bond.atom1.atnum != 1 and bond.atom2.atnum != 1]
    for bond in bond_list:
        neighbors1 = plams_mol.neighbors_mod(bond.atom1, exclude_atnum=exclude_atnum)
        neighbors2 = plams_mol.neighbors_mod(bond.atom2, exclude_atnum=exclude_atnum)
        if len(neighbors1) > 2 and len(neighbors2) > 2 and not bond.in_ring():
            plams_mol.split_bond(bond)
    properties = plams_mol.properties
    mol_list = plams_mol.separate()
    for mol in mol_list:
        mol.properties = properties

    return mol_list


def recombine_mol(mol_list, mol_old=False):
    """
    Recombine a list of molecules into a single molecule.
    A list of 4-tuples of plams.Atoms will be read from mol_list[0].
    A bond will be created between tuple[0] & tuple[2] and atoms tuple[1] and tuple[3]
        will be deleted.

    mol_list <list>[<plams.Molecule>, ...]: A list of n plams molecules with the
        atom.properties.mark atribute.

    return <plams.Molecule>: The (re-)merged PLAMS molecule.
    """
    # Create a set of all mol.mark values in mol_list
    tup_list = mol_list[0].properties.mark
    for tup in tup_list:
        mol1, mol2 = tup[0].mol, tup[2].mol
        mol1.merge_mol(rotate_ligand(mol1, mol2, tup, bond_length=1.5))
        mol1.delete_atom(tup[1])
        mol1.delete_atom(tup[2])
        mol1.add_bond(tup[0], tup[2])
        tup_list.remove(tup)

    return mol_list[0]


def optimize_ligand(ligand, database, opt=True):
    """
    Pull the structure if a match has been found or alternatively optimize a new geometry.

    ligand <plams.Molecule>: The ligand molecule.
    database <pd.DataFrame>: Database of previous calculations.
    opt <bool>: If the geometry of the ligand (RDKit UFF) should be optimized (True) or not (False).

    return <plams.Molecule>: The optimized ligand molecule.
    """
    # Searches for matches between the input ligand and the database; imports the structure
    if database is not None:
        ligand, match, pdb = compare_database(ligand, database)
    else:
        match, pdb = False, False

    # Optimize the ligand if no match has been found with the database
    ligand.properties.entry = False
    if not match or not pdb:
        # Export the unoptimized ligand to a .pdb and .xyz file
        export_mol(ligand, message='Ligand:\t\t\t\t')

        # If ligand optimization is enabled: Optimize the ligand,
        # set pdb_info and export the result
        if opt:
            ligand_opt = molkit.global_minimum(ligand, n_scans=2, no_h=True)
            ligand_opt.properties.name = ligand.properties.name + '.opt'
            ligand_opt.properties.path = ligand.properties.path
            export_mol(ligand_opt, message='Optimized ligand:\t\t')
            for i, atom in enumerate(ligand):
                atom.move_to(ligand_opt[i + 1])

        # Create an entry for in the database if no previous entries are present
        # or prints a warning if a structure is present in the database but
        # the .pdb file is missing
        if not match and not pdb:
            ligand.properties.entry = True
        else:
            print('\ndatabase entry exists for ' + ligand_opt.properties.name +
                  ' yet the corresponding .pdb file is absent. The geometry has been reoptimized.')

    return ligand


def find_substructure(ligand, split=True):
    """
    Identify the ligand functional groups.

    ligand <plams.Molecule>: The ligand molecule.
    split <bool>: If a functional group should be split from the ligand (True) or not (False).

    return <list>[<plams.Molecule>]: A copy of the ligand for each identified functional group.
    """
    ligand_rdkit = molkit.to_rdmol(ligand)

    # Creates a list containing predefined functional groups, each saved as an rdkit molecule
    # IMPORTANT: The first atom should ALWAYS be the atom that should attach to the core
    if split:
        functional_group_list = ['[N+]C.[-]',
                                 'O(C)[H]',
                                 'S(C)[H]',
                                 'N(C)[H]',
                                 'P(C)[H]',
                                 '[O-]C.[+]',
                                 '[S-]C.[+]',
                                 '[N-]C.[+]',
                                 '[P-]C.[+]']
    else:
        functional_group_list = ['[N+]C',
                                 'O[H]C',
                                 'S[H]C',
                                 'N[H]C',
                                 'P[H]C',
                                 '[O-]C',
                                 '[S-]C',
                                 '[N-]C',
                                 '[P-]C']

    functional_group_list = [Chem.MolFromSmarts(smarts) for smarts in functional_group_list]

    # Searches for functional groups (defined by functional_group_list) within the ligand
    # Duplicates are removed
    get_match = ligand_rdkit.GetSubstructMatches
    matches = [get_match(mol) for mol in functional_group_list]
    matches = list(itertools.chain(*matches))

    # Remove all duplicate matches, each heteroatom (match[0]) should have <= 1 entry
    ligand_indices = []
    for match in matches:
        if match[0] not in [item[0] for item in ligand_indices]:
            ligand_indices.append(match)

    if ligand_indices:
        ligand_list = [ligand.copy() for match in ligand_indices]
        ligand_list = [find_substructure_split(ligand, ligand_indices[i], split) for i, ligand in
                       enumerate(ligand_list)]
    else:
        print('No functional groups were found for ' + str(ligand.get_formula()))
        ligand_list = []

    return ligand_list


def find_substructure_split(ligand, ligand_index, split=True):
    """
    Delete the hydrogen or mono-/polyatomic counterion attached to the functional group.
    Sets the charge of the remaining heteroatom to -1 if split=True.

    ligand <plams.Molecule>: The ligand molecule.
    ligand_index <list>[<int>, <int>]: A list of atomic indices associated with a functional group.
    split <bool>: If a functional group should be split from the ligand (True) or not (False).

    return <plams.Molecule>: The ligand molecule.
    """
    at1 = ligand[ligand_index[0] + 1]
    at2 = ligand[ligand_index[-1] + 1]
    ligand.properties.group = at1.symbol + str(ligand_index[0] + 1)

    if split:
        if len(ligand.separate()) == 1:
            ligand.delete_atom(at2)
        else:
            mol1, mol2 = ligand.separate()
            if str(at1) in [str(atom) for atom in mol1]:
                ligand = mol1
            else:
                ligand = mol2

        # Check if the ligand heteroatom has a charge assigned, assigns a charge if not
        if not at1.properties.charge or at1.properties.charge == 0:
            at1.properties.charge = -1

    # Update the index of the ligand heteroatom
    ligand.properties.dummies = at1

    # Set the molecular charge
    ligand.properties.charge = sum([atom.properties.charge for atom in ligand
                                    if atom.properties.charge])

    return ligand


def rotate_ligand(core, ligand, atoms, bond_length=False, residue_number=False):
    """
    Connects two molecules by alligning the vectors of two bonds.

    core <plams.Molecule>: The core molecule.
    ligand <plams.Molecule>: The ligand molecule.
    i <int>: The residue number that is to be assigned to the ligand.
    core_dummy <plams.Atom>: A dummy atom at the core center of mass.

    return <plams.Molecule>, <plams.Atom>: The rotated ligand and the ligand atom that will be
        attached to the core.
    """

    # Defines first atom on coordinate list (hydrogen),
    # The atom connected to it and vector representing bond between them
    if isinstance(atoms[0], Atom):
        core_at1, core_at2 = atoms[0], atoms[1]
        lig_at1, lig_at2 = atoms[2], atoms[3]
    else:
        core_at1, core_at2 = core[atoms[0]], core[atoms[1]]
        lig_at1, lig_at2 = ligand[atoms[2]], ligand[atoms[3]]

    # Create two vectors defining two bonds
    core_vector = core_at1.vector_to(core_at2)
    lig_vector = lig_at2.vector_to(lig_at1)

    # Rotation of ligand - aligning the ligand and core vectors
    rotmat = rotate_ligand_rotation(lig_vector, core_vector)
    ligand.rotate(rotmat)
    ligand.translate(lig_at1.vector_to(core_at1))

    # Translation of the ligand
    if bond_length:
        vec = np.array(core_at1.vector_to(core_at2))
        vec = vec*(bond_length/np.linalg.norm(vec))
        ligand.translate(vec)

    # Update the residue numbers
    if residue_number:
        for atom in ligand:
            atom.properties.pdb_info.ResidueNumber = residue_number + 1

    ligand.properties.anchor = lig_at1

    return ligand


def rotate_ligand_rotation(vec1, vec2):
    """
    Calculates the rotation matrix rotating *vec1* to *vec2*.
    Vectors can be any containers with 3 numerical values. They don't need to be normalized.
    Returns 3x3 numpy array.
    """
    a = np.array(vec1) / np.linalg.norm(vec1)
    b = np.array(vec2) / np.linalg.norm(vec2)
    v1, v2, v3 = np.cross(a, b)
    M = np.array([[0, -v3, v2], [v3, 0, -v1], [-v2, v1, 0]])

    return np.identity(3) + M + np.dot(M, M)/(1+np.dot(a, b))


def combine_qd(core, ligand_list):
    """
    Combine the rotated ligands with the core, creating a bond bewteen the core and ligand.

    core <plams.Molecule>: The core molecule.
    ligand_list <list>[<plams.Molecule>]: A list of the rotated ligands.

    return <plams.Molecule>: The quantum dot (core + n*ligands).
    """
    qd = core.copy()

    # Create a list of ligand atoms and intraligand bonds
    ligand_bonds = np.concatenate([ligand.bonds for ligand in ligand_list])
    ligand_atoms = np.concatenate(ligand_list)

    # Combined the ligand bond and atom list with the core
    for atom in ligand_atoms:
        qd.add_atom(atom)
    for bond in ligand_bonds:
        qd.add_bond(bond)

    # Delete redundant atoms
    indices = [atom.get_index() for atom in core.properties.dummies]
    indices += [atom.get_index() for atom in qd if atom.atnum == 0]
    indices.sort(reverse=True)
    for index in indices:
        qd.delete_atom(qd[index])

    return qd


def adf_connectivity(plams_mol):
    """
    Create an ADF-compatible connectivity list.

    plams_mol <plams.Molecule>: A PLAMS molecule.

    return <list>[<str>]: An ADF-compatible connectivity list.
    """
    # Create list of indices of all aromatic bonds
    rdmol = molkit.to_rdmol(plams_mol)
    aromatic = [Bond.GetIsAromatic(bond) for bond in rdmol.GetBonds()]
    aromatic = [i for i, bond in enumerate(aromatic) if bond]

    # Create a connectivity list; aromatic bonds get a bond order of 1.5
    at1 = [plams_mol.atoms.index(bond.atom1) + 1 for bond in plams_mol.bonds]
    at2 = [plams_mol.atoms.index(bond.atom2) + 1 for bond in plams_mol.bonds]
    bonds = [bond.order for bond in plams_mol.bonds]
    for i, bond in enumerate(plams_mol.bonds):
        if i in aromatic:
            bonds[i] = 1.5
    bonds = [str(at1[i]) + ' ' + str(at2[i]) + ' ' + str(bond) for i, bond in enumerate(bonds)]

    return bonds


def fix_h(plams_mol):
    """
    If a C=C-H angle is smaller than 20.0 degrees, set it back to 120.0 degrees.

    plams_mol <plams.Molecule>: A PLAMS molecule.

    return <plams.Molecule>: A PLAMS molecule without C=C-H angles smaller than 20.0 degrees.
    """
    H_list = [atom for atom in plams_mol if atom.atnum is 1 and 2.0 in
              [bond.order for bond in plams_mol.neighbors(atom)[0].bonds]]

    rdmol = molkit.to_rdmol(plams_mol)
    idx = plams_mol.atoms.index
    set_angle = rdMolTransforms.SetAngleDeg
    get_angle = rdMolTransforms.GetAngleDeg

    update = []
    for atom in H_list:
        at1 = atom
        at2 = plams_mol.neighbors(at1)[0]
        at3 = [atom for atom in plams_mol.neighbors(at2) if atom != at1]
        if get_angle(rdmol.GetConformer(), idx(at3[0]), idx(at2), idx(at1)) <= 20.0:
            set_angle(rdmol.GetConformer(), idx(at3[0]), idx(at2), idx(at1), 120.0)
            update.append(True)
        elif get_angle(rdmol.GetConformer(), idx(at3[1]), idx(at2), idx(at1)) <= 20.0:
            set_angle(rdmol.GetConformer(), idx(at3[1]), idx(at2), idx(at1), 120.0)
            update.append(True)

    if update:
        return molkit.from_rdmol(rdmol)
    else:
        return plams_mol


def qd_int(plams_mol):
    """
    Perform an activation-strain analyses (RDKit UFF) on the ligands in the absence of the core.

    plams_mol <plams.Molecule>: A PLAMS molecule.
    job <str>: The to be executed AMS job (see qd_ams.py).

    return <plams.Molecule>: A PLAMS molecule with the int and int_mean properties.
    """
    mol_copy = plams_mol.copy()
    uff = AllChem.UFFGetMoleculeForceField

    # Calculate the total energy of all perturbed ligands in the absence of the core
    atom_list = [atom for atom in mol_copy if atom.properties.pdb_info.ResidueName == 'COR']
    for atom in atom_list:
        mol_copy.delete_atom(atom)
    rdmol = molkit.to_rdmol(mol_copy)
    E_no_frag = uff(rdmol, ignoreInterfragInteractions=False).CalcEnergy()

    # Calculate the total energy of the isolated perturbed ligands in the absence of the core
    mol_frag = mol_copy.separate()
    E_frag = 0.0
    for mol in mol_frag:
        rdmol = molkit.to_rdmol(mol)
        E_frag += uff(rdmol, ignoreInterfragInteractions=False).CalcEnergy()

    # Calculate the total energy of the optimized ligand
    uff(rdmol, ignoreInterfragInteractions=False).Minimize()
    E_opt = uff(rdmol, ignoreInterfragInteractions=False).CalcEnergy()

    # Calculate E, Eint and Estrain
    plams_mol.properties.Eint = float(E_no_frag - E_frag)
    plams_mol.properties.Estrain = float(E_frag - (E_opt * len(mol_frag)))
    plams_mol.properties.E = plams_mol.properties.Eint + plams_mol.properties.Estrain

    return plams_mol
