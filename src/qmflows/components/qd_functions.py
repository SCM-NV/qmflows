__all__ = ['optimize_ligand', 'find_substructure', 'find_substructure_split', 'rotate_ligand',
           'merge_mol', 'qd_int', 'adf_connectivity', 'fix_h', 'fix_carboxyl', 'update_coords']

import itertools
import numpy as np
import time

from scm.plams import Atom, Molecule, Bond
from scm.plams.core.functions import add_to_class
from scm.plams.core.errors import MoleculeError
from scm.plams.tools.units import Units
import scm.plams.interfaces.molecule.rdkit as molkit
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

from .qd_database import compare_database
from .qd_import_export import export_mol


time_print = '[' + time.strftime('%H:%M:%S') + '] '


@add_to_class(Molecule)
def global_minimum_scan(self, indices):
    """
    Optimize the molecule with 3 different values for the given dihedral angle and
        find the lowest energy conformer.

    :parameter self: PLAMS molecule
    :type self: plams.Molecule
    :parameter tuple indices: indices of two atoms defining a bond
    """
    # Define a number of variables and create 3 copies of the ligand
    uff = AllChem.UFFGetMoleculeForceField
    angles = (-120, 0, 120)
    mol_list = [self.copy() for i in range(angles)]
    for angle, plams_mol in zip(angles, mol_list):
        bond = plams_mol[indices]
        atom = plams_mol[indices[0]]
        plams_mol.rotate_bond(bond, atom, angle, unit='degree')

    # Optimize the geometry for all dihedral angles in angle_list
    # The geometry that yields the minimum energy is returned
    mol_list = [molkit.to_rdmol(plams_mol) for plams_mol in mol_list]
    for rdmol in mol_list:
        uff(rdmol).Minimize()
    energy_list = [uff(rdmol).CalcEnergy() for rdmol in mol_list]
    minimum = energy_list.index(min(energy_list))
    self.update_coords(mol_list[minimum], obj='rdkit')


@add_to_class(Molecule)
def to_array(self):
    """
    Convert the cartesian coordinates of atoms within a PLAMS molecule into a 3*n numpy array.
    self <plams.Molecule>: A PLAMS molecule.
    return <np.array>: A 3*n numpy array with atomic coordinates.
    """
    x, y, z = zip(*[atom.coords for atom in self])
    return np.array((x, y, z)).T


@add_to_class(Molecule)
def translate_np(self, vector, unit='angstrom'):
    """
    Translate a PLAMS molecule according to a vector using numpy.
    self <plams.Molecule>: A PLAMS molecule.
    vector: An iterable containing three floats or integers.
    """
    ratio = Units.conversion_ratio(unit, 'angstrom')
    array = self.to_array()
    array += np.array(vector)*ratio
    self.update_coords(array, obj='array')


def to_atnum(item):
    """
    Turn an atomic symbol into an atomic number.
    item <str> or <int>: An atomic symbol or number.
    return <int>: An atomic number.
    """
    if isinstance(item, str):
        return Atom(symbol=item).atnum
    return item


def to_symbol(item):
    """
    Turn an atomic number into an atomic symbol.
    item <str> or <int>: An atomic symbol or number.
    return <str>: An atomic symbol.
    """
    if isinstance(item, int):
        return Atom(atnum=item).symbol
    return item


@add_to_class(Molecule)
def update_coords(self, iterable, obj='plams'):
    """
    Update the atomic coordinates of self with coordinates from an iterable.
    self <plams.Molecule>: A PLAMS molecule.
    iterable <list>, <tuple>, <np.array>, <plams.Molecule>, etc.: A (nested) iterable containg
        x, y & z coordinates as floats.
    obj <str>: The nature of the iterable argument; accepted values:
        'plams': An iterable with PLAMS atoms as elements.
        'rdkit': An RDKit molecule.
        'array': A n*3 numpy array consisting of floats.
        'iterable': A nested iterable, each element consisting of an iterable with 3 floats.
    """
    def from_plams(self, iterable):
        for at1, at2 in zip(self, iterable):
            at1.coords = at2.coords

    def from_rdkit(self, iterable):
        conf = iterable.GetConformer()
        for at1, at2 in zip(self, iterable.GetAtoms()):
            pos = conf.GetAtomPosition(at2.GetIdx())
            at1.coords = (pos.x, pos.y, pos.z)

    def from_array(self, iterable):
        iterable = iterable.T
        for at1, x, y, z in zip(self, iterable[0], iterable[1], iterable[2]):
            at1.coords = (x, y, z)

    def from_iterable(self, iterable):
        for at1, coords in zip(self, iterable):
            at1.coords = (coords[0], coords[1], coords[2])

    obj_dict = {'plams': from_plams, 'rdkit': from_rdkit,
                'array': from_array, 'iterable': from_iterable}
    obj_dict[obj](self, iterable)


@add_to_class(Atom)
def get_atom_index(self):
    """
    Return the index of an atom (numbering starts with 1).
    self <plams.Atom>: A PLAMS atom.
    return <int>: An atomic index.
    """
    return self.mol.atoms.index(self) + 1


@add_to_class(Bond)
def get_bond_index(self):
    """
    Return a tuple of two atomic indices defining a bond (numbering starts with 1).
    self <plams.Bond>: A PLAMS bond.
    return <tuple>[<int>, <int>]: A tuple of two atomic indices defining a bond.
    """
    return self.atom1.get_atom_index(), self.atom2.get_atom_index()


@add_to_class(Atom)
def in_ring(self):
    """
    Check if this |Atom| is part of a ring. Returns a boolean.
    """
    mol = self.mol.copy()
    atom = mol[self.get_atom_index()]
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
def separate_mod(self):
    """
    Modified PLAMS function: seperates a molecule instead of a copy of a molecule.
    Separate the molecule into connected components.
    Returned is a list of new |Molecule| objects (all atoms and bonds are disjoint with
        the original molecule).
    Each element of this list is identical to one connected component of the base molecule.
    A connected component is a subset of atoms such that there exists a path
        (along one or more bonds) between any two atoms.
    """
    frags = []
    for at in self:
        at._visited = False

    def dfs(v, mol):
        v._visited = True
        v.mol = mol
        for e in v.bonds:
            e.mol = mol
            u = e.other_end(v)
            if not u._visited:
                dfs(u, mol)

    for src in self.atoms:
        if not src._visited:
            m = Molecule()
            dfs(src, m)
            frags.append(m)

    for at in self.atoms:
        del at._visited
        at.mol.atoms.append(at)
    for b in self.bonds:
        b.mol.bonds.append(b)

    return frags


@add_to_class(Molecule)
def split_bond(self, bond, element='H', length=1.1):
    """
    Delete a bond and cap the resulting fragments.
    A link to the two atoms previously defining the bond & the two capping atoms is stored under
        self.properties.mark in a list of 4-tuples.

    self <plams.Molecule>: A PLAMS molecule.
    bond <plams.Bond>: A PLAMS bond.
    element <str> or <int>: The atomic symbol or number of the two to be created capping atoms.
    length <float>: The length of the two new bonds in angstrom.
    """
    element = to_symbol(element)
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


@add_to_class(Molecule)
def neighbors_mod(self, atom, exclude=1):
    """
    A modified PLAMS function: Allows the exlucison of specific elements from the return list.
    Return a list of neighbors of *atom* within the molecule. Atoms with
    *atom* has to belong to the molecule. Returned list follows the same order as the ``bonds``
        attribute of *atom*.
    """
    if atom.mol != self:
        raise MoleculeError('neighbors: passed atom should belong to the molecule')
    return [b.other_end(atom) for b in atom.bonds if b.other_end(atom).atnum != exclude]


def split_mol(plams_mol):
    """
    Split a molecule into multiple smaller fragments for every branch within the molecule.
    plams_mol <plams.Molecule>: The input molecule with the properties.dummies attribute.
    return <list>[<plams.Molecule>] A list of one or more plams molecules.
    """
    # Remove undesired bonds
    h_mol = Molecule()
    h_mol.atoms, h_mol.bonds = zip(*[(atom, atom.bonds[0]) for atom in plams_mol.atoms if
                                     atom.atnum == 1])
    for atom, bond in zip(h_mol.atoms, h_mol.bonds):
        plams_mol.atoms.remove(atom)
        plams_mol.bonds.remove(bond)

    bond_list = [bond for bond in plams_mol.bonds if not bond.atom1.in_ring() and not
                 bond.atom1.in_ring()]

    # Remove even more undesired bonds
    for bond in reversed(bond_list):
        n1, n2 = plams_mol.neighbors_mod(bond.atom1), plams_mol.neighbors_mod(bond.atom2)
        if not (len(n1) >= 3 and len(n2) >= 2) and not (len(n1) >= 2 and len(n2) >= 3):
            bond_list.remove(bond)

    def find_dummy(mol_list):
        for mol in mol_list:
            for atom in mol:
                if plams_mol.properties.dummies.coords == atom.coords:
                    return len(mol)

    # Fragment the molecule such that the functional group is on the largest fragment
    plams_mol.merge_mol(h_mol)
    atom_list = list(itertools.chain.from_iterable((bond.atom1, bond.atom2) for bond in bond_list))
    atom_set = {atom for atom in atom_list if atom_list.count(atom) >= 3}
    atom_dict = {atom: [bond for bond in atom.bonds if bond in bond_list] for atom in atom_set}
    for at in atom_dict:
        for i in atom_dict[at][2:]:
            len_atom = []
            for bond in atom_dict[at]:
                idx = bond.get_bond_index()
                mol = plams_mol.copy()
                mol.delete_bond(mol[idx])
                mol_list = mol.separate()
                len_atom.append(find_dummy(mol_list))
            idx = len_atom.index(max(len_atom))
            bond = atom_dict[at][idx]
            plams_mol.split_bond(bond)
            atom_dict[at].remove(bond)

    # Copy the properties attribute to all fragment molecules
    properties = plams_mol.properties
    mol_list = plams_mol.separate_mod()
    for mol in mol_list:
        mol.properties = properties

    return mol_list


def recombine_mol(mol_list):
    """
    Recombine a list of molecules into a single molecule.
    A list of 4-tuples of plams.Atoms will be read from mol_list[0].properties.mark.
    A bond will be created between tuple[0] & tuple[2]; tuple[1] and tuple[3] will be deleted.
    mol_list <list>[<plams.Molecule>]: A list of on or more plams molecules with the
        properties.mark atribute.
    return <plams.Molecule>: The (re-)merged PLAMS molecule.
    """
    if len(mol_list) == 1:
        return mol_list[0]
    tup_list = mol_list[0].properties.mark
    if not tup_list:
        error = 'No PLAMS atoms specified in mol_list[0].properties.mark, aborting recombine_mol()'
        raise IndexError(error)

    for tup in tup_list:
        mol1, mol2 = tup[0].mol, tup[2].mol
        mol1.merge_mol(rotate_ligand(mol1, mol2, tup, bond_length=1.5))
        mol1.delete_atom(tup[1])
        mol1.delete_atom(tup[3])
        mol1.add_bond(tup[0], tup[2])
        bond_index = mol1.bonds[-1].get_bond_index()
        mol1.global_minimum_scan(bond_index)
    del mol1.properties.mark

    return mol1


def get_dihed(atoms, unit='degree'):
    """
    Returns the dihedral angle defined by four atoms.
    atoms <tuple>: An iterable consisting of 4 PLAMS atoms
    unit <str>: The output unit..
    return <float>: A dihedral angle.
    """
    vec1 = -1*np.array(atoms[0].vector_to(atoms[1]))
    vec2 = np.array(atoms[1].vector_to(atoms[2]))
    vec3 = np.array(atoms[2].vector_to(atoms[3]))

    v1v2, v2v3 = np.cross(vec1, vec2), np.cross(vec3, vec2)
    v1v2_v2v3 = np.cross(v1v2, v2v3)
    v2_norm_v2 = vec2/np.linalg.norm(vec2)
    epsilon = np.arctan2(np.dot(v1v2_v2v3, v2_norm_v2), np.dot(v1v2, v2v3))

    return Units.convert(epsilon, 'radian', unit)


@add_to_class(Molecule)
def set_dihed(self, angle, unit='degree'):
    """
    Change a dihedral angle into a specific value.
    self <plams.Molecule>: A PLAMS molecule.
    """
    angle = Units.convert(angle, unit, 'degree')
    bond_list = [bond for bond in self.bonds if bond.atom1.atnum != 1 and bond.atom2.atnum != 1
                 and bond.order == 1 and not bond.in_ring()]

    for bond in bond_list:
        n1, n2 = self.neighbors_mod(bond.atom1), self.neighbors_mod(bond.atom2)
        n1 = [atom for atom in n1 if atom != bond.atom2]
        n2 = [atom for atom in n2 if atom != bond.atom1]
        if len(n1) > 1:
            n1 = [atom for atom in n1 if len(self.neighbors_mod(atom)) > 1]
        if len(n2) > 1:
            n2 = [atom for atom in n2 if len(self.neighbors_mod(atom)) > 1]
        if n1 and n2:
            dihed = get_dihed((n1[0], bond.atom1, bond.atom2, n2[0]))
            self.rotate_bond(bond, bond.atom1, angle - dihed, unit='degree')

    rdmol = molkit.to_rdmol(self)
    AllChem.UFFGetMoleculeForceField(rdmol).Minimize()
    self.update_coords(rdmol, obj='rdkit')


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
        export_mol(ligand, message='Ligand:\t\t\t')

        # If ligand optimization is enabled: Optimize the ligand,
        # set pdb_info and export the result
        if opt:
            mol_list = split_mol(ligand)
            for mol in mol_list:
                mol.set_dihed(180.0)
            ligand = recombine_mol(mol_list)
            ligand = fix_carboxyl(ligand)
            ligand.properties.name = ligand.properties.name + '.opt'
            export_mol(ligand, message='Optimized ligand:\t\t')
            ligand.properties.name = ligand.properties.name.split('.opt')[0]

        # Create an entry for in the database if no previous entries are present
        # or prints a warning if a structure is present in the database but
        # the .pdb file is missing
        if not match and not pdb:
            ligand.properties.entry = True
        else:
            print(time_print + 'database entry exists for ' + ligand.properties.name +
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
        functional_group_list = ['[N+]C.[-]', '[n+]C.[-]', '[N+]c.[-]', '[n+]c.[-]',
                                 'O(C)[H]', 'O(c)[H]',
                                 'S(C)[H]', 'S(c)[H]',
                                 'N(C)[H]', 'N(c)[H]',
                                 'P(C)[H]', 'P(c)[H]',
                                 '[O-]C.[+]', '[O-]c.[+]',
                                 '[S-]C.[+]', '[S-]c.[+]',
                                 '[N-]C.[+]', '[N-]c.[+]',
                                 '[P-]C.[+]', '[P-]c.[+]']
    else:
        functional_group_list = ['[N+]C', '[n+]C', '[N+]c', '[n+]c',
                                 'OC', 'Oc', 'oC', 'oc',
                                 'SC', 'Sc', 'sC', 'sc',
                                 'NC', 'Nc', 'nC', 'nc',
                                 'PC', 'Pc', 'pC', 'pc',
                                 '[O-]C', '[O-]c',
                                 '[S-]C', '[S-]c',
                                 '[N-]C', '[N-]c', '[n-]C', '[n-]c',
                                 '[P-]C', '[P-]c', '[p-]C', '[p-]c']

    functional_group_list = [Chem.MolFromSmarts(smarts) for smarts in functional_group_list]

    # Searches for functional groups (defined by functional_group_list) within the ligand
    # Duplicates are removed
    get_match = ligand_rdkit.GetSubstructMatches
    matches = [get_match(mol) for mol in functional_group_list]
    matches = list(itertools.chain(*matches))

    # Remove all duplicate matches, each heteroatom (match[0]) should have <= 1 entry
    ligand_indices = []
    ref = []
    for match in matches:
        if match[0] not in ref:
            ligand_indices.append(match)
            ref.append(match[0])

    if ligand_indices:
        ligand_list = [ligand.copy() for match in ligand_indices]
        ligand_list = [find_substructure_split(ligand, ligand_indices[i], split) for i, ligand in
                       enumerate(ligand_list)]
    else:
        print(time_print + 'No functional groups were found for ' + str(ligand.get_formula()))
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
            mol1, mol2 = ligand.separate_mod()
            if at1 in mol1:
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
        atoms = [at.get_atom_index() for at in atoms]
    else:
        core_at1, core_at2 = core[atoms[0]], core[atoms[1]]
        lig_at1, lig_at2 = ligand[atoms[2]], ligand[atoms[3]]

    # Create two vectors defining two bonds
    core_vector = core_at1.vector_to(core_at2)
    lig_vector = lig_at2.vector_to(lig_at1)

    # Rotation of ligand - aligning the ligand and core vectors
    rotmat = create_rotmat(lig_vector, core_vector)
    xyz_array = rotmat.dot(ligand.to_array().T).T
    xyz_array += np.array(core_at1.coords)-xyz_array[atoms[2] - 1]
    if bond_length:
        vec = np.array(core_at1.vector_to(core_at2))
        xyz_array += vec*(bond_length/np.linalg.norm(vec))
    ligand.update_coords(xyz_array, obj='array')

    # Update the residue numbers
    if residue_number:
        for atom in ligand:
            atom.properties.pdb_info.ResidueNumber = residue_number + 1

    ligand.properties.anchor = lig_at1

    return ligand


def create_rotmat(vec1, vec2):
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


@add_to_class(Molecule)
def merge_mol(self, mol_list):
    """
    Merge two or more molecules into a single molecule.
    No new copies of atoms/bonds are created, all atoms/bonds are moved from mol_list to plams_mol.
    plams_mol <plams.Molecule>: A PLAMS molecule.
    mol_list <plams.Molecule> or <list>[<plams.Molecule>]: A PLAMS molecule or list of
        PLAMS molecules.
    return <plams.Molecule>: The new combined PLAMS molecule
    """
    if isinstance(mol_list, Molecule):
        mol_list = [mol_list]

    for mol in mol_list:
        self.properties.soft_update(mol.properties)

    atom_list = list(itertools.chain.from_iterable(mol.atoms for mol in mol_list))
    bond_list = list(itertools.chain.from_iterable(mol.bonds for mol in mol_list))
    for atom in atom_list:
        atom.mol = self
    for bond in bond_list:
        bond.mol = self
    self.atoms += atom_list
    self.bonds += bond_list


def adf_connectivity(plams_mol):
    """
    Create an ADF-compatible connectivity list.
    plams_mol <plams.Molecule>: A PLAMS molecule.
    return <list>[<str>]: An ADF-compatible connectivity list.
    """
    # Create list of indices of all aromatic bonds
    rdmol = molkit.to_rdmol(plams_mol)
    aromatic = [bond.GetIsAromatic() for bond in rdmol.GetBonds()]

    # Create a list of bond orders; aromatic bonds get a bond order of 1.5
    plams_mol.set_atoms_id()
    bond_orders = [bond.order for bond in plams_mol.bonds]
    for i, ar in enumerate(aromatic):
        if ar:
            bond_orders[i] = 1.5
    bonds = [str(bond.atom1.id) + ' ' + str(bond.atom2.id) + ' ' + str(order) for
             bond, order in zip(plams_mol.bonds, bond_orders)]
    plams_mol.unset_atoms_id()

    return bonds


def fix_carboxyl(plams_mol):
    """
    Resets carboxylate OCO angles if smaller than 60 degrees.
    return <plams.Molecule>: A PLAMS molecule without OCO angles smaller than 60.0 degrees.
    """
    rdmol = molkit.to_rdmol(plams_mol)
    carboxylate = Chem.MolFromSmarts('[O-]C(C)=O')
    matches = rdmol.GetSubstructMatches(carboxylate)

    if matches:
        get_angle = rdMolTransforms.GetAngleDeg
        set_angle = rdMolTransforms.SetAngleDeg
        for idx in matches:
            if get_angle(rdmol.GetConformer(), idx[3], idx[1], idx[0]) < 60:
                set_angle(rdmol.GetConformer(), idx[2], idx[1], idx[3], 180.0)
                set_angle(rdmol.GetConformer(), idx[0], idx[1], idx[3], 120.0)
        plams_mol.update_coords(rdmol, obj='rdkit')
    return plams_mol


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
        plams_mol.update_coords(rdmol, obj='rdkit')
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
