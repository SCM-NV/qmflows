__all__ = ['optimize_ligand']

import itertools
import numpy as np

from scm.plams.core.basemol import (Molecule, Atom, Bond)
from scm.plams.core.errors import MoleculeError
from scm.plams.core.functions import add_to_class
from scm.plams.tools.units import Units
import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit.Chem import AllChem

from .qd_functions import (to_symbol, fix_carboxyl, rotate_ligand, get_time)
from .qd_database import compare_database
from .qd_import_export import export_mol


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
            print(get_time() + 'database entry exists for ' + ligand.properties.name +
                  ' yet the corresponding .pdb file is absent. The geometry has been reoptimized.')

    return ligand


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
    mol_list = [self.copy() for i in angles]
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
def in_ring(self, arg):
    """Check if an atom or a bond belonging to this |Molecule| forms a ring.
    *arg* should be an instance of |Atom| or |Bond| belonging to this |Molecule|.
    """
    if (not isinstance(arg, (Atom, Bond))) or arg.mol != self:
        raise MoleculeError('in_ring: Argument should be a Bond or an Atom and it should be',
                            'a part of the Molecule')

    def dfs(v, depth=0):
        v._visited = True
        for bond in v.bonds:
            if bond is not arg:
                u = bond.other_end(v)
                if u is arg and depth > 1:
                    u._visited = 'cycle'
                if not u._visited:
                    dfs(u, depth+1)
    for at in self:
        at._visited = False
    if isinstance(arg, Atom):
        dfs(arg)
        ret = (arg._visited == 'cycle')
    else:
        dfs(arg.atom1)
        ret = arg.atom2._visited
    for at in self:
        del at._visited
    return ret


@add_to_class(Molecule)
def count_frags(self):
    """
    Modified PLAMS seperate() function; equavalent to len(self.separate()).
    Returns the number of fragments if the molecule were to be separated into connected components.
    """
    frags = 0
    for at in self:
        at._visited = False

    def dfs(at1):
        at1._visited = True
        for bond in at1.bonds:
            at2 = bond.other_end(at1)
            if not at2._visited:
                dfs(at2)

    for atom in self.atoms:
        if not atom._visited:
            frags += 1
            dfs(atom)
    for atom in self.atoms:
        del atom._visited

    return frags


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
    # Temporary remove hydrogen atoms
    h_atoms = []
    h_bonds = []
    for atom in reversed(plams_mol.atoms):
        if atom.atnum == 1:
            h_atoms.append(atom)
            h_bonds.append(atom.bonds[0])
            plams_mol.delete_atom(atom)

    # Remove undesired bonds
    bond_list = [bond for bond in plams_mol.bonds if not plams_mol.in_ring(bond.atom1) and not
                 plams_mol.in_ring(bond.atom2)]

    # Remove even more undesired bonds
    for bond in reversed(bond_list):
        n1, n2 = plams_mol.neighbors(bond.atom1), plams_mol.neighbors_mod(bond.atom2)
        if not (len(n1) >= 3 and len(n2) >= 2) and not (len(n1) >= 2 and len(n2) >= 3):
            bond_list.remove(bond)

    # Add the hydrogen atoms and bonds back to the molecule
    for atom, bond in zip(h_atoms, h_bonds):
        plams_mol.add_atom(atom)
        plams_mol.add_bond(bond)

    def get_atom_dict(bond_list):
        """
        Returns a dictionary with atoms as keys and all its connected bonds,
            if they are in bond_list, as values.
        """
        atom_list = list(itertools.chain.from_iterable((bond.atom1, bond.atom2) for
                                                       bond in bond_list))
        atom_set = {atom for atom in atom_list if atom_list.count(atom) >= 3}
        return {atom: [bond for bond in atom.bonds if bond in bond_list] for atom in atom_set}

    # Fragment the molecule such that the functional group is on the largest fragment
    atom_dict = get_atom_dict(bond_list)
    for at in atom_dict:
        for i in atom_dict[at][2:]:
            len_atom = [plams_mol.get_frag_size(bond, plams_mol.properties.dummies) for
                        bond in atom_dict[at]]
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


@add_to_class(Molecule)
def get_frag_size(self, bond, atom):
    """
    self <plams.Molecule>: A PLAMS molecule.
    bond <plams.Bond>: A PLAMS bond. The
    atom <plams.Atom>: A PLAMS atom. The size of the fragment containg this atom will be returned.
    """
    if bond not in self.bonds:
        error = 'get_frag_size: The argument bond should be of type plams.Bond and be part'
        error += ' of the Molecule'
        raise MoleculeError(error)
    elif atom not in self.atoms:
        error = 'get_frag_size: The argument atom should be of type plams.Atom and be part'
        error += ' of the Molecule'
        raise MoleculeError(error)

    for at in self:
        at._visited = False

    def dfs(at1, size=0, has_atom=0, atom=atom):
        at1._visited = True
        size += 1
        if at1 is atom:
            has_atom = 1
        for bond in at1.bonds:
            at2 = bond.other_end(at1)
            if not at2._visited:
                ret1, ret2 = dfs(at2)
                size += ret1
                has_atom += ret2
        return size, has_atom

    bond.atom1._visited = bond.atom2._visited = True
    size1, has_atom1 = dfs(bond.atom1)
    size2, has_atom2 = dfs(bond.atom2)
    for at in self.atoms:
        del at._visited

    if has_atom1:
        return size1
    return size2


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
                 and bond.order == 1 and not self.in_ring(bond)]

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
