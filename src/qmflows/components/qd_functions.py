__all__ = ['optimize_ligand', 'find_substructure', 'find_substructure_split', 'rotate_ligand',
           'combine_qd', 'qd_opt', 'qd_int']

import itertools
import os
import numpy as np

from scm.plams import (Atom)
import scm.plams.interfaces.molecule.rdkit as molkit
from rdkit import Chem
from rdkit.Chem import AllChem

from .qd_database import compare_database
from .qd_import_export import export_mol
from .qd_ams import ams_job


def optimize_ligand(ligand, opt, database):
    """
    Pull the structure if a match has been found or alternatively optimize a new geometry.
    """
    # Searches for matches between the input ligand and the database; imports the structure
    if not isinstance(database, bool):
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
            ligand_opt.properties.source_folder = ligand.properties.source_folder
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

    ligand.properties.source = os.path.join(ligand.properties.source_folder,
                                            ligand.properties.name + '.opt.pdb')
    return ligand


def find_substructure(ligand, split=True):
    """
    Identify the ligand functional groups.
    """
    ligand_rdkit = molkit.to_rdmol(ligand)

    # Creates a list containing predefined functional groups, each saved as an rdkit molecule
    # IMPORTANT: The first atom should ALWAYS be the atom that should attach to the core
    if split:
        functional_group_list = ['[N+].[-]',
                                 'O[H]',
                                 'S[H]',
                                 'N[H]',
                                 'P[H]',
                                 '[O-].[+]',
                                 '[S-].[+]',
                                 '[N-].[+]',
                                 '[P-].[+]']
    else:
        functional_group_list = ['[N+]',
                                 'O[H]',
                                 'S[H]',
                                 'N[H]',
                                 'P[H]',
                                 '[O-]',
                                 '[S-]',
                                 '[N-]',
                                 '[P-]']

    functional_group_list = [Chem.MolFromSmarts(smarts) for smarts in functional_group_list]

    # Searches for functional groups (defined by functional_group_list) within the ligand
    # Duplicates are removed
    get_match = ligand_rdkit.GetSubstructMatches
    matches = [get_match(mol) for mol in functional_group_list]
    matches = list(itertools.chain(*matches))

    # Remove all duplicate matches
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
    Delete the hydrogen or mono-/polyatomic counterion attached to the functional group
    Sets the charge of the remaining heteroatom to -1 if split=True
    """
    at1 = ligand[ligand_index[0] + 1]
    at2 = ligand[ligand_index[1] + 1]
    ligand.properties.name += '@' + at1.symbol + str(ligand_index[0] + 1)

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
    ligand.properties.ligand_dummy = at1
    ligand.add_atom(Atom(atnum=0, coords=ligand.get_center_of_mass()))

    return ligand


def rotate_ligand(core, ligand, i, core_dummy):
    """
    Connects two molecules by alligning the vectors of two bonds.
    """
    ligand = ligand.copy()
    ligand.properties.ligand_dummy = ligand.closest_atom(ligand.properties.ligand_dummy.coords)

    # Defines first atom on coordinate list (hydrogen),
    # The atom connected to it and vector representing bond between them
    core_at1 = core_dummy         # core dummy atom
    core_at2 = core[-1]                 # core center of mass
    core_vector = core_at1.vector_to(core_at2)
    lig_at1 = ligand.properties.ligand_dummy  	# ligand heteroatom
    lig_at2 = ligand[-1]                # ligand center of mass
    lig_vector = lig_at2.vector_to(lig_at1)

    # Rotation of ligand - aligning the ligand and core vectors
    rotmat = rotate_ligand_rotation(lig_vector, core_vector)
    ligand.rotate(rotmat)
    ligand.translate(lig_at1.vector_to(core_at1))

    # Translation of the ligand
    hc_vec = lig_at1.vector_to(core_at1)
    ligand.translate(hc_vec)

    # Update the residue numbers
    for atom in ligand:
        atom.properties.pdb_info.ResidueNumber = i + 2

    # Deletes the core dummy atom and ligand center of mass
    ligand.delete_atom(lig_at2)
    core.delete_atom(core_at1)

    return ligand, lig_at1


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

    return qd


def qd_int(plams_mol, job='qd_sp'):
    """
    Perform an activation-strain analyses (RDKit UFF) on the ligands in the absence of the core.

    plams_mol <plams.Molecule>: A PLAMS molecule
    job <str>: The to be executed AMS job (see qd_ams.py)
    return <plams.Molecule>: A PLAMS molecule with the int and int_mean properties.
    """
    mol_copy = plams_mol.copy()
    uff = AllChem.UFFGetMoleculeForceField

    # Calculate the total energy of all ligands in the absence of the core
    atom_list = [atom for atom in mol_copy if atom.properties.pdb_info.ResidueName is 'COR']
    for atom in atom_list:
        mol_copy.delete_atom(atom)
    rdmol = molkit.to_rdmol(mol_copy)
    E_no_frag = uff(rdmol, ignoreInterfragInteractions=False).CalcEnergy()

    # Calculate the total energy of the isolated ligands in the absence of the core
    mol_frag = mol_copy.separate()

    E_frag = 0.0
    for mol in mol_frag:
        rdmol = molkit.to_rdmol(mol)
        E = uff(rdmol, ignoreInterfragInteractions=False).CalcEnergy()
        print(E)
        E_frag += E

    # Calculate the total energy of an optimized ligand
    rdmol = molkit.to_rdmol(mol_frag[0])
    uff(rdmol, ignoreInterfragInteractions=False).Minimize()
    E_opt = uff(rdmol, ignoreInterfragInteractions=False).CalcEnergy()

    # Calculate the total and mean (pair-wise) interaction between between fragments
    plams_mol.properties.int = float(E_no_frag - E_frag)
    plams_mol.properties.strain = float(E_frag - (E_opt * len(mol_frag)))
    plams_mol.properties.energy = plams_mol.properties.int + plams_mol.properties.strain

    return plams_mol
