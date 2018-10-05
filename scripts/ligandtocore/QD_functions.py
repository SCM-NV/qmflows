import copy
import os
import itertools

from scm.plams import (Molecule, MoleculeError, add_to_class, Settings, AMSJob, init, finish)
from qmflows import molkit
from rdkit import Chem
from rdkit.Chem import Bond, AllChem, rdMolTransforms

import numpy as np


def create_dir(dir_name, path=os.getcwd()):
    """
    Creates a new directory if this directory does not yet exist.
    """
    dir_path = os.path.join(path, str(dir_name))
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

    return dir_path


def read_mol(folder_path, file_name, column=0, row=0, smiles_extension='.txt'):
    """
    First checks if the argument 'mol_name' is a string or a list.
    Then checks if 'mol_name' consists of .xyz/.pdb files, SMILES strings or .txt files containing
    SMILES strings.
    Returns a list of PLAMS molecules.
    """
    # Check if filename is a string or a list, returns an error if it is neither
    if isinstance(file_name, str):
        file_name = [file_name]
    if not isinstance(file_name, str) and not isinstance(file_name, list):
        raise MoleculeError("the argument 'mol_name' " + str(type(file_name)) +
                            " is not recognized as a <class 'str'> or <class 'list'>")

    # Determine the nature of filename
    input_mol_list = [os.path.join(folder_path, name) for name in file_name]
    for i, item in enumerate(input_mol_list):
        # If file_name is an .xyz file
        if file_name[i].find('.xyz') != -1:
            mol_list = [Molecule(mol) for mol in input_mol_list]

        # If file_name is a .pdb file
        elif file_name[i].find('.pdb') != -1:
            mol_list = [molkit.readpdb(mol) for mol in input_mol_list]

        # If file_name is a .mol file
        elif file_name[i].find('.mol') != -1:
            mol_list = [molkit.from_rdmol(Chem.MolFromMolFile(mol)) for mol in input_mol_list]

        # If file_name is a plain text file with smile strings
        elif file_name[i].find(smiles_extension) != -1:
            mol_list = []
            for file in input_mol_list:
                with open(file, 'r') as file_open:
                    tmp_list = file_open.read().splitlines()
                    mol_list.append(tmp_list[row:])
            mol_list = list(itertools.chain(*mol_list))
            mol_list = [line.split() for line in mol_list if bool(line)]
            mol_list = [molkit.from_smiles(mol[column]) for mol in mol_list]

        # If file_name is none of the above it is assumed to be a smile string
        else:
            mol_list = [molkit.from_smiles(mol) for mol in file_name]

    return mol_list


def set_pdb(mol, residue_name, is_core=True):
    """
    Set a number of atomic properties
    """
    # Prepare a list of letters for pdb_info.Name
    alphabet = list('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')
    alphabet = [i + j for i in alphabet for j in alphabet]

    for i, atom in enumerate(mol):
        symbol = atom.symbol + alphabet[i] + '  '

        # Add a number of properties to atom
        atom.properties.pdb_info.ResidueName = residue_name
        atom.properties.pdb_info.Occupancy = 1.0
        atom.properties.pdb_info.TempFactor = 0.0
        atom.properties.pdb_info.ResidueNumber = 1
        atom.properties.pdb_info.Name = symbol[:4]
        atom.properties.pdb_info.ChainId = 'A'

        # Changes hydrogen and carbon from heteroatom to atom
        if atom.symbol == 'H' or atom.symbol == 'C':
            atom.properties.pdb_info.IsHeteroAtom = False

        # Sets the formal atomic charge
        if is_core:
            if atom.symbol in ['Li', 'Na', 'K', 'Rb', 'Cs']:
                atom.properties.charge = 1
            elif atom.symbol in ['Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Cd', 'Pb']:
                atom.properties.charge = 2
            elif atom.symbol in ['N', 'P', 'As', 'Sb', 'Bi']:
                atom.properties.charge = -3
            elif atom.symbol in ['O', 'S', 'Se', 'Te', 'Po']:
                atom.properties.charge = -2
            elif atom.symbol in ['H', 'F', 'Cl', 'Br', 'I', 'At']:
                atom.properties.charge = -1

    return mol


def read_database(ligand_folder, database_name='ligand_database.txt'):
    """
    Open the database.
    If the database does not exist, create the database.
    """
    # Checks if database_name exists, if not creates database_name
    if not os.path.exists(os.path.join(ligand_folder, database_name)):
        with open(os.path.join(ligand_folder, database_name), 'w') as database:
            database.write('{0:6} {1:19} {2:30} {3:34} {4:}'.format(
                'Index', 'Molecular_formula', 'pdb_filename', 'pdb_opt_filename', 'SMILES_string'))

        database = [[], [], [], [], []]

    # Else opens the database
    else:
        with open(os.path.join(ligand_folder, database_name), 'r') as database:
            database = database.read().splitlines()
        database = [line.split() for line in database if line]
        database = list(zip(*database[1:]))

        if database:
            database[4] = [Chem.MolFromSmiles(smiles) for smiles in database[4]]
            database[4] = [Chem.AddHs(mol, addCoords=True) for mol in database[4]]
        else:
            database = [[], [], [], [], []]

    return database


def create_entry(ligand, opt):
    """
    Create a new entry for the database.
    """
    # compare the structures provided in ligand_list with the ones in mol_database
    ligand_smiles = Chem.MolToSmiles(Chem.RemoveHs(molkit.to_rdmol(ligand)))

    # Create the new database_entry
    a = []
    b = ligand.get_formula()
    c = 'ligand_' + b + '.pdb'
    if opt:
        d = 'ligand_' + b + '.opt.pdb'
    else:
        d = 'ligand optimization disabled'
    e = ligand_smiles

    database_entry = [a, b, c, d, e]

    return database_entry


def write_database(database_entries, ligand_folder, database, database_name='ligand_database.txt'):
    """
    Write the new database entries to the database.
    """
    # Check if any previous indices are present in database[0]
    if not database or not database[0]:
        j = -1
    else:
        j = int(database[0][-1])

    # Format the database entries
    spacing = '{0:6} {1:19} {2:30} {3:34} {4:}'
    database_entries = [spacing.format(str(j + i + 1), item[1], item[2], item[3], item[4]) for
                        i, item in enumerate(database_entries) if item]

    # Write the database entries to the database
    with open(os.path.join(ligand_folder, database_name), 'a') as database:
        for entry in database_entries:
            database.write('\n' + entry)


@add_to_class(Molecule)
def neighbors_mod(self, atom, exclude=''):
    """
    Modified PLAMS function, returns all connected atom with the exception of 'exclude'.
    Exclude can be either an atom or list of atoms.
    No atoms are excluded by default.
    """
    if not isinstance(exclude, list):
        exclude = [exclude]
    if atom.mol != self:
        raise MoleculeError('neighbors: passed atom should belong to the molecule')
    return [b.other_end(atom) for b in atom.bonds if b.other_end(atom) not in exclude]


def global_minimum(ligand, ligand_folder):
    """
    Find the glibal minimum of the ligand with RDKit UFF.
    """
    # Delete all hydrogens and create a list of all bond indices [0], bond orders [1]
    # and dihedral indices [2, 3, 4 and 5](i.e. indices of four atoms defining a dihedral angle)
    for atom in ligand:
        if atom.atnum == 1:
            ligand.delete_atom(atom)
    dihedral_list = [dihedral_index(ligand, item) for item in ligand.bonds]

    # Find the global minimum by systematically varying a select number of dihedral angles
    # All bonds are scanned that meet the following four requirements:
    # They are single bonds, non-terminal, not part of a ring and do not contain hydrogen.
    # 3 dihedral angles are checked for all abovementioned bonds
    ligand = molkit.to_rdmol(ligand)
    n_scans = 2
    for i in range(n_scans):
        if i > 0:
            ligand = Chem.AddHs(ligand, addCoords=True)
        for item in dihedral_list:
            InRing = ligand.GetBondWithIdx(item[0]).IsInRing()
            if item[2] != 'skip' and item[1] == 1.0 and not InRing:
                ligand = dihedral_scan(ligand, item)

    return molkit.from_rdmol(ligand)


def dihedral_index(ligand, bond):
    """
    Create a list of bond indices [0], bond orders [1] and dihedral indices [2, 3, 4 and 5].
    """
    # The two atoms associated with a given bond
    at1 = bond.atom1
    at2 = bond.atom2

    # The atoms bonded to at1 and at2
    at1_bonds = ligand.neighbors_mod(at1, exclude=at2)
    at2_bonds = ligand.neighbors_mod(at2, exclude=at1)

    # A list of bond indices [0], bond orders [1]
    # and dihedral indices [2, 3, 4 and 5](i.e. indices of four atoms defining a dihedral angle)
    # Only the indices of non-terminal bonds are added to this list, the rest is skipped
    if len(at1_bonds) >= 1 and len(at2_bonds) >= 1:
        bond_idx = ligand.bonds.index
        atom_idx = ligand.atoms.index
        dihedral_list = [bond_idx(bond), bond.order, atom_idx(at1_bonds[0]),
                         atom_idx(at1), atom_idx(at2), atom_idx(at2_bonds[0])]

    else:
        dihedral_list = [ligand.bonds.index(bond), bond.order, 'skip']

    return dihedral_list


def dihedral_scan(ligand, dihedral_list):
    """
    Scan a dihedral angle and find the lowest energy conformer.
    """
    # Define a number of variables and create 4 copies of the ligand
    a = dihedral_list
    uff = AllChem.UFFGetMoleculeForceField
    ligand = [copy.deepcopy(ligand) for i in range(3)]

    # Create a list of all dihedral angles for which the geometry will be optimized (rdkit uff)
    get_dihed = rdMolTransforms.GetDihedralDeg
    angle = get_dihed(ligand[0].GetConformer(), a[2], a[3], a[4], a[5])
    angle_list = [angle, angle + 120, angle - 120]

    # Optimized the geometry for all dihedral angles in angle_list
    # The geometry that yields the lowest energy is returned
    set_dihed = rdMolTransforms.SetDihedralDeg
    for i, item in enumerate(angle_list):
        set_dihed(ligand[i].GetConformer(), a[2], a[3], a[4], a[5], item)

    for item in ligand:
        uff(item).Minimize()
    energy_list = [uff(item).CalcEnergy() for item in ligand]
    minimum = energy_list.index(min(energy_list))

    return ligand[minimum]


def find_substructure(ligand):
    """
    Identify the ligand functional groups.
    """
    ligand_rdkit = molkit.to_rdmol(ligand)

    # Creates a list containing predefined functional groups, each saved as an rdkit molecule
    functional_group_list = ['[F-].[N+]', '[Cl-].[N+]', '[Br-].[N+]', '[I-].[N+]', '[H]O', '[H]S',
                             '[H]N', '[H]P']
    functional_group_list = [Chem.MolFromSmarts(smarts) for smarts in functional_group_list]

    # Searches for functional groups (defined by functional_group_list) within the ligand
    # Duplicates are removed
    get_match = ligand_rdkit.GetSubstructMatches
    matches = [get_match(mol) for mol in functional_group_list]
    matches = list(itertools.chain(*matches))

    # Remove all duplicate matches
    unique_matches = []
    for match in matches:
        if not any(match[1] in item for item in unique_matches):
            unique_matches.append(match)

    # Rotates the functional group hydrogen atom
    # In addition to returning the indices of various important ligand atoms
    ligand_list = [copy.deepcopy(ligand) for match in unique_matches]

    # Delete the hydrogen attached to the functional group and remove its index from unique_matches
    for i, ligand in enumerate(ligand_list):
        ligand.delete_atom(ligand[unique_matches[i][0] + 1])
    unique_matches = [item[1] for item in unique_matches]

    if not ligand_list:
        print('No functional groups were found for ' + str(ligand.get_formula()))

    return ligand_list, unique_matches


def rotate_ligand(core, ligand, core_index, ligand_index, i):
    """
    Connects two molecules by alligning the vectors of two bonds.
    """
    ligand = copy.deepcopy(ligand)

    # Defines first atom on coordinate list (hydrogen),
    # The atom connected to it and vector representing bond between them
    core_at1 = core[core_index]         # core dummy atom
    core_at2 = core[-1]                 # core center of mass
    core_vector = core_at1.vector_to(core_at2)
    lig_at1 = ligand[ligand_index + 1]  # ligand heteroatom
    lig_at2 = ligand[-1]                # ligand center of mass
    lig_vector = lig_at2.vector_to(lig_at1)

    # Rotation of ligand - aligning the ligand and core vectors
    rotmat = rotation_matrix(lig_vector, core_vector)
    ligand.rotate(rotmat)
    ligand.translate(lig_at1.vector_to(core_at1))

    # Translation of the ligand
    hc_vec = lig_at1.vector_to(core_at1)
    ligand.translate(hc_vec)

    # Deletes the core dummy atom and ligand center of mass
    ligand.delete_atom(lig_at2)
    core.delete_atom(core_at1)

    # Update the residue numbers
    for atom in ligand:
        atom.properties.pdb_info.ResidueNumber = i + 2

    # Check if the ligand heteroatom has a charge assigned, assigns a charge if not
    if not lig_at1.properties.charge or lig_at1.properties.charge == 0:
        lig_at1.properties.charge = -1

    return ligand, lig_at1


def rotation_matrix(vec1, vec2):
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
    qd = copy.deepcopy(core)

    # Create a list of ligand atoms and intraligand bonds
    ligand_bonds = np.concatenate([ligand.bonds for ligand in ligand_list])
    ligand_atoms = np.concatenate(ligand_list)

    # Combined the ligand bond and atom list with the core
    for atom in ligand_atoms:
        qd.add_atom(atom)
    for bond in ligand_bonds:
        qd.add_bond(bond)

    return qd


def run_ams_job(qd, pdb_name, qd_folder, qd_indices, maxiter=1000, opt=False):
    """
    Converts PLAMS connectivity into adf .run script connectivity and optimizes the qd.
    """
    # Create a list of aromatic bond indices
    qd_rdkit = molkit.to_rdmol(qd)
    aromatic = [Bond.GetIsAromatic(bond) for bond in qd_rdkit.GetBonds()]
    aromatic = [i for i, item in enumerate(aromatic) if item]
    
    # Create a connectivity list; aromatic bonds get a bond order of 1.5
    at1 = [qd.atoms.index(bond.atom1) + 1 for bond in qd.bonds]
    at2 = [qd.atoms.index(bond.atom2) + 1 for bond in qd.bonds]
    bonds = [bond.order for bond in qd.bonds]
    for i, bond in enumerate(qd.bonds):
        if i in aromatic:
            bonds[i] = 1.5
    bonds = [str(at1[i]) + ' ' + str(at2[i]) + ' ' + str(bond) for i, bond in enumerate(bonds)]

    # General AMS settings
    s = Settings()
    s.input.ams.Task = 'GeometryOptimization'
    s.input.ams.Constraints.Atom = qd_indices
    s.input.ams.System.BondOrders._1 = bonds
    s.input.ams.GeometryOptimization.MaxIterations = maxiter
    s.input.ams.Properties.Gradients = 'Yes'

    # Settings specific to UFF
    s.input.uff.Library = 'UFF'

    # Run the job if opt = True
    init(path=qd_folder, folder=pdb_name)
    job = AMSJob(molecule=qd, settings=s, name=pdb_name)
    if opt:
        results = job.run()
        output_mol = results.get_main_molecule()
        for i, atom in enumerate(qd):
            atom.move_to(output_mol[i + 1])

        qd.write(os.path.join(qd_folder, pdb_name + '.opt.xyz'))
        molkit.writepdb(qd, os.path.join(qd_folder, pdb_name + '.opt.pdb'))
        print('Optimized core + ligands:\t\t' + pdb_name + '.opt.pdb')

        return qd

    finish()
