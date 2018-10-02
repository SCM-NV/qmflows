from scm.plams import *
from qmflows import molkit
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms, rdForceFieldHelpers
from rdkit.Chem import Atom as At
import numpy as np
import copy
import os
import re


def check_database(ligand_list, ligand_folder, ligand_opt, database_name='ligand_database.txt'):
    """
    Check if the input ligand is already in the database.
    If yet: append the .pdb file.
    If not: create a new entry in the database
    """
    if not os.path.exists(os.path.join(ligand_folder, database_name)):
        with open(os.path.join(ligand_folder, database_name), 'w') as database:
                database.write("{0:6} {1:19} {2:30} {3:34} {4:}".format('Index', 'Molecular_formula', 'pdb_filename', 'pdb_opt_filename', 'SMILES_string'))
    
    with open(os.path.join(ligand_folder, database_name), 'r') as database:
        database = database.read().splitlines()
    database = [line.split() for line in database if line]
    database = np.transpose(database)

    index_database = database[0]
    pdb_opt_database = database[3]
    smiles_database = database[4]
    mol_database = [Chem.MolFromSmiles(smiles) for smiles in smiles_database[2:]]
    mol_database = [Chem.AddHs(mol, addCoords=True) for mol in mol_database if mol]

    if not isinstance(index_database[-1], int):
        index_database = np.append(index_database, -1)

    smiles_list = [Chem.MolToSmiles(Chem.RemoveHs(molkit.to_rdmol(ligand))) for ligand in ligand_list]
    matches = [[molkit.to_rdmol(ligand).HasSubstructMatch(mol) for mol in mol_database] for ligand in ligand_list]

    for i,item in enumerate(matches):
        if any(item):
            formula = ligand_list[i].get_formula()
            print('{0:12} present in '.format(formula) + database_name + ', appending ligand_' + formula + '.opt.pdb')
    
            index = np.where(smiles_database == smiles_list[i])[0]
            read_pdb = molkit.readpdb(os.path.join(ligand_folder, pdb_opt_database[index][0]))
            if ligand_opt:
                ligand_list[i] = [read_pdb, False]
            else:
                ligand_list[i] = read_pdb

        else:
            a = str(int(index_database[-1]) + i + 1)
            b = molkit.from_rdmol(ligand_list[i]).get_formula()
            c = 'ligand_' + b + '.pdb'
            d = 'ligand_' + b + '.opt.pdb'
            e = smiles_list[i]
            entry = "{0:6} {1:19} {2:30} {3:34} {4:}".format(a, b, c, d, e)

            with open(os.path.join(ligand_folder, database_name), 'a') as database:
                database.write('\n' + entry)

            if ligand_opt:
                ligand_list[i] = [ligand_list[i], True]
            else:
                ligand_list[i] = ligand_list[i]

    return ligand_list


def create_dir(dir_name, path=os.getcwd()):
    """
    Creates a new directory if this directory does not yet exist
    """
    dir_path = os.path.join(path, str(dir_name))
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

    return dir_path



def read_mol(folder_path, file_name, smiles_column=0, smiles_extension='.txt'):
    """
    first checks if the argument 'mol_name' is a string or a list
    then checks if 'mol_name' consists of .xyz/.pdb files, SMILES strings or .txt files containing SMILES strings
    returns a list of PLAMS molecules
    """
    # check if filename is a string or a list, returns an error if it is neither
    if isinstance(file_name, str):
        file_name = [file_name]
    if not isinstance(file_name, str) and not isinstance(file_name, list):
        raise MoleculeError("the argument 'mol_name' " + str(type(file_name)) + " is not recognized as a <class 'str'> or <class 'list'>")

    # determine the nature of filename
    input_mol_list = [os.path.join(folder_path, name) for name in file_name]
    for i,item in enumerate(input_mol_list):
        # if file_name is an .xyz file 
        if file_name[i].find('.xyz') != -1:
            mol_list = [Molecule(mol) for mol in input_mol_list]
            [mol.guess_bonds for mol in mol_list]
        
        # if file_name is a .pdb file
        elif file_name[0].find('.pdb') != -1:
            mol_list = [molkit.readpdb(mol) for mol in input_mol_list]

        # if file_name is a .mol file
        elif file_name.find('.mol') != -1:
            mol_list = [molkit.from_rdmol(Chem.MolFromMolFile(mol)) for mol in input_mol_list]

        # if file_name is a plain text file with smile strings
        elif file_name[0].find(smiles_extension) != -1:
            mol_list = []
            for file in input_mol_list:
                with open(file, 'r') as file_open:
                    mol_list.append(file_open.read())
            mol_list = mol_list.splitlines()
            mol_list = [line.split() for line in mol_list if bool(line)]
            mol_list = [molkit.from_smiles(mol[smiles_column]) for mol in mol_list]

        # if file_name is none of the above it is assumed to be a smile string
        else:
            mol_list = [molkit.from_smiles(mol) for mol in file_name]

    return mol_list


@add_to_class(Molecule)
def neighbors_mod(self, atom, exclude=[]):
    """
    Modified PLAMS function, returns all connected atom with the exception of 'exclude'
    Exclude can be either an atom or list of atoms
    No atoms are excluded by default
    """
    if not isinstance(exclude, list):
        exclude = [exclude]
    if atom.mol != self:
        raise MoleculeError('neighbors: passed atom should belong to the molecule')
    atom_list = [b.other_end(atom) for b in atom.bonds if b.other_end(atom) not in exclude]
    
    return atom_list



def global_minimum(ligand, ligand_folder):
    """
    optimize the geometry of the ligand with uff
    """
    ligand_name = 'ligand_' + ligand.get_formula()
    molkit.writepdb(ligand, os.path.join(ligand_folder, ligand_name + '.pdb'))
    print('Ligand:\t\t\t\t' + str(ligand_name) + '.pdb')
    
    # delete all hydrogens and create a list of all bond indices [0], bond orders [1] and dihedral indices [2, 3, 4 & 5](i.e. the indices of the four atoms defining a dihedral angle)
    [ligand.delete_atom(atom) for atom in ligand.atoms if atom.atnum == 1]
    dihedral_list = [dihedral_index(ligand, item) for item in ligand.bonds]

    # find the global minimum by systematically varying a select number of dihedral angles in a serial manner; this process is carried out twice.
    # all bonds are scanned that meet the following four requirements: they are single bonds, non-terminal, not part of a ring and do not contain hydrogen.
    # 3 dihedral angles are checked for all abovementioned bonds
    ligand = molkit.to_rdmol(ligand)
    n_scans = 1
    for i in range(n_scans):
        for item in dihedral_list:
            if item[2] != 'skip' and item[1] == 1.0 and not ligand.GetBondWithIdx(item[0]).IsInRing():     
                ligand = dihedral_scan(ligand, item)
    
    # reatatch all hydrogens and optimize the resulting structure
    ligand = Chem.AddHs(ligand, addCoords=True)
    uff = AllChem.UFFGetMoleculeForceField
    uff(ligand).Minimize()
    
    # export the ligand to a .pdb and .xyz file
    ligand = molkit.from_rdmol(ligand)
    molkit.writepdb(ligand, os.path.join(ligand_folder, ligand_name + '.opt.pdb'))
    print('Optimized ligand:\t\t' + str(ligand_name) + '.opt.pdb')

    return ligand


def dihedral_index(ligand, bond):
    """
    create a list of bond indices [0], bond orders [1] and dihedral indices [2, 3, 4 & 5]
    """
    # the two atoms associated with a given bond
    at1 = bond.atom1
    at2 = bond.atom2
    
    # the atoms bonded to at1 and at2
    at1_bonds = ligand.neighbors_mod(at1, exclude=at2)
    at2_bonds = ligand.neighbors_mod(at2, exclude=at1)
    
    # a list of bond indices [0], bond orders [1] and dihedral indices [2, 3, 4 & 5](i.e. the indices of the four atoms defining a dihedral angle)
    # only the indices of non-terminal bonds are added to this list, the rest is skipped
    if len(at1_bonds) >= 1 and len(at2_bonds) >= 1:
        dihedral_list = [ligand.bonds.index(bond), bond.order, ligand.atoms.index(at1_bonds[0]), ligand.atoms.index(at1), ligand.atoms.index(at2), ligand.atoms.index(at2_bonds[0])]
    else:
        dihedral_list = [ligand.bonds.index(bond), bond.order, 'skip']
        
    return dihedral_list


def dihedral_scan(ligand, dihedral_list):
    """
    Scan a dihedral angle and find the lowest energy conformer
    """
    # define a number of variables and create 4 copies of the ligand
    a = dihedral_list
    uff = AllChem.UFFGetMoleculeForceField
    ligand = [copy.deepcopy(ligand), copy.deepcopy(ligand), copy.deepcopy(ligand)]

    # create a list of all dihedral angles for which an uff geometry optimization will be carried out
    angle = rdMolTransforms.GetDihedralDeg(ligand[0].GetConformer(), a[2], a[3], a[4], a[5])
    angle_list = [angle, angle + 120, angle - 120]
    
    # optimized the geometry for all dihedral angles in angle_list; the one geometry that yields the lowest energy is returned
    [rdMolTransforms.SetDihedralDeg(ligand[i].GetConformer(), a[2], a[3], a[4], a[5], item) for i,item in enumerate(angle_list)]
    [uff(item).Minimize() for item in ligand]
    energy_list = [uff(item).CalcEnergy() for item in ligand]
    minimum = energy_list.index(min(energy_list))
    
    return ligand[minimum]



def find_substructure(ligand):
    """
    Identify the ligand functional groups
    """
    ligand_rdkit = molkit.to_rdmol(ligand)
    
    # creates a list containing predefined functional groups, each saved as an rdkit molecule
    functional_group_list = []
    functional_group_list.append(Chem.MolFromSmarts('[F-].[N+]CC'))    # ammonium halide
    functional_group_list.append(Chem.MolFromSmarts('[Cl-].[N+]CC'))    # ammonium halide
    functional_group_list.append(Chem.MolFromSmarts('[Br-].[N+]CC'))    # ammonium halide
    functional_group_list.append(Chem.MolFromSmarts('[I-].[N+]CC'))    # ammonium halide
    functional_group_list.append(Chem.MolFromSmarts('[H]OCC'))          # hydroxides
    functional_group_list.append(Chem.MolFromSmarts('[H]OCc'))          # benzylic hydroxides
    functional_group_list.append(Chem.MolFromSmarts('[H]SCC'))          # sulfides
    functional_group_list.append(Chem.MolFromSmarts('[H]SCc'))          # benzylic sulfides
    
    # searches for functional groups (defined by functional_group_list) within the ligand; duplicates are removed 
    matches = [list(ligand_rdkit.GetSubstructMatches(item)) for item in functional_group_list]
    matches = [j for i in matches for j in i]
    matches = [item for i,item in enumerate(matches) if item[1] != matches[i - 1][1] or i == 0]
    
    # rotates the functional group hydrogen atom in addition to returning the indices of various important ligand atoms
    ligand_copies = [find_substructure_part2(ligand, item) for item in matches if item]
    
    if not ligand_copies:
        print('No functional groups were found for ' + str(ligand.get_formula()))

    return ligand_copies


def find_substructure_part2(ligand, matches):
    """
    Identify import atoms within a functional group and manipulate bond angles associated with aforementioned atoms
    """
    # the indices of various import ligand atoms
    H = matches[0] + 1
    X = matches[1] + 1
    Cb = matches[3] + 1
    
    # various important ligand bond lengths
    Cb_X = ligand[Cb].distance_to(ligand[X])
    X_H = ligand[X].distance_to(ligand[H])

    # orientation of the H atom: sets the angle between H, X and Cb to 180 degrees
    vec = np.asarray(ligand[Cb].vector_to(ligand[X])) * ((Cb_X + X_H) / Cb_X)
    ligand[H].move_to(ligand[Cb])
    ligand[H].translate(vec)

    matches = [index + 1 for index in matches]

    return [ligand] + matches


def rotate_ligand(core, ligand, index):
    """
    Connects two molecules by alligning the vectors of two bonds
    """
    # Define the indices of important atoms
    H = ligand[1]
    X = ligand[2]
    Cb = ligand[4]
    ligand = copy.deepcopy(ligand[0])

    # Defines first atom on coordinate list (hydrogen), atom connected to it and vector representing bond between them
    core_h = core[index]
    core_other = core[-1]
    core_vector = core_h.vector_to(core_other)
    lig_h = ligand[H]
    lig_other = ligand[-1]
    lig_vector = lig_other.vector_to(ligand[X])

    # Rotation of ligand - aligning diresctions of two vectors
    rotmat = rotation_matrix(lig_vector, core_vector)
    ligand.rotate(rotmat)
    ligand.translate(ligand[X].vector_to(core_h))

    # Translation of the ligand
    hc_vec = ligand[X].vector_to(core_h)
    ligand.translate(hc_vec)

    # Deletes atom in ligand
    ligand.delete_atom(lig_other)
    ligand.delete_atom(lig_h)
    
    # Deletes atom in ligand
    core.delete_atom(core_h)

    return ligand, ligand[X]


# Function for matrix rotation
def rotation_matrix(vec1, vec2):
    """
    Calculates the rotation matrix rotating *vec1* to *vec2*. Vectors can be any containers with 3 numerical values. They don't need to be normalized. Returns 3x3 numpy array.
    """
    a = np.array(vec1) / np.linalg.norm(vec1)
    b = np.array(vec2) / np.linalg.norm(vec2)
    v1, v2, v3 = np.cross(a, b)
    M = np.array([[0, -v3, v2], [v3, 0, -v1], [-v2, v1, 0]])
    
    return (np.identity(3) + M + np.dot(M, M)/(1+np.dot(a, b)))


def combine_core_ligand(core, ligand_list):
    """
    combine the rotated ligands with the core, creating a bond bewteen the core and ligand in the process
    """
    core_ligand = copy.deepcopy(core)
    
    # create a list of ligand atoms and intraligand bonds
    ligand_bonds = np.concatenate([ligand.bonds for ligand in ligand_list])
    ligand_atoms = np.concatenate(ligand_list)

    # Combined the ligand bond & atom list with the core
    [core_ligand.add_atom(atom) for atom in ligand_atoms]
    [core_ligand.add_bond(bond) for bond in ligand_bonds]

    return core_ligand


def optimize_core_ligand(core_ligand, core_ligand_indices, maxiter=200):
    """
    optimize the combined core and ligands with the core frozen
    """
    uff = AllChem.UFFGetMoleculeForceField(core_ligand, ignoreInterfragInteractions=False)
    [[uff.AddFixedPoint(index) for index in index_list] for index_list in core_ligand_indices]
    uff.Initialize()
    
    if not rdForceFieldHelpers.UFFHasAllMoleculeParams(core_ligand):
        print('Warning: uff parameters unavailable for one or more atoms, possibly due to incorrect valency or formal atomic charges')
    
    print('\nCore + ligands optimization:')
    for i in range(int(maxiter / 10)):
        uff.Minimize(maxIts = 10)
        print(str((i + 1) * 10) + '/' + str(maxiter) + '\tEnergy:\t' + "{0:.4f}".format(uff.CalcEnergy()) + '\t|Grad|: ' + "{0:.4f}".format(np.linalg.norm(uff.CalcGrad())))   

    return core_ligand


def prepare_pdb(core_ligand, core, ligand, core_ligand_indices):  
    """
    add residue names and formal atomic charges to the molecule
    """
    # define the number of atoms in the core and the number of ligands
    len_core = len(core.atoms)
    len_ligand = len(ligand[0].atoms) - 1
    n_ligands = int((len(core_ligand.GetAtoms()) - len_core) / len_ligand)
    
    # add formal atomic charges to the carboxylate oxygens, Se and Cd 
    charge_list = [[1, 'Li', 'Na', 'K', 'Rb', 'Cs'][2, 'Be', 'Mg', 'Ca', 'Sr', 'Ba'][-2, 'O', 'S', 'Se', 'Te', 'Po'][-1, 'F', 'Cl', 'Br', 'I', 'At']]
    charge = [+1, +1, -1, +2]
    [[At.SetFormalCharge(core_ligand.GetAtoms()[index], charge[i]) for index in index_list] for i,index_list in enumerate(core_ligand_indices)]
    
    # Assigns the core to the "QD " residual
    pdb = Chem.MolToPDBBlock(core_ligand)
    pdb = pdb.splitlines()
    pdb_residue = [item.replace("UNL", "QD ") for i, item in enumerate(pdb) if i + 1 <= len_core]
    
    # assign the ligands to the "###" residues, where ### is an integer consisting of three numbers
    offset = len_core
    for i in range(n_ligands):
        pdb_residue += [pdb[j + offset].replace("UNL", str(i).zfill(3)) for j in range(len_ligand)]
        offset += len_ligand
   
    # export the new residue names back to the molecule
    pdb_residue += pdb[offset:-1]
    pdb_residue.append('END')
    pdb_residue = '\n'.join(pdb_residue)
    core_ligand = Chem.MolFromPDBBlock(pdb_residue, removeHs = False, proximityBonding = False)
    
    return core_ligand


def run_ams_job(core_ligand):
    """
    converts PLAMS connectivity into adf .run script connectivity
    """
    at1 = [core_ligand.atoms.index(bond.atom1) + 1 for bond in core_ligand.bonds]
    at2 = [core_ligand.atoms.index(bond.atom2) + 1 for bond in core_ligand.bonds]
    bonds = [bond.order for bond in core_ligand.bonds]
    bonds = [str(at1[i]) + ' ' + str(at2[i]) + ' ' + str(bond) for i,bond in enumerate(bonds)]

    s = Settings()
    s.input.ams.Task = 'GeometryOptimization'
    s.input.ams.Constraints.Atom = [core_ligand.atoms.index(atom) for atom in core_ligand if atom.atnum == 7 or atom.atnum == 35 or atom.atnum == 82 or atom.atnum == 55]
    s.input.ams.System.BondOrders._1 = bonds
    s.input.ams.GeometryOptimization.MaxIterations = 100
    s.input.ams.Properties.Gradients = 'Yes'
    s.input.uff.Library = 'UFF'
    
    j = AMSJob(molecule=core_ligand, settings=s)
    
    with open('run.run', 'w') as run:
        run.write(j.get_input())


def update_adf_pdb(core_ligand):
    """
    update a .pdb file using coordinates provided by an .xyz file
    """
    with open('core_ligand/ams.1258.xyz', 'r') as cube:
        xyz = cube.read().splitlines()
    xyz = [item.split() for i,item in enumerate(xyz) if i > 1]
    xyz = [[(subitem)[:6] for i,subitem in enumerate(item) if i > 0] for item in xyz]
    xyz = ['     1      ' + item[0] + '  ' + item[1] + '  ' + item[2] + '  1.00  0.00           ' for item in xyz]
    
    with open('core_Br450Cs157Pb125__&__ligand_C26H56N1_@_N13.pdb', 'r') as cube_read:
        cube_read = cube_read.read().splitlines()
        cube_read = [re.sub('     1      .*?  1.00  0.00           ', xyz[i], item) for i, item in enumerate(cube_read) if i < 4301]
        with open('core_Br450Cs157Pb125__&__ligand_C26H56N1_@_N13.opt.pdb', 'w') as cube_write:
            for item in cube_read:
                cube_write.write("%s\n" % item)
