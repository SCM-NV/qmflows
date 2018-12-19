__all__ = ['get_topology_dict', 'dissociate_ligand', 'diss_list_to_pd']

import copy
import itertools
import numpy as np
import pandas as pd

from .qd_functions import qd_int


def get_topology_dict(mol, dist=4.5):
    """
    Returns the topology of all ligands on a quantum dot (vertice, edge or face) based on
        the amount of neighbouring ligands.
    mol <plams.Molecule>: A PLAMS molecule.
    dist <float>: The maximum distance between a pair of ligands to be considered neighbours.
    return <dict>: A dictionary with residue numbers of ligands as keys and
        the topology of the ligand (vertice, edge, face or None) as value.
    """
    dist_dict = {0: None, 1: None, 2: 'vertice', 3: 'edge', 4: 'face'}
    topology_dict = {}
    residue_dict = get_residue_dict(mol)
    del residue_dict[1]
    idx = list(residue_dict.keys())[-1]
    idx = residue_dict[idx].index(mol[mol.properties.indices[-1]])
    for residue in residue_dict:
        at = residue_dict[residue][idx]
        dist_list = [at.distance_to(residue_dict[residue][idx]) for residue in residue_dict if
                     at.distance_to(residue_dict[residue][idx]) <= dist]
        topology_dict[residue] = dist_dict[len(dist_list) - 1]
    topology_dict[None] = None
    return topology_dict


def dissociate_ligand(plams_mol, n=2, res_old=False):
    """
    Create all possible combinations of quantum dots by removing n ligands.
    mol <plams.Molecule>: A PLAMS molecule.
    n <int>: The number of to be removed ligands.
    res_old False or <list>: A list of residue numbers of previously removed residues.
    return <list>[<plams.Molecule>]: A list of PLAMS molecules which have between
        1 and n ligands removed.
    """
    residue_dict = get_residue_dict(plams_mol)
    del residue_dict[1]
    plams_mol.set_atoms_id()
    mol_list, residue_list = zip(*(delete_ligand(plams_mol.copy(), residue_dict[residue]) for
                                   residue in residue_dict))
    mol_list, residue_list = list(mol_list), list(residue_list)
    for mol in mol_list:
        qd_int(mol)
    if res_old:
        residue_list = [res_old.copy() + res for res in residue_list]
    if n > 1:
        for mol, res in zip(reversed(mol_list), reversed(residue_list)):
            mol_new, res_new = dissociate_ligand(mol, n=n-1, res_old=res)
            residue_list += res_new
            mol_list += mol_new
    return mol_list, residue_list


def delete_ligand(mol, atom_list):
    """
    Delete all atoms in a molecule whose id attribute matches the one provided by atom_list.
    mol <plams.Molecule>: A PLAMS molecule.
    atom <list>[<plams.Atom>]: A list of PLAMS atoms with the id attribute.
    return <plams.Molecule>, <list>[<int>]: A copy of mol with all atoms from atom_list removed and
        the residue number of the first atom removed.
    """
    for atom in reversed(atom_list):
        mol.delete_atom(mol[atom.id])
    return mol, [atom_list[0].properties.pdb_info.ResidueNumber]


def get_residue_dict(mol, mark=True):
    """
    Creates a dictionary of atom residue numbers and their corresponding atoms.
    mol <plams.Molecule>: A PLAMS molecule with the properties.pdb_info.ResidueNumber attribute.
    return <dict>: A dictionary with residue numbers as keys and
        a list of corresponding PLAMS atoms as values.
    """
    residue_dict = {}
    for atom in mol:
        try:
            residue_dict[atom.properties.pdb_info.ResidueNumber].append(atom)
        except KeyError:
            residue_dict[atom.properties.pdb_info.ResidueNumber] = [atom]
    if mark:
        for i in mol.properties.indices:
            mol[i].properties.mark = True
        mol.properties.mark = []
    return residue_dict


def diss_list_to_pd(diss_list, residue_list, top_dict):
    """
    Converts specific entries from diss_list and residue_list to a pandas dataframe.
    """
    # Make all sublists of equal length
    max_depth = len(residue_list[-1])
    for res in residue_list:
        len_res = len(res)
        if len_res < max_depth:
            res.extend([None for i in range(max_depth - len_res)])

    # Create a dictionary containg E, Eint and Estrain
    gen1 = ((qd.properties.E, qd.properties.Eint, qd.properties.Estrain) for qd in diss_list)
    dict1 = dict(zip(('E', 'Eint', 'Estrain'), zip(*gen1)))

    # Create a dictionary containing residue numbers and toplogy
    gen2 = ((itertools.chain(*(res, [top_dict[i] for i in res]))) for res in residue_list)
    keys = ('Residue numbers', 'Topology')
    keys = [key + '_' + str(i + 1) for key in keys for i in range(max_depth)]
    dict2 = dict(zip(keys, zip(*gen2)))

    dict1.update(dict2)
    return pd.DataFrame(dict1)


def average_energy(df):
    """
    Calculate the average of E, Eint & Estrain for the removal of each ligand with a given topology.
    """
    topology = [item for item in df.keys() if 'Topology' in item]
    topology = np.array(df[topology])

    vertice = np.array([(i, j, k) for i, j, k, top in
                        zip(df['E'], df['Eint'], df['Estrain'], topology) if 'vertice' in top])
    edge = np.array([(i, j, k) for i, j, k, top in
                     zip(df['E'], df['Eint'], df['Estrain'], topology) if 'edge' in top])
    face = np.array([(i, j, k) for i, j, k, top in
                     zip(df['E'], df['Eint'], df['Estrain'], topology) if 'face' in top])

    ret = pd.DataFrame({'vertice': np.mean(vertice, axis=0),
                        'edge': np.mean(edge, axis=0),
                        'face': np.mean(face, axis=0),
                        'Energy': ('E', 'Eint', 'Estrain')})
    ret.set_index('Energy')

    return ret


def delete_ligand_cd(plams_mol, atom_list, get_point=False):
    mol = plams_mol.copy()
    mol.properties = copy.deepcopy(plams_mol.properties)
    for atom in reversed(atom_list):
        mol.delete_atom(mol[atom.id])
        if atom.properties.mark:
            mol.properties.mark.append(atom.properties.pdb_info.ResidueNumber)
            if get_point:
                point = atom.coords
    if get_point:
        return mol, point
    return mol


def delete_cd(xyz_array, mol, point):
    """
    Create a copy of mol and delete the two atoms closest to the atom in mol.properties.mark.
    xyz_array <np.ndarray>: A n*3 numpy array.
    mol <plams.Molecule>: A PLAMS molecule with the mol.properties.mark attribute.
    return <plams.Molecule>, <plams.Molecule>: Two copies of mol with the first and second
        closest atom to mol.properties.mark removed, respectively.
    """
    idx1, idx2 = closest_atoms(xyz_array, point)
    mol1, mol2 = mol.copy(), mol.copy()
    mol1.properties, mol2.properties = copy.deepcopy(mol.properties), copy.deepcopy(mol.properties)
    mol1.properties.mark.append(idx1)
    mol1.delete_atom(mol1[idx1])
    mol2.properties.mark.append(idx2)
    mol2.delete_atom(mol2[idx2])
    return mol1, mol2


def closest_atoms(xyz_array, point, n=2):
    """
    Returns the indices of the n rows in xyz_array closest to point (i.e. smallest norm).
    xyz_array <np.ndarray>: A n*3 numpy array.
    point <tuple>: A 3-tuple consisting of floats.
    n <int>: How many of the smallest norms should be returned.
    return <list>[<int>]: The indices of the n rows yielding the smallest norms.
    """
    dist_array = np.linalg.norm(np.array(point) - xyz_array, axis=1)
    return [int(np.argpartition(dist_array, i)[i] + 1) for i in range(n)]


def dissociate_ligand_cd(plams_mol):
    """
    Create all possible combinations of quantum dots by removing 2 ligands and one of the two Cd
        atoms closest to the first ligand.
    mol <plams.Molecule>: A PLAMS molecule.
    return <generator>: A generator that yields n*2*(n-1) molecules.
    """
    plams_mol.set_atoms_id()
    from_iter = itertools.chain.from_iterable
    res_dict = get_residue_dict(plams_mol)
    core_array = np.array([at.coords for at in res_dict[1]])
    del res_dict[1]

    def dissociate_ligand_cd_2(plams_mol):
        plams_mol.set_atoms_id()
        res_dict = get_residue_dict(plams_mol, mark=False)
        del res_dict[1]
        return (delete_ligand_cd(plams_mol, res_dict[key]).properties for key in res_dict)

    # Remove the first ligand: n possibilities
    mol_gen, points = zip(*(delete_ligand_cd(plams_mol, res_dict[key], get_point=True) for
                            key in res_dict))

    # Remove one of the two closest cadmium atoms: n*2 possibilities
    mol_gen = from_iter(delete_cd(core_array, mol, point) for mol, point in zip(mol_gen, points))

    # Remove the second ligand: n*2*(n-1) posibilities
    mol_gen = from_iter(dissociate_ligand_cd_2(mol) for mol in mol_gen)
    prop_list = [mol.properties.mark for mol in mol_gen]

    plams_mol.unset_atoms_id()
    del plams_mol.properties.mark

    return prop_list
