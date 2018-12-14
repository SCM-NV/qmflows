__all__ = ['get_topology_dict', 'dissociate_ligand']

from .qd_functions import qd_int


def get_topology_dict(mol, dist=4.5):
    """
    Returns the topology of all ligands on a quantum dot (vertice, edge or face) based on
        the amount of neighbouring ligands.
    mol <plams.Molecule>: A PLAMS molecule.
    dist <float>: The maximum distance between a pair of ligands to be considered neighbours.
    return <dict>: A dictionary with residue numbers as keys and
        the topology of the ligand (vertice, edge or face) as value.
    """
    dist_dict = {0: 'None', 1: 'None', 2: 'vertice', 3: 'edge', 4: 'face'}
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
    return topology_dict


def dissociate_ligand(plams_mol, n=2, res_old=False):
    """
    Create all possible combinations of quantum dots by removing n ligands.
    mol <plams.Molecule>: A PLAMS molecule.
    n <int>: The number of to be removed ligands.
    res_old False or <list>: A list of residue numbers of previously removed residues.
    return <list>[<plams.Molecule>]: A list of PLAMS molecules with between 1 and n ligands removed.
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


def get_residue_dict(mol):
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
    return residue_dict
