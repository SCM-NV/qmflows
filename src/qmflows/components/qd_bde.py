__all__ = ['get_bde']

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from scm.plams import (Molecule, Atom)
from scm.plams.core.functions import (init, finish, config)
from scm.plams.interfaces.adfsuite.ams import AMSJob

from .qd_functions import to_atnum
from .qd_ligand_rotate import rot_mol_angle
from .qd_ams import ams_job_mopac_opt, ams_job_mopac_sp, ams_job_uff_opt
from .qd_functions import merge_mol
from .qd_dissociate import dissociate_ligand

from ..templates.templates import get_template


def init_bde(mol, job1=None, job2=None, s1=None, s2=None):
    lig = get_cdx2(mol)
    core = dissociate_ligand(mol)

    res_list1, idx_list, res_list2 = zip(*[mol.properties.mark for mol in core])
    df = pd.DataFrame({'Ligand Residue Num #1': res_list1,
                       'Cd Topology': get_topology(mol, idx_list),
                       'Ligand Residue Num #2': res_list2})
    df['dE'] = get_bde_dE(mol, lig, core, job=job1, s=s1)
    df['dG'], df['ddG'] = get_bde_ddG(mol, lig, core, job=job2, s=s2)

    return df


def get_bde_dE(tot, lig, core, job=None, s=None):
    """ Calculate the bond dissociation energy: dE = dE(mopac) + (dG(uff) - dE(uff))
    """
    init(path=tot.properties.path, folder='BDE_dE')
    config.default_jobmanager.settings.hashing = None

    if job is None and s is None:
        job = AMSJob
        s = get_template('qd.json')['MOPAC']
    elif job is None or s is None:
        finish()
        raise TypeError('job & s should neither or both be None')

    # Perform single points
    tot.job_single_point(job, s)
    for mol in core:
        mol.job_single_point(job, s)
    lig.job_geometry_opt(job, s)

    # Extract total energies
    E_lig = lig.properties.energy.E
    E_core = np.array([mol.properties.energy.E for mol in core])
    E_tot = tot.properties.energy.E

    # Calculate and return dE
    dE = (E_lig + E_core) - E_tot
    finish()

    return dE


def get_bde_ddG(tot, lig, core, job=None, s=None):
    """ Calculate the bond dissociation energy: dE = dE(mopac) + (dG(uff) - dE(uff))
    """
    init(path=tot.properties.path, folder='BDE_ddG')
    config.default_jobmanager.settings.hashing = None

    if job is None and s is None:
        job = AMSJob
        s = get_template('qd.json')['UFF']
    elif job is None or s is None:
        finish()
        raise TypeError('job & s should neither or both be None')

    # Perform a constrained geometry optimizations + frequency analyses
    s.input.ams.Constraints.Atom = lig.properties.indices
    lig.job_freq(job, s)
    for mol in core:
        s.input.ams.Constraints.Atom = mol.properties.indices
        mol.job_freq(job, s)
    s.input.ams.Constraints.Atom = mol.properties.indices
    tot.job_freq(job, s)

    # Extract total Gibbs free energies
    G_lig = lig.properties.energy.G
    G_core = np.array([mol.properties.energy.G for mol in core])
    G_tot = tot.properties.energy.G

    # Extract total energies
    E_lig = lig.properties.energy.E
    E_core = np.array([mol.properties.energy.E for mol in core])
    E_tot = tot.properties.energy.E

    # Calculate and return dG and ddG
    dG = (G_lig + G_core) - G_tot
    ddG = dG - ((E_lig + E_core) - E_tot)
    finish()

    return np.array([dG, ddG])


def get_cdx2(mol, ion='Cd'):
    """ Takes a quantum dot with ligands (X) and an ion (Y) and turns it into YX2.
    Returns the total energy of YX2 at the MOPAC level of theory. """
    def get_anchor(mol):
        for i, at in enumerate(mol.atoms):
            if at.properties.anchor:
                return i, at

    def get_ligand(mol):
        at_list = []
        res = mol.atoms[-1].properties.pdb_info.ResidueNumber
        for at in reversed(mol.atoms):
            if at.properties.pdb_info.ResidueNumber == res:
                at_list.append(at)
            else:
                ret = Molecule()
                ret.atoms = at_list
                return ret.copy()

    # Translate the ligands to their final position
    lig1 = get_ligand(mol)
    lig2 = lig1.copy()
    idx1, anchor1 = get_anchor(lig1)
    idx2, anchor2 = get_anchor(lig2)
    target = np.array([2.2, 0.0, 0.0])
    lig1.translate(anchor1.vector_to(target))
    lig2.translate(anchor2.vector_to(-target))

    # Define vectors for the ligand rotation
    vec1_1 = np.array(anchor1.vector_to(lig1.get_center_of_mass()))
    vec2_1 = -1 * np.array(anchor1.vector_to(np.zeros(3)))
    vec1_2 = np.array(anchor2.vector_to(lig2.get_center_of_mass()))
    vec2_2 = -1 * np.array(anchor2.vector_to(np.zeros(3)))

    # Rotate the ligands
    lig1_ar = rot_mol_angle(lig1, vec1_1, vec2_1, idx=idx1, atoms_other=anchor1, bond_length=False)
    lig2_ar = rot_mol_angle(lig2, vec1_2, vec2_2, idx=idx2, atoms_other=anchor2, bond_length=False)
    lig1.from_array(lig1_ar)
    lig2.from_array(lig2_ar)

    # Construct the CdX2 molecule
    CdX2 = Molecule()
    CdX2.add_atom(Atom(atnum=to_atnum(ion)))
    CdX2.merge_mol([lig1, lig2])
    CdX2.properties.name = 'YX2'
    CdX2.properties.path = mol.properties.path
    CdX2.properties.indices = [1, 1 + idx1, 2 + len(lig2) + idx2]

    return CdX2


def get_topology(mol, idx_list, max_dist=5.0):
    """ Return the topology of all atoms *idx_list*, a list of atomic indices. The returned topology
    is based on the number of atoms with a radius *max_dist* from a reference atom. Only atoms with
    the same atomic number as those in *idx_list* will be considered, which in turn should only
    contain atoms of a single element.

    mol <plams.Molecule>: A PLAMS molecule.
    idx_list <list>: A list of atomic indices.
    max_dist <float>: The maximum interatomic distance which is to be considered in topology
        determination.
    return <list>: A list of the topologies in idx_list. Can be either a <str> or <None> if the
        topology is not recognized.
    """
    # Create a dictionary which translates the number of neighbours to a topology
    neighbour_dict = {0: None, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None,
                      7: 'Vertice', 8: 'Edge', 9: 'Face', 10: None, 11: None}
    atnum = mol[idx_list[0]].atnum

    # Create an array with the number of neighbouring atoms in idx_list
    array = mol.as_array(atom_subset=[mol[idx] for idx in idx_list])
    array_other = mol.as_array(atom_subset=[at for at in mol.atoms if at.atnum == atnum])
    dist = cdist(array, array_other)
    dist_count = np.bincount(np.where(dist <= max_dist)[0])

    return [neighbour_dict[i] for i in dist_count]

