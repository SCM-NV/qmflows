__all__ = ['get_bde']

import numpy as np
import pandas as pd

from scm.plams import (Molecule, Atom)
from scm.plams.tools.units import Units

from .qd_functions import to_atnum
from .qd_ligand_rotate import rot_mol_angle
from .qd_ams import ams_job_mopac_opt, ams_job_mopac_sp, ams_job_uff_opt
from .qd_functions import merge_mol
from .qd_dissociate import dissociate_ligand


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
    CdX2.properties.path = mol.properties.path
    CdX2.properties.indices = [1, 1 + idx1, 2 + len(lig2) + idx2]
    CdX2 = ams_job_mopac_opt(CdX2)

    return CdX2


def get_bde(mol, get_ddG=True, unit='kcal/mol'):
    """ Calculate the bond dissociation energy: dE = dE(mopac) + (dG(uff) - dE(uff))
    """
    # Perform MOPAC single points
    core = dissociate_ligand(mol)
    lig = get_cdx2(mol)
    lig.properties.name = 'CdX2'
    tot = ams_job_mopac_sp(mol)

    # Extract MOPAC total energies
    E_lig_mopac = lig.properties.energy.E
    E_core_mopac = np.array([ams_job_mopac_sp(mol).properties.energy.E for mol in core])
    E_tot_mopac = tot.properties.energy.E

    # Calculate dE
    dE_mopac = (E_lig_mopac + E_core_mopac) - E_tot_mopac
    dE_mopac *= Units.conversion_ratio('Hartree', 'kcal/mol')

    res_list1, idx_list, res_list2 = zip(*[mol.properties.mark for mol in core])

    if get_ddG:
        # Perform UFF constrained geometry optimizations + frequency analyses
        lig = ams_job_uff_opt(lig, get_freq=True, fix_angle=False)
        core = [ams_job_uff_opt(mol, get_freq=True, fix_angle=False) for mol in core]
        tot = ams_job_uff_opt(tot, get_freq=True, fix_angle=False)

        # Extract UFF total Gibbs free energies
        G_lig_uff = lig.properties.energy.G
        G_core_uff = np.array([mol.properties.energy.G for mol in core])
        G_tot_uff = tot.properties.energy.G

        # Extract UFF total energies
        E_lig_uff = lig.properties.energy.E
        E_core_uff = np.array([mol.properties.energy.E for mol in core])
        E_tot_uff = tot.properties.energy.E

        # Calculate ddG
        dG_uff = (G_lig_uff + G_core_uff) - G_tot_uff
        dE_uff = (E_lig_uff + E_core_uff) - E_tot_uff
        ddG_uff = dG_uff - dE_uff

        return pd.DataFrame({'First ligand': res_list1,
                             'Cd index': idx_list,
                             'Second ligand': res_list2,
                             'dE': dE_mopac,
                             'ddG': ddG_uff,
                             'dG': dE_mopac+ddG_uff})

    return pd.DataFrame({'First ligand': res_list1,
                         'Cd index': idx_list,
                         'Second ligand': res_list2,
                         'dE': dE_mopac})


def get_y_axis():
    pass
