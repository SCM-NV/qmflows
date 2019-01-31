from os.path import join
import pdb

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from scm.plams.core.basemol import Molecule
from scm.plams.tools.units import Units

from qmflows.components.qd_functions import to_array


def get_radial_distr(array1, array2, dr=0.01, r_max=10.0):
    """ Calculate the radial distribution function between *array1* and *array2*: g(r_ij).

    array1 <np.ndarray>: A n*3 array representing the cartesian coordinates of the reference atoms.
    array2 <np.ndarray>: A m*3 array representing the cartesian coordinates.
    dr <float>: The integration step-size in Angstrom.
    r_max <float>: The maximum to be returned interatomic distance.
    """
    dist = cdist(array1, array2)
    dist_int = np.array(dist / dr, dtype=np.int64).flatten()

    r = np.arange(dr, r_max + dr, dr)
    dens_mean = len(array2) / ((4/3) * np.pi * (0.5 * dist.max())**3)

    dens = np.bincount(dist_int)[1:1+int(r_max/dr)]
    dens = dens / len(array1)
    try:
        dens /= (4 * np.pi * r**2 * dr)
    except ValueError:
        zeros = np.zeros(len(r))
        zeros[0:len(dens)] = dens
        dens = zeros / (4 * np.pi * r**2 * dr)
    dens /= dens_mean

    return dens


def get_rel_error(g_QM, g_MM, T=298.15, unit='kcal/mol'):
    """ Return the relative error defined as dF_rij = dXi = Xi_QM - Xi_MM.

    g_QM & G_MM <np.ndarray>: A m*n numpy arrays of m radial distribution functions of QM & MM
        calculations, respectively.
    T <float>: The temperature in Kelvin.
    unit <str>: The unit of the to be returned energy
    return <np.array>: The relative error dXi.
    """
    RT = 0.00198720 * T * Units.conversion_ratio('kcal/mol', unit)
    return -RT * np.log1p((g_MM - g_QM) / g_QM)


def get_aux_error(g_QM, g_MM):
    """ Return the auxiliary error defined as dEps = Eps_QM - Eps_MM.

    g_QM & G_MM <np.ndarray>: A m*n numpy arrays of m radial distribution functions of QM & MM
        calculations, respectively.
    return <float>: The auxiliary error dEps.
    """
    return np.linalg.norm(g_QM - g_MM, axis=0).sum()


def get_increment(phi_old, a_old, a, y=2.0):
    """ Return the incremental factor Phi for the (k)th block.
    phi_old <float>: Phi_Omega produced in the (k - 1)th.
    a_old <np.ndarray>: The accepatance rates from all Omega iterations in the (k - 1)th block.
    a <float>: The target accepatance rate.
    y <float>: Regulates the correction of Phi_old, must be larger than or equal to 1.0.
    return <float>: The new incremental factor Phi.
    """
    return phi_old * y**np.sign(a - a_old.mean())


def get_all_radial(mol, dr=0.05, r_max=12.0, elements=('Cd', 'Se', 'O')):
    """ Return all radial distribution functions for all possible atom-pairs in *elements* as dataframe.
    mol <plams.Molecule>: A PLAMS molecule.
    elements <tuple>[<str>]: A tuple of strings representing atomic symbols.
    return <pd.DataFrame>: A Pandas dataframe of radial distribution functions (numpy arrays).
    """
    xyz_dict = {atom: to_array([at for at in mol if at.symbol == atom]) for atom in elements}
    df = pd.DataFrame(index=np.arange(dr, r_max + dr, dr))
    for i in xyz_dict:
        for j in xyz_dict:
            if j + '_' + i not in df.keys():
                df[i + '_' + j] = get_radial_distr(xyz_dict[i], xyz_dict[j], dr=dr, r_max=r_max)

    return df.rename_axis('r(ij) / A')


path = '/Users/basvanbeek/Documents/CdSe/Week_5/QD'
mol = 'Core_Cd360Se309__102_Ligand_CCCCCCCCCCCCCCCO@O16.pdb'
mol = Molecule(join(path, mol))

# Plot it
dr, r_max = 0.05, 12.0
df = get_all_radial(mol, dr=dr, r_max=r_max)
phi_init = 1.0
w = 100
df.plot()
