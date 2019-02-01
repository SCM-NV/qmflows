__all__ = []

from os.path import (dirname, join)
import time
import pdb

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from scm.plams.core.basemol import Molecule
from scm.plams.tools.units import Units

from qmflows.components.qd_functions import to_array, from_iterable, get_time


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


def get_radial_distr(array1, array2, dr=0.05, r_max=12.0):
    """ Calculate the radial distribution function between *array1* and *array2*: g(r_ij).

    array1 <np.ndarray>: A n*3 array representing the cartesian coordinates of the reference atoms.
    array2 <np.ndarray>: A m*3 array representing the cartesian coordinates of (non-reference) atoms.
    dr <float>: The integration step-size in Angstrom.
    r_max <float>: The maximum to be returned interatomic distance.
    """
    idx_max = 1 + int(r_max / dr)
    dist = cdist(array1, array2)
    dist_int = np.array(dist / dr, dtype=int).flatten()

    # Calculate the average particle density N / V
    # The diameter of the spherical volume (V) is defined by the largest inter-particle distance: max(r_ij)
    dens_mean = len(array2) / ((4/3) * np.pi * (0.5 * dist.max())**3)

    # Count the number of occurances of each (rounded) distance; the first element (0 A) is skipped
    dens = np.bincount(dist_int)[1:idx_max]

    # Correct for the number of reference atoms
    dens = dens / len(array1)

    # Convert the particle count into a partical density
    r = np.arange(dr, r_max + dr, dr)
    try:
        dens /= (4 * np.pi * r**2 * dr)
    except ValueError:
        # Plan b: Pad the array with zeros if r_max is larger than dist.max()
        zeros = np.zeros(len(r))
        zeros[0:len(dens)] = dens
        dens = zeros / (4 * np.pi * r**2 * dr)

    # Normalize and return the particle density
    return dens / dens_mean


def get_all_radial(mol_list, dr=0.05, r_max=12.0, atoms=('Cd', 'Se', 'O')):
    """ Return the radial distribution functions for all possible atom-pairs in *elements*
    as dataframe. Accepts both single molecules and list of molecules as input, the later allowing
    for conformational and/or configurational averaging.

    mol_list <list>[<plams.Molecule>]: A PLAMS molecule or list of PLAMS molecules.
    atoms <tuple>[<str>]: A tuple of strings representing atomic symbols; RDF's will be calculated
        for all unique atomp pairs.
    return <pd.DataFrame>: A Pandas dataframe of radial distribution functions.
    """
    # Make sure we're dealing with a list of molecules
    if isinstance(mol_list, Molecule):
        mol_list = [mol_list]

    # Create a dictionary consisting of {atomic symbol: [atomic indices]}
    idx_dict = {}
    for i, at in enumerate(mol_list[0]):
        try:
            idx_dict[at.symbol].append(i)
        except KeyError:
            idx_dict[at.symbol] = [i]

    # Create a dataframe of RDF's, summed over all conformations in mol_list
    df = pd.DataFrame(index=np.arange(dr, r_max + dr, dr))
    for mol in mol_list:
        xyz_array = to_array(mol)
        for i, at1 in enumerate(atoms):
            for at2 in atoms[i:]:
                try:
                    df[at1 + '_' + at2] += get_radial_distr(xyz_array[idx_dict[at1]],
                                                            xyz_array[idx_dict[at2]],
                                                            dr=dr, r_max=r_max)
                except KeyError:
                    df[at1 + '_' + at2] = get_radial_distr(xyz_array[idx_dict[at1]],
                                                           xyz_array[idx_dict[at2]],
                                                           dr=dr, r_max=r_max)

    # Average the RDF's over all conformations in mol_list
    df /= len(mol_list)
    return df.rename_axis('r(ij) / A')


path = dirname(__file__)
name = 'Cd360Se309X102.pdb'
mol = Molecule(join(path, name))

# Create copies of mol and randomly displace all atoms between 0.00 and 0.25 Angstrom
# Skip this step if the variable mol_list is alreayd assigned and is of length n
n = 100
print('')
try:
    mol_list == True
    if len(mol_list) != n:
        name_error == True
except NameError:
    print(get_time() + 'generating mol_list')
    mol_list = [mol.copy() for i in range(n)]
    for mol in mol_list:
        ar = to_array(mol)
        ar += np.random.rand(len(ar), 3) / 4
        from_iterable(mol, ar, obj='array')
        del ar
    print(get_time() + 'mol_list has been generated')

# Run the actual script and plot the results
dr, r_max = 0.05, 12.0
start = time.time()
df = get_all_radial(mol_list, dr=dr, r_max=r_max)
print(get_time() + 'run time:', '%.2f' % (time.time() - start), 'sec')
df.plot()
