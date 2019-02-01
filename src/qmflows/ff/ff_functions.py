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


def get_all_radial(xyz_array, idx_dict, dr=0.05, r_max=12.0, atoms=('Cd', 'Se', 'O')):
    """ Return the radial distribution functions for all possible atom-pairs in *elements*
    as dataframe. Accepts both single molecules and list of molecules as input, the later allowing
    for conformational and/or configurational averaging.

    xyz_array <np.ndarray>: A m*n*3 or n*3 numpy array of cartesian coordinates.
    atoms <tuple>[<str>]: A tuple of strings representing atomic symbols; RDF's will be calculated
        for all unique atomp pairs.
    return <pd.DataFrame>: A Pandas dataframe of radial distribution functions.
    """
    # Make sure we're dealing with a list of molecules
    if len(xyz_array.shape) == 2:
        xyz_array = xyz_array[None, :, :]

    # Create a dataframe of RDF's, summed over all conformations in mol_list
    df = pd.DataFrame(index=np.arange(dr, r_max + dr, dr))
    for xyz in xyz_array:
        for i, at1 in enumerate(atoms):
            for at2 in atoms[i:]:
                try:
                    df[at1 + '_' + at2] += get_radial_distr(xyz[idx_dict[at1]],
                                                            xyz[idx_dict[at2]],
                                                            dr=dr, r_max=r_max)
                except KeyError:
                    df[at1 + '_' + at2] = get_radial_distr(xyz[idx_dict[at1]],
                                                           xyz[idx_dict[at2]],
                                                           dr=dr, r_max=r_max)

    # Average the RDF's over all conformations in mol_list
    df /= xyz_array.shape[0]
    return df.rename_axis('r(ij) / A')


def read_multi_xyz(file, ret_idx_dict=True):
    """ Return a m*n*3 array of cartesian coordinates extracted from a multi-xyz file.
    file <str>: The path + filename to a multi-xyz file.
    ret_idx_dict <bool>: Return a dictionary consisting of {atomic symbols: [atomic indices]} as
        derived from the first molecule in *file*.
    return <np.ndarray> (and <dict>): A m*n*3 numpy array of cartesian coordinates and, optionally,
        a dictionary consisting of {atomic symbols: [atomic indices]}.
    """
    # Read the multi-xyz file
    with open(join(path, name)) as file:
        file = file.read().splitlines()
    mol_size = int(file[0])

    # Create a list of atomic symbols and xyz coordinates
    xyz = []
    at_symbol = []
    for item in file:
        item = item.split()
        if len(item) == 4:
            at_symbol.append(item[0])
            xyz.append(item[1:])

    # Captilize the first letter of each atomic symbol
    for i, at in enumerate(at_symbol[0:mol_size]):
        at_symbol[i] = at.capitalize()

    # Turn the xyz list into a m*n*3 numpy array
    shape = int(len(xyz) / mol_size), mol_size, 3
    xyz = np.array(xyz, dtype=float, order='F').reshape(shape)

    # Turn the atomic symbols list into a dictionary: {atomic symbols: [atomic indices]}
    if ret_idx_dict:
        idx_dict = {}
        for i, at in enumerate(at_symbol[0:mol_size]):
            try:
                idx_dict[at].append(i)
            except KeyError:
                idx_dict[at] = [i]
        return xyz, idx_dict
    return xyz


# Define variables
dr = 0.05
r_max = 12.0
path = dirname(__file__)
name = 'Cd68Se55_26COO_MD_trajec.xyz'

# Run the actual script and plot the results
print('')
start = time.time()
xyz_array, idx_dict = read_multi_xyz(join(path, name))
df = get_all_radial(xyz_array, idx_dict, dr=dr, r_max=r_max)
print(get_time() + 'run time:', '%.2f' % (time.time() - start), 'sec')
df.plot()
