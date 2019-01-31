from os.path import join
import pdb

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist

from scm.plams.core.basemol import Molecule

from qmflows.components.qd_functions import to_array


def get_radial_distr(array1, array2, dr=0.01, r_max=10.0):
    """ Calculate the radial distribution function between *array1* and *array2*.
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
    dens /= (4 * np.pi * r**2 * dr)
    dens /= dens_mean

    return dens


path = '/Users/basvanbeek/Documents/CdSe/Week_5/QD'
mol = 'Core_Cd360Se309__102_Ligand_CCCCCCCCCCCCCCCO@O16.pdb'
mol = Molecule(join(path, mol))
mol.properties.name = 'QD'

# Create xyz arrays with Cd, Se, O
Cd_array = to_array([at for at in mol if at.symbol == 'Cd'])[None, :, :]
Se_array = to_array([at for at in mol if at.symbol == 'Se'])[None, :, :]
O_array = to_array([at for at in mol if at.symbol == 'O'])[None, :, :]

# Plot it
dr, r_max = 0.02, 10.0
r = np.arange(dr, r_max + dr, dr)
dens = np.array([get_radial_distr(i, j, dr=dr) for i, j in zip(Cd_array, O_array)]).mean(axis=0)
plt.plot(r, dens, 'r')
