from os.path import join
import pdb

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist

from scm.plams.core.basemol import Molecule

from qmflows.components.qd_functions import to_array


path = '/Users/bvanbeek/Documents/CdSe/Week_5/QD'
mol = 'Core_Cd68Se55__26_Ligand_OC[CCCCCF]CCCCCCCCCC@O1.pdb'
mol = Molecule(join(path, mol))
mol.properties.name = 'QD'

# Create xyz arrays with Cd, Se, O
Cd_array = to_array([at for at in mol if at.symbol == 'Cd'])
Se_array = to_array([at for at in mol if at.symbol == 'Se'])
O_array = to_array([at for at in mol if at.symbol == 'O'])


def get_radial_distr(array1, array2, dr=0.01, r_max=10.0):
    """ Calculate the radial distribution function between *array1* and *array2*.
    dr <float>: The integration step-size in Angstrom.
    r_max <float>: The maximum to be returned interatomic distance.
    """
    dist = cdist(array1, array2)
    dist_int = np.array(dist / dr, dtype=np.int64).flatten()

    r = np.arange(dr, r_max + dr, dr)
    dens_mean = (len(array2) / (4 * np.pi * (dist.max() / 2)**3))
    dens = np.bincount(dist_int)[1:1+int(r_max/dr)]
    dens = dens / len(array1)
    dens /= (4 * np.pi * r**2 * dr)
    dens /= dens_mean

    return dens

# Plot that shit
dr, r_max = 0.01, 10.0
r = np.arange(dr, r_max + dr, dr)
dens = get_radial_distr(Cd_array, Cd_array, dr=dr, r_max=r_max)
plt.plot(r, dens, 'r')

