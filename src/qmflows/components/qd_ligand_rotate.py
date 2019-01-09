__all__ = ['rotation_check']

import numpy as np
from scipy.spatial.distance import cdist

from .qd_functions import to_array, update_coords


def rotate(vec, array, step=0.5, anchor=0):
    """
    Rotate an array along an axis defined by vec, yielding all possible rotated arrays given a
        rotation between 0*pi and 2*pi with a user-specified rotation stepsize (radian).
    vec <np.ndarray>: A vector with 3 elements.
    array <np.ndarray>: A n*3 xyz array.
    step <float>: The rotation stepsize in radian, yielding m = int(2 * step**-1) angles.
    anchor <int>: The index of the anchor; the array will be translated to ensure that the
        coordinates of array[anchor] remain unchanged upon rotation.
    return <np.ndarray>: A m*n*3 xyz array, each 2d array along the axis m representing an array
        rotated at a different angle.
    """
    vec /= np.linalg.norm(vec)

    W = np.array([[0, -vec[2], vec[1]],
                 [vec[2], 0, -vec[0]],
                 [-vec[1], vec[0], 0]])

    r = np.arange(0.0, 2.0, step)
    a1 = np.sin(r)[:, None, None]
    a2 = (np.sin(0.5 * r)**2)[:, None, None]
    rotmat = np.identity(3) + a1 * W + a2 * np.dot(W, W)

    ret = np.matmul(rotmat, array.T)
    ret += (array[anchor] - ret.T[anchor].T)[:, :, None]
    ret = np.swapaxes(ret, 1, 2)
    a, b, c = ret.shape
    return ret.reshape(a*b, c)


def rotation_check(mol, step=0.5, anchor=0, radius=10.0):
    """
    mol <plams.Molecule>: A PLAMS molecule with each atom marked by a residue number.
    step <float>: The rotation stepsize in radian, yielding m = int(2 * step**-1) angles.
    anchor <int>: The index of the anchor within a ligand residue; ligand arrays will be translated
        to ensure that the coordinates of array[anchor] remain unchanged upon rotation.
    radius <float>: The radius
    """
    xyz = mol.to_array()

    def count_res(iterable, res=1):
        ret = 0
        for at in iterable:
            if at.properties.pdb_info.ResidueNumber == res:
                ret += 1
            else:
                return ret

    # cor/lig: the number of atoms in the core/ligands;
    # a: the number of ligand residues;
    # r: the number of rotations per ligand
    cor = count_res(mol.atoms, res=1)
    lig = count_res(mol.atoms[cor:], res=2)
    a = len({at.properties.pdb_info.ResidueNumber for at in mol}) - 1
    r = int(2 * step**-1)

    dist = cdist(xyz[cor:], xyz[cor:]).reshape(a, lig, len(xyz[cor:]))

    for i, dist_slice in enumerate(dist):
        # Set all intra-ligand distances to >100.0 A, ignoring them
        idx1 = cor + i*lig
        idx2 = cor + (i+1)*lig
        dist_slice[:, idx1:idx2] = radius + 100.0

        # Create a vector from the reference ligand anchor to the center of the reference ligand
        # The reference ligand is rotated along this vector in user-defined steps
        vec = xyz[idx1 + anchor] - xyz[idx2]
        xyz2 = rotate(vec, xyz[idx1:idx2], step=step, anchor=anchor)

        # Create an array with all ligands within a radius of 10.0 A from the reference ligand
        idx_other = np.unique(np.where(dist_slice < radius)[1] // lig)
        xyz_other = np.concatenate([xyz[j*lig:(j+1)*lig] for j in idx_other] + [xyz[0:cor]])

        dist2 = cdist(xyz2, xyz_other).reshape(r, lig, len(xyz_other))

        # Find the conformation with the largest distance between atoms in
        # the reference ligand and all other ligands. Weighting sum: sum(e^-r(ij))
        idx_min = np.argmin(np.sum(np.exp(-dist2), axis=(1, 2)))
        xyz[idx1:idx2] = xyz2[idx_min*lig:(idx_min+1)*lig]

    mol.update_coords(xyz, obj='array')
