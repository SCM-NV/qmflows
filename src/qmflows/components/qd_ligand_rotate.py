__all__ = ['ligand_to_qd']

from scipy.spatial.distance import cdist
import numpy as np

from scm.plams.core.settings import Settings

from .qd_import_export import export_mol
from .qd_functions import (to_array, from_iterable, get_center_of_mass_np, merge_mol,
                           get_atom_index)


def sanitize_dim_2(arg):
    """
    Convert a PLAMS atom or iterable consisting of n PLAMS atoms into a n*3 array.
    In addition, 3 arrays are converted into n*3 arrays.
    """
    if not isinstance(arg, np.ndarray):
        try:
            return np.array(arg.coords)[None, :]
        except AttributeError:
            return to_array(arg)
    else:
        if len(arg.shape) == 1:
            return arg[None, :]
        return arg


def sanitize_dim_3(arg, padding=np.nan):
    """
    Convert an iterable consisting of n PLAMS atoms or a nested iterable consisting of m*(≤n)
    PLAMS atoms into a m*n*3 array, padding the array n is not constant. In addition, n*3 arrays
    are converted into m*n*3 arrays.
    """
    if not isinstance(arg, np.ndarray):
        try:
            return to_array(arg)[None, :, :]
        except AttributeError:
            max_at = max(len(mol) for mol in arg)
            ret = np.empty((len(arg), max_at, 3), order='F')
            ret[:] = padding
            for i, mol in enumerate(arg):
                ret[i, 0:len(mol)] = to_array(mol)
            return ret
    else:
        if len(arg.shape) == 2:
            return arg[None, :, :]
        return arg


def rot_mol_angle(xyz_array, vec1, vec2, idx=0, atoms_other=None, bond_length=False):
    """
    Define a m*3*3 rotation matrix using vec1 and vec2.
    Depending on the dimension of xyz_array, vec1 & vec2 the following operations can be conducted:
        Rotate 1 molecule by 1 rotation matrix     >  1 rotated molecule
        Rotate 1 molecule by m rotation matrices   >  m copies of the molecule, each rotated differently
        Rotate m molecules by 1 rotation matrix    >  m molecules rotated in an identical fashion
        Rotate m molecules by m rotation matrices  >  m molecules, each rotated differently

    Numpy arrays and (nested) iterable consisting of PLAMS atoms can be used interchangeably for
    xyz_array and mol_other; Atomic coordinates will be extracted if necessary and cast into an
    appropriately shaped array.

    xyz_array <np.ndarray>: A m*n*3 array, a n*3 numpy array or a (nested) iterable consisting of
        the to be rotated PLAMS atoms.
    vec1 & vec2 <np.ndarray>: Two n*3 and/or 3 arrays representing one or more vectors.
        vec1 defines the initial vector(s) and vec2 defines the final vector(s) adopted by the to
        be rotated molecule(s).
    atoms_other <None> or <plams.Molecule>: None or a n*3 array, a 3 numpy array, a PLAMS atom or
        iterable consisting of PLAMS atoms. All molecules will be translated to ensure that the
        atom with the index idx adopts the same position as mol_other.
    idx <int>: An atomic index or iterable consisting of multiple atomic indices.
    bond_length <float>: A float or iterable consisting of floats. After translation from
        xyz_array[:, idx] to mol_other, xyz_array will be further translated along vec2 by
        a user-specified bond length.
    return <np.ndarray>: A m*n*3 array or n*3 array representing xyz coordinates of m rotated
        molecules, each consisting of n atoms.
    """
    def get_rotmat(vec1, vec2):
        """ Calculate the rotation matrix for rotating m1 vectors in vec1 to m2 vectors in vec2.
        The following values of m1 & m2 are acceptable:
            m1 = m2
            m1 = 1, m2 > 1
            m1 > 1, m2 = 1 """
        a = vec1 / np.linalg.norm(vec1)
        b = vec2 / np.linalg.norm(vec2, axis=1)[:, None]
        v1, v2, v3 = np.cross(a, b).T
        zero = np.zeros(len(b))
        M = np.array([[zero, -v3, v2],
                      [v3, zero, -v1],
                      [-v2, v1, zero]]).T
        return np.identity(3) + M + ((M@M).T / (1 + b@a.T).T).T

    # Turn all arguments into numpy arrays with appropiate dimensions
    xyz_array = sanitize_dim_3(xyz_array)
    vec1 = sanitize_dim_2(vec1)
    vec2 = sanitize_dim_2(vec2)

    # Rotate and translate all n ligands; readjust bond lengths if bond_length is set
    rotmat = get_rotmat(vec1, vec2)
    xyz_array = xyz_array@rotmat

    # Translate the the molecules in xyz_array
    if atoms_other is not None:
        atoms_other = sanitize_dim_2(atoms_other)
        xyz_array += (atoms_other - xyz_array[:, idx])[:, None, :]
        if bond_length:
            mult = (np.asarray(bond_length) / np.linalg.norm(vec2, axis=1))[:, None]
            xyz_array -= (vec2 * mult)[:, None, :]

    # Return a n*3 or m*n*3 array
    if xyz_array.shape[0] == 1:
        return xyz_array[0]
    return xyz_array


def rot_mol_axis(xyz_array, vec, dist_to_self=True, atoms_other=None, step=(1/16), idx=0):
    """
    Rotates a m*n*3 array representing the cartesian coordinates of m molecules, each containging
    n atoms. All m coordinates are rotated along an axis defined by vec, a m*3 array, yielding
    k = (2 / step) possible conformations.
    Returns all conformations if dist_to_self = True and atoms_other = None, resulting in a
    m*k*n*3 array. Alternatively, returns the conformation of each molecule which maximizes* the
    inter-moleculair distance, yielding a m*n*3 array.

    * More specifically, the conformation which minimizes sum(e^-r(ij)), where i loops over a set
    of n atoms in a single molecule and j over a set of all atoms within mol_other and all
    previously rotated <m molecules.

    xyz_array <np.ndarray>: A m*n*3 array, a n*3 numpy array or a (nested) iterable consisting of
        the to be rotated PLAMS atoms.
    vec <np.ndarray>: An m*3 or 3 numpy array representing m vectors. All m molecules in xyz_array
        will be rotated along an axis defined by their respective vector.
    dist_to_self <bool>: If true, consider the inter-moleculair distance(s) between a molecule in
        xyz_array and all other (previously iterated) molecules in xyz_array when determining the
        conformation that maximizes the inter-moleculair distance.
    atoms_other <None> or <plams.Molecule>: None or a n*3 array, a 3 numpy array or a PLAMS atom
         or iterable consisting of PLAMS atoms. mol_other will be taken into consideration when
         picking the conformation(s) in xyz_array that maximize the inter-moleculair distance(s).
    step <float>: Rotation stepsize in radian, yielding k = (2 / step) possible rotations.
    idx <int>: An atomic index; the molecule will be rotated to ensure it's cartesian coordinates
        remain unchanged during the rotation.
    ret_array <np.ndarray>: A m*n*3 or n*3 numpy array representing the xyz coordinates of
        m molecules with n atoms, all m molecules rotated to maximize the distance among
        themselves and with mol_other.
    """
    def get_rotmat(vec, step=(1/16)):
        """ Calculate the rotation matrix for rotating m vectors along their axes, each vector
        yielding k = (2 / step) possible rotations. """
        v = vec / np.linalg.norm(vec, axis=1)[:, None]
        v1, v2, v3 = v.T
        zero = np.zeros(len(vec))
        W = np.array([[zero, -v3, v2],
                      [v3, zero, -v1],
                      [-v2, v1, zero]]).T

        step_range = np.arange(0.0, 2.0, step)
        a1 = np.sin(step_range)[:, None, None, None]
        a2 = (np.sin(0.5 * step_range)**2)[:, None, None, None]
        return np.identity(3) + a1 * W + a2 * W@W

    # Turn all arguments into numpy arrays with appropiate dimensions
    xyz_array = sanitize_dim_3(xyz_array)
    vec = sanitize_dim_2(vec)
    if atoms_other is None:
        atoms_other = np.array([np.nan, np.nan, np.nan])[None, :]
    else:
        atoms_other = sanitize_dim_2(atoms_other)

    # Create all k possible rotations of all m ligands
    rotmat = get_rotmat(vec, step=step)
    xyz_array = np.swapaxes(xyz_array@rotmat, 0, 1)
    xyz_array += (xyz_array[:, 0, idx][:, None, :] - xyz_array[:, :, idx])[:, :, None, :]

    # Returns the conformation of each molecule that maximizes the inter-moleculair distance
    # Or return all conformations if dist_to_self = False and atoms_other = None
    if dist_to_self or atoms_other is not None:
        a, b, c, d = xyz_array.shape
        ret_array = np.empty((a, c, d), order='F')
        for i, xyz in enumerate(xyz_array):
            dist_array = cdist(xyz.reshape(b*c, d), atoms_other).reshape(b, c, len(atoms_other))
            idx_min = np.nansum(np.exp(-dist_array), axis=(1, 2)).argmin()
            if dist_to_self:
                atoms_other = np.concatenate((atoms_other, xyz[idx_min]))
            ret_array[i] = xyz[idx_min]

        # Return a n*3 or m*n*3 array
        if ret_array.shape[0] == 1:
            return ret_array[0]
        return ret_array

    # Return a n*3 or m*n*3 array
    if xyz_array.shape[0] == 1:
        return xyz_array[0]
    return xyz_array


def array_to_qd(mol, xyz_array, mol_other=False):
    """
    Takes a template molecule with n atoms and create m copies of this molecule, updating its
        cartesian coordinates based on a m*n*3 array.
    Has the option to directly combine all copied molecules with a second molecule (mol_other)
        instead of returning the rotated molecules.

    mol <plams.Molecule>: A template PLAMS molecule consisting of n atoms.
    xyz_array <np.ndarray>: A m*n*3 or n*2 numpy array representing the cartesian coordinates of
        m molecules each with n atoms.
    mol_other <plams.Molecule> or False: Add all atoms and bonds from the to be returned molecule
        list to this molecule, instead of returning the molecule list.
    return <list> or None: Returns a list of rotated molecules or concatenates aforementioned
        list and add its atoms and bonds to mol_other, returning nothing.
    """
    mol_list = []
    xyz_array = sanitize_dim_3(xyz_array)
    for i, xyz in enumerate(xyz_array, 2):
        mol_cp = mol.copy()
        from_iterable(mol_cp, xyz, obj='array')
        for at in mol_cp:
            at.properties.pdb_info.ResidueNumber = i
        mol_list.append(mol_cp)

    # return a list of molecules or merge the list with an existing molecule
    if not mol_other:
        return mol_list
    mol_other.merge_mol(mol_list)


def ligand_to_qd(core, ligand, qd_folder):
    """
    Function that handles quantum dot (qd, i.e. core + all ligands) operations.
    Combine the core and ligands and assign properties to the quantom dot.

    core <plams.Molecule>: The core molecule.
    ligand <plams.Molecule>: The ligand molecule.
    qd_folder <str>: The quantum dot export folder.

    return <plams.Molecule>: The quantum dot (core + n*ligands).
    """
    # Define vectors and indices used for rotation and translation the ligands
    vec1 = sanitize_dim_2(ligand.properties.dummies) - sanitize_dim_2(ligand.get_center_of_mass_np())
    vec2 = sanitize_dim_2(core.get_center_of_mass_np()) - sanitize_dim_2(core.properties.dummies)
    idx = ligand.properties.dummies.get_atom_index() - 1
    ligand.properties.dummies.properties.anchor = True

    # Attach the rotated ligands to the core, returning the resulting strucutre (PLAMS Molecule).
    lig_array = rot_mol_angle(ligand, vec1, vec2, atoms_other=core.properties.dummies, idx=idx)
    lig_array = rot_mol_axis(lig_array, vec2, atoms_other=core.properties.dummies, idx=idx)
    qd = core.copy()
    array_to_qd(ligand, lig_array, mol_other=qd)

    # indices of all the atoms in the core and the ligand heteroatom anchor.
    qd_indices = [i for i, at in enumerate(qd, 1) if at.properties.pdb_info.ResidueName == 'COR' or
                  at.properties.anchor]

    qd_name = core.properties.name + '__' + str(qd[-1].properties.pdb_info.ResidueNumber - 1)
    qd_name += '_' + ligand.properties.name + '@' + ligand.properties.group
    qd.properties = Settings()
    qd.properties.indices = qd_indices
    qd.properties.name = qd_name
    qd.properties.path = qd_folder
    export_mol(qd, message='core + ligands:\t\t')

    return qd
