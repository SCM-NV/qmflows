"""Functions for reading CP2K orbitals."""

from __future__ import annotations

import os
import re
import fnmatch
import types
from itertools import islice
from pathlib import Path
from collections import defaultdict

import numpy as np

from ..common import CP2KInfoMO, MO_metadata

__all__ = ["read_cp2k_coefficients"]


def read_cp2k_coefficients(
    path_mos: str | os.PathLike[str],
    plams_dir: None | str | os.PathLike[str] = None,
) -> CP2KInfoMO | tuple[CP2KInfoMO, CP2KInfoMO]:
    """Read the MO's from the CP2K output.

    First it reads the number of ``Orbitals`` and ``Orbital`` functions from the
    cp2k output and then read the molecular orbitals.

    Returns
    -------
        NamedTuple containing the Eigenvalues and the Coefficients
    """
    plams_dir = Path(plams_dir) if plams_dir is not None else Path(os.getcwd())
    file_out = fnmatch.filter(os.listdir(plams_dir), '*out')[0]
    path_out = plams_dir / file_out

    orbitals_info = read_cp2k_number_of_orbitals(path_out)
    return read_log_file(path_mos, orbitals_info)


def read_log_file(
    path: str | os.PathLike[str],
    orbitals_info: MO_metadata,
) -> CP2KInfoMO | tuple[CP2KInfoMO, CP2KInfoMO]:
    """
    Read the orbitals from the Log file.

    Notes
    -----
    IF IT IS A UNRESTRICTED CALCULATION, THERE ARE TWO SEPARATED SET OF MO
    FOR THE ALPHA AND BETA

    Parameters
    ----------
    path : path-like object
        Path to the file containing the MO coefficients
    orbitals_info : MO_metadata
        A named tuple with the number of orbital functions and the number of MOs.
        Note that the latter is an upper bound to the number of MOs printed in ``path``.

    Returns
    -------
    CP2KInfoMO or tuple[CP2KInfoMO, CP2KInfoMO]
        Molecular orbitals and orbital energies.
        Returns a tuple of two ``CP2KInfoMO`` objects if its an unrestricted calculation.

    """
    n_orb_list = orbitals_info.nOrbitals
    n_orb_funcs = orbitals_info.nOrbFuns

    unrestricted = orbitals_info.nspinstates == 2
    start = _find_mo_start(path, unrestricted=unrestricted)
    if not unrestricted:
        return read_coefficients(path, n_orb_list[0], n_orb_funcs, start)
    else:
        return (
            read_coefficients(path, n_orb_list[0], n_orb_funcs, start[0]),
            read_coefficients(path, n_orb_list[1], n_orb_funcs, start[1]),
        )


def _find_mo_start(path: str | os.PathLike[str], unrestricted: bool) -> int | tuple[int, int]:
    """Find the start of the final-most MO range in ``path``.

    Returns
    -------
    int or tuple[int, int]
        The start of the MO range (excluding the header).
        Returns a two-tuple for unrestricted calculations (for the alpha and beta component)
        and a single integer otherwise.

    """
    with open(path, "r") as f:
        # Find all headers, but skip restarts
        enumerator = enumerate((h.strip("MO| \n") for h in f), start=1)
        headers = [(i, h) for i, h in enumerator if "EIGENVALUES" in h]
        headers = [(i, h) for i, h in headers if "AFTER SCF STEP" not in h]
    err = "Failed to identify the start of the {prefix}MO range in {path!r}"

    if not unrestricted:
        if not len(headers):
            raise ValueError(err.format(prefix="", path=os.fspath(path)))
        return headers[-1][0]

    alphas = [i for i, h in headers if "ALPHA" in h]
    betas = [i for i, h in headers if "BETA" in h]
    if not len(alphas):
        raise ValueError(err.format(prefix="alpha ", path=os.fspath(path)))
    if not len(betas):
        raise ValueError(err.format(prefix="beta ", path=os.fspath(path)))
    return alphas[-1], betas[-1]


def read_coefficients(
    path: str | os.PathLike[str],
    n_orb: int,
    n_orb_funcs: int,
    start: int = 0,
) -> CP2KInfoMO:
    """Read the coefficients from the plain text output.

    MO coefficients are stored in Column-major order.
    CP2K molecular orbitalsoutput (CP2K <8.2) looks like:

    .. code-block::

        MO EIGENVALUES, MO OCCUPATION NUMBERS, AND SPHERICAL MO EIGENVECTORS

                                  5                      6
                              -0.2590267204166110    -0.1785544120250688

                               2.0000000000000000     2.0000000000000000

        1      1  C  2s        0.0021482361354044     0.0000000235522485
        2      1  C  3s       -0.0007100065367389     0.0000000102096730
        3      1  C  3py      -0.1899052318987045     0.0000000059435027
        4      1  C  3pz       0.0000000178537720    -0.5500605729231620
        5      1  C  3px       0.3686765614894165     0.0000000228716009
        6      1  C  4py       0.0014072130025389     0.0000000019199413
        7      1  C  4pz      -0.0000000014121887     0.0293850516018881
        8      1  C  4px      -0.0028383911872079     0.0000000042372601
        9      1  C  4d-2      0.0311183981707317     0.0000000014108937
        10     1  C  4d-1      0.0000000095952723     0.0253837978837068
        11     1  C  4d0       0.0005419630026689     0.0000000391888080
        12     1  C  4d+1     -0.0000000210955114     0.0147105486663415
        13     1  C  4d+2      0.0534202997324328     0.0000000021056315

    """
    energies = np.empty(n_orb, dtype=np.float64)
    coefs = np.empty((n_orb_funcs, n_orb), dtype=np.float64)
    orb_index = np.empty(n_orb, dtype=np.int64)
    occupation = np.empty(n_orb, dtype=np.float64)

    j0 = 0
    j1 = 0
    with open(path, "r") as f:
        iterator = filter(None, (i.strip("MO| \n") for i in islice(f, start, None)))
        for h in iterator:
            # Each MO pair is preceded by their indices,
            # a lack thereof meaning the end of the file has been reached
            headers = h.split()
            if not all(i.isnumeric() for i in headers):
                break

            j1 += len(headers)
            orb_index[j0:j1] = headers
            energies[j0:j1] = next(iterator).split()
            occupation[j0:j1] = next(iterator).split()
            coefs[..., j0:j1] = [i.split()[4:] for i in islice(iterator, 0, n_orb_funcs)]
            j0 += len(headers)

    # `n_orb` is an upper bound to the actual printed number of orbitals,
    # so some trimming might be required
    return CP2KInfoMO(energies[:j0], coefs[:, :j0], orb_index[:j0], occupation[:j0])


_ORBITAL_PATTERN = re.compile((
    r"(?P<key>Number of occupied orbitals|Number of molecular orbitals|Number of orbital functions):"  # noqa: E501
    r"\s+"
    r"(?P<value>[0-9]+)"
))

# Map CP2K field names to qmflows `MO_metadata` field names
_ORBITAL_NAMES = types.MappingProxyType({
    "Number of occupied orbitals": "nOccupied",
    "Number of molecular orbitals": "nOrbitals",
    "Number of orbital functions": "nOrbFuns",
})

_REQUIRED_KEYS = frozenset(_ORBITAL_NAMES.values())


def read_cp2k_number_of_orbitals(file_name: str | os.PathLike[str]) -> MO_metadata:
    """Extract the number of (occupied) orbitals and basis functions."""
    with open(file_name, "r") as f:
        kwargs: "dict[str, list[int]]" = defaultdict(list)
        for line in f:
            match = _ORBITAL_PATTERN.search(line)
            if match is None:
                continue

            orb_name = _ORBITAL_NAMES[match["key"]]
            orb_count = int(match["value"])
            kwargs[orb_name].append(orb_count)

            # Check if all keys are present and return
            if _REQUIRED_KEYS == kwargs.keys():
                n_orb_funcs = kwargs.pop("nOrbFuns")
                assert len(n_orb_funcs) == 1
                return MO_metadata(**kwargs, nOrbFuns=n_orb_funcs[0])
        else:
            missing_fields = sorted(k for k, v in _ORBITAL_NAMES.items() if v not in kwargs)
            missing = "".join(f"\n  - {i!r}" for i in missing_fields)
            raise ValueError(
                f"Failed to parse one or more fields in {os.fspath(file_name)!r}:{missing}"
            )
