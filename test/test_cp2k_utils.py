"""Test PSF interface."""
from io import BytesIO, StringIO
from itertools import repeat

import numpy as np
import pandas as pd
from assertionlib import assertion

from qmflows import Settings
from qmflows.cp2k_utils import (CP2K_KEYS_ALIAS, LengthError, _construct_df,
                                _get_key_path, _parse_unit, _validate_unit,
                                map_psf_atoms, set_prm, prm_to_df)
from qmflows.test_utils import PATH_MOLECULES

PSF_STR = """
PSF EXT

    10 !NATOM
     1 MOL1     1        LIG      C        C331   -0.272182       12.010600        0
     2 MOL1     1        LIG      C        C321   -0.282182       12.010600        0
     3 MOL1     1        LIG      C        C2O3    0.134065       12.010600        0
     4 MOL1     1        LIG      O        O2D2   -0.210848       15.999400        0
     5 MOL1     1        LIG      O        O2D2   -0.210848       15.999400        0
     6 MOL1     1        LIG      H        HGA2    0.087818        1.007980        0
     7 MOL1     1        LIG      H        HGA2    0.087818        1.007980        0
     8 MOL1     1        LIG      H        HGA3    0.087818        1.007980        0
     9 MOL1     1        LIG      H        HGA3    0.087818        1.007980        0
    10 MOL1     1        LIG      H        HGA3    0.087818        1.007980        0
"""


def test_CP2K_KEYS_ALIAS() -> None:
    """Tests for :data:`CP2K_KEYS_ALIAS`."""
    key_exclude = {'lennard_jones'}

    assertion.isinstance(CP2K_KEYS_ALIAS, dict)

    for k, v in CP2K_KEYS_ALIAS.items():
        if k.endswith('14'):
            k = k[:-2]
        if k in key_exclude:
            continue

        assertion.eq(k, v[-1])
        assertion.isinstance(v, tuple)
        assertion.eq(v[:2], ('specific', 'cp2k'))


def test_map_psf_atoms() -> None:
    """Tests for :func:`map_psf_atoms`."""
    ref = {'C331': 'C', 'C321': 'C', 'C2O3': 'C', 'O2D2': 'O', 'HGA2': 'H', 'HGA3': 'H'}

    path_like = PATH_MOLECULES / 'mol.psf'
    file_like = StringIO(PSF_STR)
    assertion.eq(map_psf_atoms(path_like), ref)
    assertion.eq(map_psf_atoms(file_like), ref)

    # Incorrect input types
    assertion.assert_(map_psf_atoms, None, exception=TypeError)
    assertion.assert_(map_psf_atoms, {1: 1}, exception=TypeError)
    assertion.assert_(map_psf_atoms, 'bob.psf', exception=FileNotFoundError)

    # Incorrect mode (bytes)
    file_bytes = BytesIO(PSF_STR.encode())
    assertion.assert_(map_psf_atoms, file_bytes, mode='rb', exception=TypeError)

    # The !NATOM will now be missing
    file_wrong1 = StringIO(PSF_STR[24:])
    assertion.assert_(map_psf_atoms, file_wrong1, exception=ValueError)

    # Rows [5:] will now be missing
    file_wrong2 = StringIO('\n'.join(i[:25] for i in PSF_STR.splitlines()))
    assertion.assert_(map_psf_atoms, file_wrong2, exception=ValueError)


def test_get_key_path() -> None:
    """Tests for :func:`_get_key_path`."""
    # Pass a string with a valid key path alias
    assertion.eq(_get_key_path('bond'), CP2K_KEYS_ALIAS['bond'])
    assertion.eq(_get_key_path('lennard-jones'), CP2K_KEYS_ALIAS['lennard-jones'])
    assertion.eq(_get_key_path('genpot14'), CP2K_KEYS_ALIAS['genpot14'])

    # Directly pass a key path
    ref = ('specific', 'cp2k', 'a', 'b', 'c')
    assertion.eq(_get_key_path(('a', 'b', 'c')), ref)
    assertion.eq(_get_key_path(('specific', 'cp2k', 'a', 'b', 'c')), ref)
    assertion.eq(_get_key_path(('input', 'a', 'b', 'c')), ref)

    # Don't know these keys
    assertion.assert_(_get_key_path, 'bob', exception=KeyError)
    assertion.assert_(_get_key_path, 5, exception=KeyError)


def test_validate_unit() -> None:
    """Tests for :func:`_validate_unit`."""
    columns = ['1', '2']

    # Valid inputs
    unit_iter1 = ['a', 'b']
    unit_iter2 = repeat('a')
    assertion.assert_(_validate_unit, unit_iter1, columns, invert=True)
    assertion.assert_(_validate_unit, unit_iter2, columns, invert=True)

    # *unit_iter3* is too short; should raise a LengthError
    unit_iter3 = ['a']
    assertion.assert_(_validate_unit, unit_iter3, columns, exception=LengthError)


def test_parse_unit() -> None:
    """Tests for :func:`_parse_unit`."""
    unit1 = None
    unit2 = 'kcalmol'
    unit3 = [None, None]
    unit4 = ['kcalmol', None]

    assertion.eq(_parse_unit(unit1), ['{}'])
    assertion.eq(_parse_unit(unit2), ['[kcalmol] {}'])
    assertion.eq(_parse_unit(unit3), ['{}', '{}'])
    assertion.eq(_parse_unit(unit4), ['[kcalmol] {}', '{}'])


def test_construct_df() -> None:
    """Tests for :func:`_construct_df`."""
    dict1 = {'Cs': (1, 1), 'Cd': (2, 2), 'O': (3, 3), 'H': (4, 4)}
    dict2 = {'Cs': 1, 'Cd': 2, 'O': 3, 'H': 4}
    dict3 = {'Cs': 1, 'Cd': 2, 'O': (3, 3), 'H': (4, 4)}
    dict4 = {k: set(v) for k, v in dict1.items()}
    c1 = ['a', 'b']
    c2 = 'a'

    ref1 = np.array([v for v in dict1.values()], dtype=str).astype(object)
    ref2 = np.array([v for v in dict2.values()], dtype=str).astype(object)[..., None]

    np.testing.assert_array_equal(_construct_df(c1, dict1), ref1)
    np.testing.assert_array_equal(_construct_df(c2, dict2), ref2)

    # No mixing of scalars and sequences allowed
    assertion.assert_(_construct_df, c1, dict2, exception=TypeError)
    assertion.assert_(_construct_df, c2, dict1, exception=TypeError)
    assertion.assert_(_construct_df, c1, dict3, exception=TypeError)
    assertion.assert_(_construct_df, c1, dict3, exception=TypeError)

    # The dict values are now either too long or too short
    assertion.assert_(_construct_df, [c2], dict1, exception=LengthError)
    assertion.assert_(_construct_df, [c2, c2, c2], dict1, exception=LengthError)

    # Collections (e.g. a set) are not allowed as dictionary values; requires a proper Sequence
    assertion.assert_(_construct_df, c1, dict4, exception=TypeError)


def test_set_prm() -> None:
    """Tests for :func:`set_prm`."""
    ref = Settings({'specific': {'cp2k': {'force_eval': {'mm': {'forcefield': {'nonbonded': {'lennard-jones': [{'epsilon': '[kcalmol] 1', 'sigma': '[angstrom] 1', 'atoms': 'Cs'}, {'epsilon': '[kcalmol] 2', 'sigma': '[angstrom] 2', 'atoms': 'Cd'}, {'epsilon': '[kcalmol] 3', 'sigma': '[angstrom] 3', 'atoms': 'O'}, {'epsilon': '[kcalmol] 4', 'sigma': '[angstrom] 4', 'atoms': 'H'}]}}}}}}})  # noqa

    key = 'lennard_jones'

    s1 = Settings()
    value1 = {
        'param': ('epsilon', 'sigma'),
        'unit': ('kcalmol', 'angstrom'),
        'Cs': (1, 1),
        'Cd': (2, 2),
        'O': (3, 3),
        'H': (4, 4)
    }

    s2 = Settings()
    value2 = [
        {'param': 'epsilon', 'unit': 'kcalmol', 'Cs': 1, 'Cd': 2, 'O': 3, 'H': 4},
        {'param': 'sigma', 'unit': 'angstrom', 'Cs': 1, 'Cd': 2, 'O': 3, 'H': 4}
    ]

    s3 = Settings()
    value3 = pd.DataFrame({
        'param': ('epsilon', 'sigma'),
        'unit': ('kcalmol', 'angstrom'),
        'Cs': (1, 1),
        'Cd': (2, 2),
        'O': (3, 3),
        'H': (4, 4)
    })

    s4 = Settings()
    value4 = [
        pd.Series({'param': 'epsilon', 'unit': 'kcalmol', 'Cs': 1, 'Cd': 2, 'O': 3, 'H': 4}),
        pd.Series({'param': 'sigma', 'unit': 'angstrom', 'Cs': 1, 'Cd': 2, 'O': 3, 'H': 4})
    ]

    set_prm(s1, key, value1, None)
    set_prm(s2, key, value2, None)
    set_prm(s3, key, value3, None)
    set_prm(s4, key, value4, None)
    assertion.eq(s1, ref)
    assertion.eq(s2, ref)
    assertion.eq(s3, ref)
    assertion.eq(s4, ref)


def test_prm_to_df() -> None:
    """Tests for :func:`prm_to_df`."""
    s = Settings()
    s.lennard_jones = {
        'param': ('epsilon', 'sigma'),
        'unit': ('kcalmol', 'angstrom'),
        'Cs': (1, 1),
        'Cd': (2, 2),
        'O': (3, 3),
        'H': (4, 4)
    }

    ref = {'param': {'epsilon': 'epsilon', 'sigma': 'sigma'},
           'unit': {'epsilon': 'kcalmol', 'sigma': 'angstrom'},
           'Cs': {'epsilon': 1.0, 'sigma': 1.0},
           'Cd': {'epsilon': 2.0, 'sigma': 2.0},
           'O': {'epsilon': 3.0, 'sigma': 3.0},
           'H': {'epsilon': 4.0, 'sigma': 4.0}}
    prm_to_df(s)
    assertion.eq(s['lennard_jones'].to_dict(), ref)

    # The 'param' key is missing
    s2 = {'lennard-jones': {'Cs': (1, 2)}}
    assertion.assert_(prm_to_df, s2, exception=KeyError)
