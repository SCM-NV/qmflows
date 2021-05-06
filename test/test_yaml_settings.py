"""Test the conversion from yaml to settings."""

from typing import Type

import yaml
import pytest
from assertionlib import assertion
from yaml.constructor import BaseConstructor, ConstructorError

from qmflows import Settings
from qmflows.fileFunctions import yaml2Settings
from qmflows.yaml_utils import UniqueFullLoader, UniqueLoader, UniqueSafeLoader, UniqueUnsafeLoader

CP2K_PBE_GUESS_DUPLICATE = """
specific:
    cp2k:
        global:
            run_type:
                energy
        force_eval:
            subsys:
                cell:
                    periodic: "None"
                    periodic: "None"
            dft:
                xc:
                    xc_functional pbe: {}
                scf:
                    eps_scf: 1e-6
"""

CP2K_PBE_GUESS = """
specific:
    cp2k:
        global:
            run_type:
                energy
        force_eval:
            subsys:
                cell:
                    periodic: "None"
            dft:
                xc:
                    xc_functional pbe: {}
                scf:
                    eps_scf: 1e-6
"""

REF = Settings()
REF.specific.cp2k.force_eval.dft.xc["xc_functional pbe"] = {}
REF.specific.cp2k.force_eval.subsys.cell.periodic = "None"
REF.specific.cp2k["global"]["run_type"] = "energy"
REF.specific.cp2k.force_eval.dft.scf.eps_scf = "1e-6"

LOADER_TYPES = (UniqueFullLoader, UniqueLoader, UniqueSafeLoader, UniqueUnsafeLoader)


def test_yaml2Settings():
    """Test the conversion from yaml to settings."""
    s1 = yaml2Settings(CP2K_PBE_GUESS)
    s2 = yaml2Settings(CP2K_PBE_GUESS.encode())

    assertion.eq(s1.specific, REF.specific)
    assertion.eq(s2.specific, REF.specific)


class TestLoader:
    """Tests for the :mod:`qmflows.yaml_utils` loaders."""

    @pytest.mark.parametrize("loader", LOADER_TYPES)
    def test_pass(self, loader: Type[BaseConstructor]) -> None:
        """Test for succesful :func:`yaml.load` calls."""
        dct = yaml.load(CP2K_PBE_GUESS, Loader=loader)
        assertion.eq(dct["specific"], REF.specific)

    @pytest.mark.parametrize("loader", LOADER_TYPES)
    def test_raise(self, loader: Type[BaseConstructor]) -> None:
        """Test for failed :func:`yaml.load` calls."""
        with pytest.raises(ConstructorError):
            yaml.load(CP2K_PBE_GUESS_DUPLICATE, Loader=loader)
