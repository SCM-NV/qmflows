import sys
import importlib
from pathlib import Path
from typing import Any, Mapping, MutableMapping

import pytest
from qmflows.test_utils import PATH, HAS_RDKIT, requires_adf, requires_orca
from assertionlib import assertion

if sys.version_info >= (3, 7):
    from builtins import dict as OrderedDict
else:
    from collections import OrderedDict


class TestDeprecations:
    MODULES = OrderedDict(
        # qmflows.parsers
        xyzParser=("qmflows.parsers.xyzParser", "qmflows.parsers"),
        generic_parsers=("qmflows.parsers.generic_parsers", "qmflows.parsers"),
        parser=("qmflows.parsers.parser", "qmflows.parsers.utils"),
        orca_parser=("qmflows.parsers.orca_parser", "qmflows.parsers.orca"),
        cp2KParser=("qmflows.parsers.cp2KParser", "qmflows.parsers.cp2k"),
        adf_parser=("qmflows.parsers.adf_parser", "qmflows.parsers.adf"),

        # qmflows.templates
        templates=("qmflows.templates.templates", "qmflows.templates"),

        # qmflows.packages
        cp2k_mm=("qmflows.packages.cp2k_mm", "qmflows.packages"),
        cp2k_package=("qmflows.packages.cp2k_package", "qmflows.packages"),
        orca=("qmflows.packages.orca", "qmflows.packages"),
        package_wrapper=("qmflows.packages.package_wrapper", "qmflows.packages"),
        packages=("qmflows.packages.packages", "qmflows.packages"),
        serializer=("qmflows.packages.serializer", None),
        SCM=("qmflows.packages.SCM", "qmflows.packages"),

        # qmflows
        settings=("qmflows.settings", "qmflows"),
    )

    @pytest.mark.parametrize("name,alternative", MODULES.values(), ids=MODULES)
    def test_modules(self, name: str, alternative: "None | str") -> None:
        if alternative is None:
            pattern = f"`{name}` is deprecated and will be removed in the future"
        else:
            pattern = f"`{name}` is a deprecated alias for `{alternative}`"

        assertion.contains(sys.modules, name, invert=True)
        with pytest.warns(DeprecationWarning, match=pattern):
            importlib.import_module(name)

    FUNCTIONS = OrderedDict(
        # qmflows.examples
        readCp2KBasis=(
            "qmflows.parsers.cp2k.readCp2KBasis",
            "qmflows.parsers.cp2k.read_cp2k_basis",
            {"file": PATH / "BASIS_MOLOPT"},
        ),

        # qmflows.packages
        registry=("qmflows.packages.registry", None, {}),
        load_properties=("qmflows.packages.load_properties", None, {"name": "CP2K"}),
        SerMolecule=("qmflows.packages.SerMolecule", None, {}),
        SerSettings=("qmflows.packages.SerSettings", None, {}),
    )

    @pytest.mark.parametrize("full_name,alternative,kwargs", FUNCTIONS.values(), ids=FUNCTIONS)
    def test_functions(
        self,
        full_name: str,
        alternative: "None | str",
        kwargs: Mapping[str, Any],
    ) -> None:
        module_name, name = full_name.rsplit(".", 1)
        if alternative is None:
            pattern = f"`{full_name}` is deprecated and will be removed in the future"
        else:
            pattern = f"`{full_name}` is a deprecated alias for `{alternative}`"

        with pytest.warns(DeprecationWarning, match=pattern):
            module = importlib.import_module(module_name)
            func = getattr(module, name)
            func(**kwargs)

    FUNCTIONS_ADF = OrderedDict(
        example_partial_geometry_opt=(
            "qmflows.example_partial_geometry_opt",
            "qmflows.examples.example_partial_geometry_opt",
            {"folder": "example_partial_geometry_opt"},
        ),
        example_H2O2_TS=(
            "qmflows.example_H2O2_TS",
            "qmflows.examples.example_H2O2_TS",
            {"folder": "example_H2O2_TS"},
        ),
        example_generic_constraints=(
            "qmflows.example_generic_constraints",
            "qmflows.examples.example_generic_constraints",
            {"folder": "example_generic_constraints"},
        ),
    )

    @pytest.mark.slow
    @pytest.mark.skipif(not HAS_RDKIT, reason="requires RDKit")
    @requires_adf
    @pytest.mark.parametrize(
        "full_name,alternative,kwargs", FUNCTIONS_ADF.values(), ids=FUNCTIONS_ADF
    )
    def test_functions_adf(
        self,
        tmp_path: Path,
        full_name: str,
        alternative: "None | str",
        kwargs: MutableMapping[str, Any],
    ) -> None:
        kwargs["path"] = tmp_path

        module_name, name = full_name.rsplit(".", 1)
        if alternative is None:
            pattern = f"`{full_name}` is deprecated and will be removed in the future"
        else:
            pattern = f"`{full_name}` is a deprecated alias for `{alternative}`"

        with pytest.warns(DeprecationWarning, match=pattern):
            module = importlib.import_module(module_name)
            func = getattr(module, name)
            func(**kwargs)

    FUNCTIONS_ADF_ORCA = OrderedDict(
        example_freqs=(
            "qmflows.example_freqs",
            "qmflows.examples.example_freqs",
            {"folder": "example_freqs"},
        ),
    )

    @pytest.mark.slow
    @pytest.mark.skipif(not HAS_RDKIT, reason="requires RDKit")
    @requires_adf
    @requires_orca
    @pytest.mark.parametrize(
        "full_name,alternative,kwargs", FUNCTIONS_ADF_ORCA.values(), ids=FUNCTIONS_ADF_ORCA
    )
    def test_functions_adf_orca(
        self,
        tmp_path: Path,
        full_name: str,
        alternative: "None | str",
        kwargs: MutableMapping[str, Any],
    ) -> None:
        kwargs["path"] = tmp_path

        module_name, name = full_name.rsplit(".", 1)
        if alternative is None:
            pattern = f"`{full_name}` is deprecated and will be removed in the future"
        else:
            pattern = f"`{full_name}` is a deprecated alias for `{alternative}`"

        with pytest.warns(DeprecationWarning, match=pattern):
            module = importlib.import_module(module_name)
            func = getattr(module, name)
            func(**kwargs)
