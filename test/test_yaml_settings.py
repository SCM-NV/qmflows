"""Test the conversion from yaml to settings."""
from assertionlib import assertion
from qmflows import Settings
from qmflows.fileFunctions import yaml2Settings


cp2k_pbe_guess = """
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


def test_yaml2Settings():
    """Test the conversion from yaml to settings."""
    s1 = yaml2Settings(cp2k_pbe_guess)
    s2 = yaml2Settings(cp2k_pbe_guess.encode())

    ref = Settings()
    ref.specific.cp2k.force_eval.dft.xc["xc_functional pbe"] = {}
    ref.specific.cp2k.force_eval.subsys.cell.periodic = "None"
    ref.specific.cp2k["global"]["run_type"] = "energy"
    ref.specific.cp2k.force_eval.dft.scf.eps_scf = "1e-6"
    assertion.eq(s1.specific, ref.specific)
    assertion.eq(s2.specific, ref.specific)
