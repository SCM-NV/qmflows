__all__ = ['example_H2O2_TS']

from qmflows import (orca, dftb, templates, Dihedral, run, Settings)
import scm.plams.interfaces.molecule.rdkit as molkit


def example_H2O2_TS():
    """
    This example generates an approximate TS for rotation in hydrogen peroxide
    using DFTB, and performs a full TS optimization in Orca.
    It illustrates using a hessian from one package, DFTB in this case,
    to initialize a TS optimization in another, i.e. Orca in this case
    """

    # Generate hydrogen peroxide molecule
    h2o2 = molkit.from_smarts('[H]OO[H]')

    # Define dihedral angle
    dihe = Dihedral(1, 2, 3, 4)

    # Constrained geometry optimization with DFTB
    # The dihedral is set to 0.0 to obtain an approximate TS
    s1 = Settings()
    s1.constraint.update(dihe.get_settings(0.0))
    dftb_opt = dftb(templates.geometry.overlay(s1), h2o2)

    # Calculate the DFTB hessian
    dftb_freq = dftb(templates.freq, dftb_opt.molecule)

    # Transition state calculation using the DFTB hessian as starting point
    s2 = Settings()
    s2.inithess = dftb_freq.hessian
    orca_ts = orca(templates.ts.overlay(s2), dftb_opt.molecule)

    # Execute the workflow
    result = run(orca_ts)

    # Analyse the result
    ts_dihe = round(dihe.get_current_value(result.molecule))
    n_optcycles = result.optcycles

    print('Dihedral angle (degrees): {:.0f}'.format(ts_dihe))
    print('Number of optimization cycles: {:d}'.format(n_optcycles))

    return ts_dihe, n_optcycles


if __name__ == "__main__":
    example_H2O2_TS()
