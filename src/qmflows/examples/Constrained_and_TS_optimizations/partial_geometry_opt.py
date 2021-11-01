"""A module containing the example :func:`example_partial_geometry_opt` function."""

__all__ = ['example_partial_geometry_opt']

import scm.plams.interfaces.molecule.rdkit as molkit
from noodles import gather
from qmflows import (adf, run, Settings, templates)


def example_partial_geometry_opt(*args, n_processes=1, **kwargs):
    """Performa partial optimization freezing the Hydrogen atoms."""
    methanol = molkit.from_smiles('CO')

    # optimize only H atoms
    s = Settings()
    s.freeze = [1, 2]
    geom_job1 = adf(templates.geometry.overlay(s), methanol, job_name='geom_job1').molecule

    # optimize only H atoms
    s = Settings()
    s.selected_atoms = ['H']
    geom_job2 = adf(templates.geometry.overlay(s), methanol, job_name='geom_job2').molecule

    geom1, geom2 = run(gather(geom_job1, geom_job2), *args, n_processes=n_processes, **kwargs)

    return geom1, geom2
