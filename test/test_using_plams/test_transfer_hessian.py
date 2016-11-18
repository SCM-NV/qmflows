from nose.plugins.attrib import attr
from qmworks import Settings, templates, molkit
from noodles import gather

# User Defined imports
from qmworks.packages.SCM import dftb
from qmworks.packages.orca import orca
from qmworks.packages import run


@attr('slow')
def test_hessian_transfer():
    """
    Test DFTB -> Orca hessian transfer
    """
    h2o = molkit.from_smiles('O')
    h2o.properties.symmetry = 'C1'

    h2o_freq = dftb(templates.freq, h2o, job_name="freq").hessian

    s = Settings()
    s.inithess = h2o_freq

    h2o_opt = orca(templates.geometry.overlay(s), h2o, job_name="opt")

    energy = h2o_opt.energy
    dipole = h2o_opt.dipole

    wf = gather(energy, dipole)
    # draw_workflow('wf.svg', wf._workflow)

    print(run(wf))
