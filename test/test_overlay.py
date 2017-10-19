from nose.plugins.attrib import attr
from qmflows import (Settings, templates)


@attr('fast')
def test_overlay_cp2k_singlepoint():
    """
    Test if the merging with the CP2K templates produces the same Settings that
    writing it by hand.
    """
    s = Settings()
    dft = s.specific.cp2k.force_eval.dft
    force = s.specific.cp2k.force_eval
    dft.scf.scf_guess = 'atomic'
    dft.scf.ot.minimizer = 'diis'
    dft.scf.ot.n_diis = 7
    dft.scf.ot.preconditioner = 'full_single_inverse'
    dft.scf.added_mos = 0
    dft.scf.eps_scf = 5e-06

    r = templates.singlepoint.overlay(s)

    dft.basis_set_file_name = ''
    dft.mgrid.cutoff = 400
    dft.mgrid.ngrids = 4

    dft.potential_file_name = ''
    dft['print']['mo']['add_last'] = 'numeric'
    dft['print']['mo']['each']['qs_scf'] = 0
    dft['print']['mo']['eigenvalues'] = ''
    dft['print']['mo']['eigenvectors'] = ''
    dft['print']['mo']['filename'] = './mo.data'
    dft['print']['mo']['ndigits'] = 36
    dft['print']['mo']['occupation_numbers'] = ''

    dft.scf.max_scf = 200
    dft.scf.scf_guess = 'atomic'
    dft.xc.xc_functional = 'pbe'

    dft.qs.method = 'gpw'
    force.subsys.cell.periodic = 'xyz'

    g = s['specific']['cp2k']['global']
    g.print_level = 'low'
    g.project = "qmflows-cp2k"
    g.run_type = "energy_force"

    print(s.specific.cp2k.force_eval)
    print(r.specific.cp2k.force_eval)

    assert s.specific.cp2k == r.specific.cp2k


@attr('fast')
def test_overlay_adf_freq():
    """
    Test if the overlay function produces the same Settings that a user would
    write by hand.
    """
    s = Settings()
    s.specific.adf.xc.gga = "bp86"
    s.specific.adf.beckegrid.quality = "good"
    s.specific.adf.scf.iterations = "99"
    s.specific.adf.scf.converge = "0.0000001"
    s.specific.adf.analyticalfreq

    r = templates.freq.overlay(s)

    s.specific.adf.basis.type = "DZP"
    s.specific.adf.basis.core = "None"

    print(s.specific.adf)
    print(r.specific.adf)
    assert s.specific.adf == r.specific.adf
