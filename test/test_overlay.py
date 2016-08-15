
from qmworks import (Settings, templates)


def test_overlay_cp2k():
    """
    """
    s = Settings()
    dft = s.specific.cp2k.force_eval.dft
    force = s.specific.cp2k.force_eval
    dft.scf.scf_guess = 'atomic'
    dft.scf.ot.minimizer = 'DIIS'
    dft.scf.ot.n_diis = 7
    dft.scf.ot.preconditioner = 'FULL_SINGLE_INVERSE'
    dft.scf.added_mos = 0
    dft.scf.eps_scf = 5e-06

    r = templates.singlepoint.overlay(s)

    dft.basis_set_file_name = ''
    dft.mgrid.cutoff = 400
    dft.mgrid.ngrids = 4

    dft.potential_file = ''
    dft['print']['mo']['add_last'] = 'NUMERIC'
    dft['print']['mo']['each']['qs_scf'] = 0
    dft['print']['eigenvalues'] = ''
    dft['print']['eigenvectors'] = ''
    dft['print']['filename'] = './mo.data'
    dft['print']['ndigits'] = 36
    dft['print']['occupation_numbers'] = ''

    dft.scf.max_scf = 200
    dft.scf.scf_guess = 'atomic'
    dft.xc.xc_functional = 'PBE'

    dft.qs.method = 'GPW'
    force.subsys.cell.periodic = 'XYZ'
    force.subsys.topology.coordinate = 'XYZ'
    force.subsys.topology.coord_file_name = ''

    g = s['specific']['cp2k']['global']
    g.print_level = 'LOW'
    g.project = "QMWORKS-CP2K"
    g.run_type = "ENERGY_FORCE"

    
    print(s.specific.cp2k.force_eval)
    print(r.specific.cp2k.force_eval)
    
    assert s.specific.cp2k.dft == r.specific.cp2k.dft
