"""Default options to call some quantum packages."""

__all__ = ['freq', 'geometry', 'singlepoint', 'ts', 'md', 'cell_opt']

import yaml

from .._settings import Settings
from ..yaml_utils import UniqueSafeLoader


#: Templates for single-point calculations.
singlepoint = Settings(yaml.load("""
specific:
    adf:
        basis:
            type: SZ
        numericalquality:
            normal
        scf:
            converge: 1e-6
            iterations: 100

    ams:
        ams:
            Task: SinglePoint

    dftb:
        dftb:
            resourcesdir:
                "DFTB.org/3ob-3-1"

    cp2k:
        force_eval:
            dft:
                mgrid:
                    cutoff: 400
                    ngrids: 4
                print:
                    mo:
                        add_last: numeric
                        each:
                            qs_scf: 0
                        eigenvalues: ""
                        eigenvectors: ""
                        filename: "./mo.data"
                        ndigits: 36
                        occupation_numbers: ""
                qs:
                    method: gpw
                scf:
                    eps_scf: 1e-06
                    max_scf: 200
                    scf_guess: restart

            subsys:
                cell:
                    periodic: xyz
        global:
            print_level: low
            project: cp2k
            run_type: energy

    cp2k_mm:
        force_eval:
            method: FIST
            mm:
                print:
                    ff_info low:
                        spline_data: .FALSE.
                        spline_info: .FALSE.
                forcefield:
                    ei_scale14: 1.0
                    vdw_scale14: 1.0
                    ignore_missing_critical_params: ''
                    do_nonbonded: ''
                    shift_cutoff: .TRUE.
                    spline:
                        r0_nb: 0.25
                poisson:
                    periodic: NONE
                    ewald:
                        ewald_type: NONE
                        o_spline: 4
            subsys:
                cell:
                    periodic: NONE
                topology:
                    coord_file_format: 'OFF'
                    center_coordinates:
                        center_point: 0.0 0.0 0.0
        global:
            print_level: low
            project: cp2k
            run_type: energy

    orca:
        method:
            method: dft
            functional: lda
        basis:
            basis: sto_sz
""", Loader=UniqueSafeLoader))

#: Templates for geometry optimization calculations.
geometry = Settings(yaml.load("""
specific:
    adf:
        basis:
            type: SZ
        numericalquality: good
        scf:
            converge: 1e-6
            iterations: 100
        geometry:
            optim: delocal

    ams:
        ams:
            Task: GeometryOptimization
            GeometryOptimization:
                MaxIterations: 500

    dftb:
        dftb:
            resourcesdir: "DFTB.org/3ob-3-1"

    cp2k:
        motion:
            geo_opt:
                type: minimization
                optimizer: bfgs
                max_iter: 500
        force_eval:
            dft:
                mgrid:
                    cutoff: 400
                    ngrids: 4
                qs:
                    method: gpw
                scf:
                    eps_scf: 1e-06
                    max_scf: 200
                    scf_guess: atomic
                    OT:
                        minimizer: DIIS
                        N_DIIS: 7
                        preconditioner: full_single_inverse
            subsys:
                cell:
                    periodic: xyz
        global:
            print_level: low
            project: cp2k
            run_type: geometry_optimization

    cp2k_mm:
        force_eval:
            method: FIST
            mm:
                print:
                    ff_info low:
                        spline_data: .FALSE.
                        spline_info: .FALSE.
                forcefield:
                    ei_scale14: 1.0
                    vdw_scale14: 1.0
                    ignore_missing_critical_params: ''
                    do_nonbonded: ''
                    shift_cutoff: .TRUE.
                    spline:
                        r0_nb: 0.25
                poisson:
                    periodic: NONE
                    ewald:
                        ewald_type: NONE
                        o_spline: 4
            subsys:
                cell:
                    periodic: NONE
                topology:
                    coord_file_format: 'OFF'
                    center_coordinates:
                        center_point: 0.0 0.0 0.0
        motion:
            geo_opt:
                type: minimization
                optimizer: lbfgs
                max_iter: 500
        global:
            print_level: low
            project: cp2k
            run_type: geometry_optimization

    orca:
        method:
            method: dft
            functional: lda
            runtyp: opt
        basis:
            basis: sto_sz
""", Loader=UniqueSafeLoader))

#: Templates for transition state calculations.
ts = Settings(yaml.load("""
specific:
    adf:
        basis:
            type: SZ
        numericalquality: good
        scf:
            converge: 1e-6
            iterations: 100
        geometry:
            transitionstate: "mode=1"

    dftb:
        dftb:
            resourcesdir: "DFTB.org/3ob-3-1"

    cp2k: Null

    cp2k_mm: Null

    orca:
        method:
            method: dft
            functional: lda
            runtyp: opt
        geom:
            ts_search: ef
        basis:
            basis: sto_sz
""", Loader=UniqueSafeLoader))

#: Templates for frequency analyses calculations.
freq = Settings(yaml.load("""
specific:
    adf:
        analyticalfreq:
        basis:
            core: None
            type: DZP

        numericalquality: good

    ams:
        ams:
            Task: SinglePoint
            Properties:
                NormalModes: Yes

    dftb:
        dftb:
            resourcesdir: "DFTB.org/3ob-3-1"

    cp2k:
        force_eval:
            dft:
                mgrid:
                    cutoff: 400
                    ngrids: 4
                qs:
                    method: gpw
                scf:
                    eps_scf: 1e-06
                    max_scf: 200
                    scf_guess: restart
            subsys:
                cell:
                    periodic: xyz
        vibrational_analysis:
            thermochemistry: .TRUE.
        global:
            print_level: low
            project: cp2k
            run_type: VIBRATIONAL_ANALYSIS

    cp2k_mm:
        force_eval:
            method: FIST
            mm:
                print:
                    ff_info low:
                        spline_data: .FALSE.
                        spline_info: .FALSE.
                forcefield:
                    ei_scale14: 1.0
                    vdw_scale14: 1.0
                    ignore_missing_critical_params: ''
                    do_nonbonded: ''
                    shift_cutoff: .TRUE.
                    spline:
                        r0_nb: 0.25
                poisson:
                    periodic: NONE
                    ewald:
                        ewald_type: NONE
                        o_spline: 4
            subsys:
                cell:
                    periodic: NONE
                topology:
                    coord_file_format: 'OFF'
                    center_coordinates:
                        center_point: 0.0 0.0 0.0
        vibrational_analysis:
            thermochemistry: .TRUE.
        global:
            print_level: low
            project: cp2k
            run_type: VIBRATIONAL_ANALYSIS

    orca:
        method:
            method: dft
            functional: lda
        basis:
            basis: sto_sz
        main: freq
""", Loader=UniqueSafeLoader))

#: Templates for molecular dynamics (MD) calculations.
md = Settings(yaml.load("""
specific:
    adf: Null

    ams: Null

    cp2k: Null

    orca: Null

    cp2k_mm:
        force_eval:
            method: FIST
            mm:
                print:
                    ff_info low:
                        spline_data: .FALSE.
                        spline_info: .FALSE.
                forcefield:
                    ei_scale14: 1.0
                    vdw_scale14: 1.0
                    ignore_missing_critical_params: ''
                    do_nonbonded: ''
                    shift_cutoff: .TRUE.
                    spline:
                        r0_nb: 0.25
                poisson:
                    periodic: NONE
                    ewald:
                        ewald_type: NONE
                        o_spline: 4
            subsys:
                cell:
                    periodic: NONE
                topology:
                    coord_file_format: 'OFF'
                    center_coordinates:
                        center_point: 0.0 0.0 0.0
        motion:
            md:
                ensemble: NVT
                temperature: 300.0
                timestep: 1.0
                steps: 10000
                thermostat:
                    type: CSVR
                    csvr:
                        timecon: 1250
        global:
            print_level: low
            project: cp2k
            run_type: MD
""", Loader=UniqueSafeLoader))

#: Templates for cell optimization calculations.
cell_opt = Settings(yaml.load("""
specific:
    adf: Null

    ams: Null

    dftb: Null

    cp2k: Null

    orca: Null

    cp2k_mm:
        force_eval:
            stress_tensor: analytical
            method: FIST
            mm:
                print:
                    ff_info low:
                        spline_data: .FALSE.
                        spline_info: .FALSE.
                forcefield:
                    ei_scale14: 1.0
                    vdw_scale14: 1.0
                    ignore_missing_critical_params: ''
                    do_nonbonded: ''
                    shift_cutoff: .TRUE.
                    spline:
                        r0_nb: 0.25
                poisson:
                    periodic: xyz
                    ewald:
                        ewald_type: spme
                        o_spline: 4
            subsys:
                cell:
                    periodic: xyz
                topology:
                    coord_file_format: 'OFF'
                    center_coordinates:
                        center_point: 0.0 0.0 0.0
        motion:
            cell_opt:
                type: direct_cell_opt
                optimizer: lbfgs
                max_iter: 500
            print:
                cell low:
                    filename: ''
        global:
            print_level: low
            project: cp2k
            run_type: cell_opt
""", Loader=UniqueSafeLoader))
