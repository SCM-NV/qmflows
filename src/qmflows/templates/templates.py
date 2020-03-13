"""Default options to call some quantum packages."""

__all__ = ['freq', 'geometry', 'singlepoint', 'ts', 'md']

import yaml

from ..settings import Settings


singlepoint = Settings(yaml.load("""
specific:
    adf:
        basis:
            type: SZ
        xc:
            __block_replace: true
            lda: ""
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
        task:
            runtype: SP

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
            subsys:
                cell:
                    abc: '[angstrom] 50.0 50.0 50.0'
                    periodic: NONE
                topology:
                    coord_file_format: 'OFF'
                    center_coordinates:
                        center_point: 0.0 0.0 0.0
        global:
            print_level: low
            project: cp2k
            run_type: energy

    gamess:
        basis:
            gbasis: sto
            ngauss: 3
        contrl:
            scftyp: rhf
            dfttyp: pbe

    orca:
        method:
            method: dft
            functional: lda
        basis:
            basis: sto_sz
""", Loader=yaml.FullLoader))

geometry = Settings(yaml.load("""
specific:
    adf:
        basis:
            type: SZ
        xc:
            __block_replace: true
            lda: ""
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
        task:
            runtype: GO
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
            subsys:
                cell:
                    abc: '[angstrom] 50.0 50.0 50.0'
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

    gamess:
        basis:
            gbasis: n21
            ngauss: 3
        contrl:
            scftyp: rhf
            dfttyp: pbe
            runtyp: optimize
""", Loader=yaml.FullLoader))

# Transition state template
ts = Settings(yaml.load("""
specific:
    adf:
        basis:
            type: SZ
        xc:
            "__block_replace": true
            lda: ""
        numericalquality: good
        scf:
            converge: 1e-6
            iterations: 100
        geometry:
            transitionstate: "mode=1"

    dftb:
        task:
            runtype: TS
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
""",    Loader=yaml.FullLoader))

# Frequencies template
freq = Settings(yaml.load("""
specific:
    adf:
        analyticalfreq:
        basis:
            core: None
            type: DZP

        xc:
            "__block_replace": true
            lda: ""
        numericalquality: good

    ams:
        ams:
            Task: SinglePoint
            Properties:
                NormalModes: Yes

    dftb:
        task:
            runtype: Frequencies
        dftb:
            resourcesdir: "DFTB.org/3ob-3-1"

    cp2k : Null

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
            subsys:
                cell:
                    abc: '[angstrom] 50.0 50.0 50.0'
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

    gamess:
        basis:
            gbasis: n21
            ngauss: 3
        contrl:
            scftyp: rhf
            dfttyp: pbe
            runtyp: hessian
            nosym: 1
        force:
            method: seminum
""", Loader=yaml.FullLoader))


# moleculair dynamics template
md = Settings(yaml.load("""
specific:
    adf: Null

    ams: Null

    cp2k: Null

    orca: Null

    gamess: Null

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
            subsys:
                cell:
                    abc: '[angstrom] 50.0 50.0 50.0'
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
""", Loader=yaml.FullLoader))
