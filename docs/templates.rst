Templates
---------
The input generations consist of two parts: chosing a template (see templates_)
for the kind of calculation to perform and adding some settings to that template. Notice
that the user can either pick a specific package template or provides only generic
keywords.

.. _templates:

YAML
~~~~
.. currentmodule:: qmflows.templates.templates

The YAML_ markdown format is used together with the :mod:`yaml` module to implement the mechanism to load the templates.

.. _YAML: https://pyyaml.org/wiki/PyYAMLDocumentation

For Example, the default parameter for a single point calculation for several package is:

.. code-block:: python

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

      dirac:
        DIRAC: WAVEFUNCTION
        HAMILTONIAN: "LEVY-LEBLOND"
        WAVE FUNCTION: SCF

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
