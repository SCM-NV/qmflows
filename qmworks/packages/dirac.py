
__all__ = ['dirac', 'DIRAC', 'DIRAC_Result']


# =======>  Standard and third party Python Libraries <======
from functools import reduce
from noodles import (Storable)
from warnings import warn
import plams

# ==================> Internal modules <====================
from qmworks.packages.packages import Package
from qmworks.settings import Settings

# ==================> <======================


class DIRAC(Package):
    """
    """
    def __init__(self):
        super().__init__("dirac")
        self.generic_dict_file = 'generic2DIRAC.json'

    def prerun(self):
        pass

    @staticmethod
    def run_job(settings, mol, job_name="dirac_job"):

        dirac_settings = Settings()
        dirac_settings.input = settings.specific.dirac
        # check_dirac_input(dirac_settings)
        result = plams.DiracJob(name=job_name, settings=dirac_settings,
                                molecule=mol).run()

        return DIRAC_Result(dirac_settings, mol, result.job.name,
                            plams_dir=result.job.path)

    def postrun(self):
        pass

    @staticmethod
    def handle_special_keywords(settings, key, value, mol):
        """
        Create the settings input for complex Dirac keywords
        """
        warn('Keyword ' + key + ' doesn\'t exist')

# Instance
dirac = DIRAC()


class DIRAC_Result(Storable):
    """
    Class to access **DIRAC** Results.
    """
    def __init__(self, settings, molecule, job_name, plams_dir, project_name=None):
        properties = 'data/dictionaries/propertiesDIRAC.json'
        super().__init__(settings, molecule, job_name=job_name, plams_dir=plams_dir,
                         project_name=project_name, properties=properties)

    def as_dict(self):
        return {
            "settings": self.settings,
            "molecule": self._molecule,
            "filename": self.result.path}

    @property
    def molecule(self):
        pass

    @classmethod
    def from_dict(cls, settings, molecule, job_name, archive, project_name):
        plams_dir = archive["plams_dir"].path
        return DIRAC_Result(settings, molecule, job_name, plams_dir, project_name)


def check_dirac_input(s):
    """
    hamiltonian:
    Possible choices are:
     - DC:    Dirac-Coulomb
     - DCG:   Dirac-Coulomb-Gaunt (available for DFT, HF)
     - MMF:   Molecular Mean Field (based on DCG Hamiltonian)
     - SFDC:  Spin-free Dirac-Coulomb
     - X2C:   eXact 2-Component
     - ECP:   Relativistic (1- or 2-component) ECP
     - Levy:  Levy-Leblond (NR limit of Dirac equation)
     - Nonr:  Non-relativistic (true 1 component) Hamiltonian
     :param hamiltonian: A string identifying the Hamiltonian.
     :type  hamiltonian: str

    method:
    Available options are: HF, DFT, MP2, CCSD, CCSDt, (FSCC, IHFSCC)

     :param s.input.method: string identifying the selected method
     :type  s.input.method: str

   properties:
    :param s.input.properties:
       A list of the properties to calculate. Possible choices
       are 'dipole', 'efg', 'nqcc'. See Dirac documentation of
       **PROPERTIES** for more.
    :type  s.input.properties: list

    functional:
    The exchange-correlation functional for DFT.

        :param s.input.functional:
            A string identifying the functional.
            See Dirac manual for available options.
        :type s.input.functional: str

    dftgrid:
    the numerical integration grid.

     Switch on/off export of density and Coulomb potential for FDE.

      :param s.input.export: Whether to export the density and Coulomb potential.
      :type  s.input.export: bool

   :domoltra:
      :param s.input.domoltra: Switch on transformation.
      :type  s.input.domoltra: bool

   :moltra:
    set options for 4-index transformation.
    :param s.input.moltra: either ['all'] or a list giving minimum and maximum orbital energy
                           and the degeneracy treshold (three numbers)
    :type s.input.moltra: list

   :nucmod:
    Choose the nuclear model used in the Dirac calculation.
    :param s.input.nucmod: 'finite' or None for finite nuclei (default) or 'point' for point nuclei
    :type  s.input.nucmod: String

    """
    pass
    # for k in s.input.keys():
    #     val = s.input[k]
    #     key = k.upper()
    #     if isinstance(val, list):
    #         val = list([x.upper() for x in val])
    #     elif isinstance(val, str):
    #         val = val.upper()
    #     s.input[key] = val

    # return check_easy_input(s)


def check_easy_input(easySett):
    """
    Checks the user input validity and translate the easy input to the
    real Dirac input structure
    :param easySett: Settings containing the user provided input
    :type  easySett: Settings
    """
    return reduce(lambda acc, f: f(acc),
                  [check_hamiltonian, check_method, check_properties,
                   check_transform, check_nucmod], easySett)


def check_hamiltonian(s):
    """
    Checks the Selected the Hamiltonian, store in ``s.input.hamiltonian``
    :param s: Settings containing the input tree
    :type  s: Settings

    """
    supported_hamiltonians = {
        'DC': 'Dirac-Coulomb',
        'DCG': 'Dirac-Coulomb-Gaunt (available for DFT, HF)',
        'SFDC': 'Spin-free Dirac-Coulomb',
        'MMF': 'DCG-based Molecular Mean Field',
        'X2C': 'eXact 2-Component',
        'ECP': 'Relativistic ECP',
        'LEVY': 'Levy-Leblond (NR limit of Dirac equation)',
        'NONR': 'Non-relativistic (true 1 component) Hamiltonian'
    }

    ham = s.input.HAMILTONIAN
    if ham not in supported_hamiltonians:
        err = 'Dirac does not support Hamiltonian:{}\
        \nSupported Hamiltonians:{}\n'.format(ham, supported_hamiltonians)
        raise RuntimeError(err)

    return s


def check_method(s):
    """
    Checks the user input list of properties stored ``in S.input.properties``.
    :param s: Settings containing the input tree
    :type  s: Settings
    """
    met = s.input.METHOD
    if any([met == x for x in ['HF', 'DFT']]):
        s.input.EXPORTFDE_LEVEL = 'DHF'
    elif any([met == x for x in ['MP2', 'CCSD', 'CCSDt', 'FSCC', 'IHFSCC']]):
        s.input.EXPORTFDE_LEVEL = 'MP2'
        s.input.DOMOLTRA = True

    return s


def check_transform(s):
    """
    Checks if the ``domoltra`` is active otherwise add some defaults.
    :param s: Settings containing the input tree
    :type  s: Settings
    """
    if s.input.DOMOLTRA:
        moltra = s.input.get('MOLTRA')
        if not moltra:
            s.input.MOLTRA = [-5.0, 10.0, 0.1]
            xs  = s.input.MOLTRA.strip('[]')
            msg = 'Default moltra parameters are:{}'.format(xs)
            raise RuntimeError(msg)
    return s


def check_properties(s):
    """
    Checks the user input list of properties stored ``in s.input.properties``.
    Possible choices  are 'dipole', 'efg', 'nqcc'. See Dirac documentation of
        **PROPERTIES for more.
    :param s: Settings containing the input tree
    :type  s: Settings

    """
    predifined_properties = ['DIPOLE', 'QUADRUPOLE', 'EFG', 'NQCC',
                             'POLARIZABILITY', 'FIRST ORDER HYPERPOLARIZABILITY',
                             'VERDET', 'TWO-PHOTON', 'NMR', 'SHIELDING',
                             'MAGNET', 'SPIN-SPIN COUPLING', 'DSO',
                             'NSTDIAMAGNETIC', 'MOLGRD', 'PVC', 'RHONUC',
                             'EFFDEN']

    def iselem(x):
        if x not in predifined_properties:
            err = 'unkown property:{}\n'.format(x)
            raise RuntimeError(err)

    properties = s.input.get('PROPERTIES')
    for p in properties:
        iselem(properties)

    return s


def check_nucmod(s):
    """
    Checks the ``nuclear model`` input validity.
     :param s: Settings containing the input tree
     :type  s: Settings

    """
    nucmod = s.input.NUCMOD

    if not nucmod:
        if any(nucmod == x for x in ['FINITE', 'POINT']):
            err = 'Unkown nuclear mode:{}\n'.format(nucmod)
            raise RuntimeError(err)

    return s


# ======================> <==================================
def build_dirac_input(s):
    """
    Transform the Setting provided by the user together with the defaults into
    a Dirac input.
     :param easySett:  Settings containing the user provided input
     :type  easySett:  Settings

    """
    inp = Settings()
    inp.input.DIRAC["WAVE FUNCTION"]
    if s.input.GEOMOPT:
        inp.input.DIRAC["WAVE FUNCTION"]["OPTIMIZE"]

    funs = [build_hamiltonian_opts, build_method_opts, build_basis_opts,
            build_transf_opts, build_integral_opts, build_properties_opts,
            build_exportlevel_opts]

    for f in funs:
        inp, s = f(*(inp, s))

    return inp


def build_hamiltonian_opts(inp, s):
    """
    Adds Hamiltonian to input structure together with some special keywords
    related to the Hamiltonians.
    :parameter inp: Real Dirac Input structure
    :type      inp: Settings
    :parameter   s: generic keywords
    :type        s: Settings

    """
    _special_names = {
        'DCG': ['DCG', 'GAUNT'], 'MMF': ['X2Cmmf', 'GAUNT'],
        'SFDC': ['SPINFREE'], 'LEVY': ['Levy-LEBLOND', 'LVCORR'],
        'NONR': ['NONREL', 'LVCORR']
    }

    ham   = s.input.HAMILTONIAN
    xs    = _special_names.get(ham)
    met   = s.input.METHOD
    funct = s.input.FUNCTIONAL
    if xs:
        for x in xs:
            inp.input.HAMILTONIAN[x]
    elif met == 'DFT':
        inp.input.HAMILTONIAN.DFT = funct
    else:
        inp.input.HAMILTONIAN[ham]

    return inp, s


def build_method_opts(inp, s):
    """
    Adds the method keywords to the input settings.
    :parameter inp: Real Dirac Input structure
    :type      inp: Settings
    :parameter   s: generic keywords
    :type        s: Settings
    """
    cc_family = ['CCSD', 'CCSDt', 'FSCC', 'IHFSCC']
    met       = s.input.METHOD

    inp.input["WAVE FUNCTION"]["SCF"]
    if met in cc_family:
        inp.input["WAVE FUNCTION"]["RELCCSD"]
    elif met == 'MP2':
        inp.input["WAVE FUNCTION"]["MP2"]

    return inp, s


def build_basis_opts(inp, s):
    """
    Adds Basis set opts
    :parameter inp: Real Dirac Input structure
    :type      inp: Settings
    :parameter   s: generic keywords
    :type        s: Settings
    """

    basis = s.input.BASIS
    inp.input["MOLECULE"]["BASIS"]["DEFAULT"] = basis
    return inp, s


def build_transf_opts(inp, s):
    """
    Options to Specify the active set of spinors in th integral transformation module.
    :parameter inp: Real Dirac Input structure
    :type      inp: Settings
    :parameter   s: generic keywords
    :type        s: Settings

    """
    moltra = s.input.MOLTRA
    if s.input.DOMOLTRA:
        if moltra[0] == "all":
            inp.input.MOLTRA.ACTIVE = "all"
        else:
            inp.input.MOLTRA.ACTIVE = moltra

    return inp, s


def build_integral_opts(inp, s):
    """
    Nuclear model specification.
    :parameter inp: Real Dirac Input structure
    :type      inp: Settings
    :parameter   s: generic keywords
    :type        s: Settings

    """
    nucmod = s.input.NUCMOD
    if nucmod == "POINT":
        inp.input.INTEGRALS.NUCMOD = 1

    return inp, s


def build_properties_opts(inp, s):
    """
    Adds Molecular property calculation to the input.
    :parameter inp: Real Dirac Input structure
    :type      inp: Settings
    :parameter   s: generic keywords
    :type        s: Settings
    """
    ps = s.input.PROPERTIES
    for p in ps:
        inp.input.PROPERTIES[p]
    return inp, s


def build_exportlevel_opts(inp, s):
    """
    """
    return inp, s
