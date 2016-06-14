
# =======>  Standard and third party Python Libraries <======
from noodles import files
import pkg_resources as pkg
import plams

# ==================> Internal Modules  <=====================
from qmworks.settings import Settings
from qmworks.packages.packages import Package, Result  # ChemResult
from qmworks.fileFunctions import json2Settings
# ========================= ADF ============================


class ADF(Package):
    """
    This class takes care of calling the *ADF* quantum package.
    it uses both Plams and the Templates module to create the input
    files, while Plams takes care of running the Job.
    It returns a ``ADF_Result`` object containing the output data.
    """
    def __init__(self):
        super(ADF, self).__init__("adf")
        self.generic_dict_file = 'generic2ADF.json'

    def prerun(self):
        pass

    def run_job(self, settings, mol, job_name='ADFjob'):
        """
        Execute ADF job.

        :param settings: user input settings.
        :type settings: |Settings|
        :param mol: Molecule to run the simulation
        :type mol: Plams Molecule
        :parameter input_file_name: The user can provide a name for the
                                    job input.
        :type input_file_name: String
        :parameter out_file_name: The user can provide a name for the
                                  job output.
        :type out_file_name: String
        :returns: :class:`~qmworks.packages.SCM.ADF_Result`
        """
        adf_settings = Settings()
        adf_settings.input = settings.specific.adf
        result = plams.ADFJob(name=job_name, molecule=mol,
                              settings=adf_settings).run()

        return ADF_Result(adf_settings, mol, result._kf)

    def postrun(self):
        pass

    def handle_special_keywords(self, settings, key, value, mol):
        """
        some keywords provided by the user do not have a straightforward
        translation to *ADF* input and require some hooks that handles the
        special behaviour of the following keywords:

        * ``freeze``
        * ``selected_atoms``
        """
        if key == "freeze":
            settings.specific.adf.geometry.optim = "cartesian"
            for a in value:
                settings.specific.adf.constraints['atom ' + str(a)] = ""
        elif key == "selected_atoms":
            settings.specific.adf.geometry.optim = "cartesian"
            if not isinstance(value, list):
                msg = 'selected_atoms ' + str(value) + ' is not a list'
                raise RuntimeError(msg)
            if isinstance(value[0], int):
                for a in range(1, len(mol) + 1):
                    if a not in value:
                         settings.specific.adf.constraints['atom ' + str(a)] = ""
            else:
                for a in range(len(mol)):
                    if mol.atoms[a].symbol not in value:
                        settings.specific.adf.constraints['atom ' + str(a+1)] = ""
        else:
            raise RuntimeError('Keyword ' + key + ' doesn\'t exist')


class ADF_Result(Result):
    """Class providing access to PLAMS ADFJob result results"""

    def __init__(self, settings, molecule, result):
        self.settings = settings
        self._molecule = molecule
        self.result = result
        properties = 'data/dictionaries/propertiesADF.json'
        xs = pkg.resource_string("qmworks", properties)
        self.prop_dict = json2Settings(xs)
        self.archive = files.Path(result.path)

    def as_dict(self):
        """
        Method to serialize as a JSON dictionary the results given
        by an ADF computation.
        """
        return {
            "settings": self.settings,
            "molecule": self._molecule,
            "filename": self.result.path}

    @classmethod
    def from_dict(cls, settings, molecule, filename):
        return ADF_Result(settings, molecule, plams.kftools.KFFile(filename))

    def get_property(self, prop, section=None):
        return self.result.read(section, prop)

    def __getattr__(self, prop):
        """Returns a section of the results.

        Example:

        ..

            dipole = result.properties.dipole

        """
        section, property = self.prop_dict[prop]
        return self.result.read(section, property)

    @property
    def molecule(self, unit='bohr', internal=False, n=1):
        """WARNING: Cheap copy from PLAMS, do not keep this!!!"""
        m = self._molecule.copy()
        natoms = len(m)
        # Find out correct location
        coords = self.result.read('Geometry', 'xyz InputOrder')
        coords = [coords[i:i + 3] for i in range(0, len(coords), 3)]

        if len(coords) > natoms:
            coords = coords[(n - 1) * natoms: n * natoms]
        if internal:
            mapping = self._int2inp()
            coords = [coords[mapping[i] - 1] for i in range(len(coords))]
        for at, coord in zip(m, coords):
            at.move_to(coord, unit)

        return m


class DFTB(Package):
    """
    Add some documentation to this class
    """
    def __init__(self):
        super(DFTB, self).__init__("dftb")
        self.generic_dict_file = 'generic2DFTB.json'

    def prerun(self):
        pass

    def run_job(self, settings, mol, job_name='DFTBjob'):
        """
        Execute an DFTB job with the *ADF* quantum package.

        :param settings: user input settings.
        :type settings: |Settings|
        :param mol: Molecule to run the simulation
        :type mol: Plams Molecule
        :parameter input_file_name: The user can provide a name for the
                                   job input.
        :type input_file_name: String
        :parameter out_file_name: The user can provide a name for the
                                 job output.
        :type out_file_name: String
        :returns: :class:`~qmworks.packages.SCM.DFTB_Result`
        """
        dftb_settings = Settings()
        dftb_settings.input = settings.specific.dftb
        result = plams.DFTBJob(name=job_name, molecule=mol,
                               settings=dftb_settings).run()

        return DFTB_Result(dftb_settings, mol, result._kf)

    def postrun(self):
        pass

    def handle_special_keywords(self, settings, key, value, mol):
        pass


class DFTB_Result(Result):
    """Class providing access to PLAMS DFTBJob result results"""

    def __init__(self, settings, molecule, result):
        self.settings = settings
        self._molecule = molecule
        self.result = result
        properties = 'data/dictionaries/propertiesDFTB.json'
        xs = pkg.resource_string("qmworks", properties)
        self.prop_dict = json2Settings(xs)
        self.properties = self.extract_properties()
        self.archive = files.Path(result.path)

    def as_dict(self):
        return {
            "settings": self.settings,
            "molecule": self._molecule,
            "filename": self.result.path}

    @classmethod
    def from_dict(cls, settings, molecule, filename):
        return DFTB_Result(settings, molecule, plams.kftools.KFFile(filename))

    def extract_properties(self):
        props = Settings()
        for i in range(self.result.read('Properties', 'nEntries')):
            type = self.result.read('Properties', 'Type(' + str(i + 1) + ')').strip()
            subtype = self.result.read('Properties', 'Subtype(' + str(i + 1) + ')').strip()
            value = self.result.read('Properties', 'Value(' + str(i + 1) + ')')
            props[type][subtype] = value
        return props

    def __getattr__(self, prop):
        """Returns a section of the results.

        Example:

        ..

            dipole = result.properties.dipole

        """
        if 'properties' in dir(self) and 'prop_dict' in dir(self):
            return self.properties[self.prop_dict[prop]]
        else:
            raise Exception("NNOETHUNTHN")
        #return '3.23'

    @property
    def molecule(self, unit='bohr', internal=False, n=1):
        """WARNING: Cheap copy from PLAMS, do not keep this!!!"""
        m = self._molecule.copy()
        natoms = len(m)
        coords = self.result.read('Molecule', 'Coords')
        coords = [coords[i: i + 3] for i in range(0, len(coords), 3)]

        if len(coords) > natoms:
            coords = coords[(n - 1) * natoms: n * natoms]
        if internal:
            mapping = self._int2inp()
            coords = [coords[mapping[i] - 1] for i in range(len(coords))]
        for at, coord in zip(m, coords):
            at.move_to(coord, 'bohr')

        return m

adf = ADF()
dftb = DFTB()
