
# =======>  Standard and third party Python Libraries <======

from noodles import files
from os.path import join
from qmworks.settings import Settings
from qmworks.packages.packages import Package, Result  # ChemResult

import plams
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
        path_t21 = result._kf.path

        return ADF_Result(adf_settings, mol, job_name, path_t21,
                          plams_dir=result.job.path)

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
                        name = 'atom ' + str(a + 1)
                        settings.specific.adf.constraints[name] = ""
        else:
            raise RuntimeError('Keyword ' + key + ' doesn\'t exist')


class ADF_Result(Result):
    """Class providing access to PLAMS ADFJob result results"""

    def __init__(self, settings, molecule, job_name, path_t21, plams_dir=None):
        properties = 'data/dictionaries/propertiesADF.json'
        super().__init__(settings, molecule, job_name, plams_dir=plams_dir,
                         properties=properties)
        self.result = plams.kftools.KFFile(path_t21)

    def as_dict(self):
        """
        Method to serialize as a JSON dictionary the results given
        by an ADF computation.
        """
        return {
            "settings": self.settings,
            "molecule": self._molecule,
            "job_name": self.job_name,
            "archive": self.archive}

    @classmethod
    def from_dict(cls, settings, molecule, job_name, archive):
        """
        Methods to deserialize an `ADF_Result` object.
        """
        plams_dir = archive["plams_dir"]
        path_t21 = join(plams_dir, '{}.t21'.format(job_name))
        return ADF_Result(settings, molecule, path_t21, job_name, plams_dir)

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

        return DFTB_Result(dftb_settings, mol, result.job.path, result.job.name)

    def postrun(self):
        pass

    def handle_special_keywords(self, settings, key, value, mol):
        pass


class DFTB_Result(Result):
    """Class providing access to PLAMS DFTBJob result results"""

    def __init__(self, settings, molecule, path, job_name):
        properties = 'data/dictionaries/propertiesDFTB.json'
        super().__init__(settings, molecule, job_name, properties=properties)
        self.path = path
        kf_filename = self.path + '/' + job_name + '.rkf'
        self.kf = plams.kftools.KFFile(kf_filename)
        self.properties = self.extract_properties()
        self.archive = files.Path(self.kf)

    def as_dict(self):
        return {
            "settings": self.settings,
            "molecule": self._molecule,
            "path": self.path,
            "job_name": self.job_name}

    @classmethod
    def from_dict(cls, settings, molecule, path, name):
        return DFTB_Result(settings, molecule, path, name)

    def extract_properties(self):
        props = Settings()
        for i in range(self.kf.read('Properties', 'nEntries')):
            type = self.kf.read('Properties', 'Type(' + str(i + 1) + ')').strip()
            subtype = self.kf.read('Properties', 'Subtype(' + str(i + 1) + ')').strip()
            value = self.kf.read('Properties', 'Value(' + str(i + 1) + ')')
            props[type][subtype] = value
        return props

    def __getattr__(self, prop):
        """Returns a section of the results.

        Example:

        ..

            dipole = result.properties.dipole

        """
        if 'properties' in dir(self) and prop in self.prop_dict:
            prop_query = self.prop_dict[prop]
            if isinstance(prop_query, str):
                if prop_query[:4] == 'awk|':
                    return self.awk_output(script=prop_query[4:])
                return self.properties[prop_query]
            else:
                print(prop_query)
                return self.kf.read(*prop_query)
        else:
            raise Exception("NNOETHUNTHN")
        #return '3.23'

    @property
    def molecule(self, unit='bohr', internal=False, n=1):
        """WARNING: Cheap copy from PLAMS, do not keep this!!!"""
        m = self._molecule.copy()
        natoms = len(m)
        coords = self.kf.read('Molecule', 'Coords')
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
