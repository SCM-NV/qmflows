__all__ = ['mfcc', 'MFCC_Result', 'adf3fde', 'ADF3FDE_Result', 'Fragment', 'adf_fragmentsjob']

from qmflows import templates
from qmflows.settings import Settings
from qmflows.packages import Result
from qmflows.packages.SCM import adf
from noodles import gather, schedule, Storable
from scm.plams import Molecule

import numpy as np


@schedule
class MFCC_Result(Result):
    def __init__(self, frags, caps):
        self.frags = [Fragment(f, f.molecule) for f in frags]
        self.caps = [Fragment(c, c.molecule) for c in caps]

    @property
    def dipole(self):
        dipole = np.zeros(3)
        # loop over fragments
        for frag in self.frags:
            dipole += frag.result.dipole
        # loop over caps
        for cap in self.caps:
            dipole -= cap.result.dipole
        return dipole


def mfcc(package, frags, caps, settings=None):
    mfcc_settings = templates.singlepoint
    if settings:
        mfcc_settings = mfcc_settings.overlay(settings)

    capped_atoms = {}
    for i, cap in enumerate(caps):
        for a in cap:
            capped_atoms[a.coords] = a

    for i, frag in enumerate(frags):
        cap_label = 'cap' + str(i + 1)
        for a in frag:
            if a.coords in capped_atoms:
                a.properties.name = cap_label
                capped_atoms[a.coords].properties.name = cap_label

    frag_jobs = [package(mfcc_settings, frag, job_name="mfcc_frag_" + str(i))
                 for i, frag in enumerate(frags)]
    cap_jobs = [package(mfcc_settings, cap, job_name="mfcc_cap_" + str(i))
                for i, cap in enumerate(caps)]
    return MFCC_Result(gather(*frag_jobs), gather(*cap_jobs))


class Fragment(Storable):
    def __init__(self, result, mol, isfrozen=True, pack_tape=False):
        self.result = result
        if pack_tape:
            pack_t21(self.result.kf.path)
        #self.mol_list = gather(*mol_list)
        self.mol = mol
        self.isfrozen = isfrozen


@schedule
def adf_fragmentsjob(settings, frags, caps=None, fragment_settings=None, job_name='fde'):
    mol_tot = Molecule()
    frag_settings = Settings()
    cap_ids = {}
    if caps:
        for i, cap in enumerate(caps):
            cap_id = 'cap' + str(i + 1)
            for a in cap.mol:
                cap_ids[a.coords] = cap_id
            path = cap.result.kf.path
            key = cap_id + ' ' + path + ' type=FDEsubstract &'
            frag_settings.specific.adf.fragments[key] = fragment_settings
    for i, frag in enumerate(frags):
        frag_id = 'frag' + str(i + 1)
        if frag.result:
            for a in frag.mol:
                a.properties.adf.fragment = frag_id
                if a.coords in cap_ids:
                    a.properties.adf.fragment += '  fs=' + cap_ids[a.coords]
            path = frag.result.kf.path
            key = frag_id + ' ' + path + ' subfrag=active'
            if frag.isfrozen:
                key += ' type=FDE'
                if fragment_settings:
                    key += ' &'
                    frag_settings.specific.adf.fragments[key] = fragment_settings
                else:
                    frag_settings.specific.adf.fragments[key] = ""
            else:
                frag_settings.specific.adf.fragments[key] = ""
        mol_tot += frag.mol
    frag_settings.specific.adf.fde.PW91k = ""
    return adf(settings.overlay(frag_settings), mol_tot, job_name=job_name)


class ADF3FDE_Result(Result):
    def __init__(self, frags, caps):
        self.frags = frags
        self.caps = caps

    @property
    def dipole(self):
        dipole = np.zeros(3)
        # loop over fragments
        for frag in self.frags:
            dipole += frag.result.dipole
        # loop over caps
        for cap in self.caps:
            dipole -= cap.result.dipole
        return dipole


def adf3fde(frags, caps, settings, fde_settings, fragment_settings, cycles=1):
    adf3fde_settings = templates.singlepoint
    adf3fde_settings.specific.adf.allow = "partialsuperfrags"
    if settings:
        adf3fde_settings = settings.overlay(adf3fde_settings)
    adf3fde_settings.specific.adf.fde = fde_settings

    for i in range(cycles):
        frags = adf3fde_cycle(
            frags, caps, adf3fde_settings, fragment_settings, job_name='fde_' + str(i))
    return schedule(ADF3FDE_Result)(frags, caps)


@schedule
def adf3fde_cycle(frags, caps, adf3fde_settings, fragment_settings, job_name='fde'):
    new_frags = []
    for i, frag in enumerate(frags):
        frag.isfrozen = False
        new_frags.append(schedule(Fragment)(adf_fragmentsjob(
            adf3fde_settings, frags, caps, fragment_settings,
            job_name=job_name + '_' + str(i)), frag.mol, pack_tape=True))
        frag.isfrozen = True

    return gather(*new_frags)


def pack_t21(fn):
    """
    Pack a TAPE21 results file of an FDE job.

    This will reduce the size of a FDE TAPE21 file by
    only keeping the ActiceFragment section and related
    essential information. All information about frozen
    fragments will be deleted. This is useful for manual
    freeze-and-thaw calculations with many subsystems,
    where this packing can reduce the size of the stored
    files significantly.

    @param fn: The filename of the TAPE21 to be packed.
    @type  fn: str
    """
    from scm.plams import KFReader
    import os
    import subprocess

    fn_orig = fn + ".orig"
    os.rename(fn, fn_orig)
    kf_orig = KFReader(fn_orig)

    keepsections = ['General', 'Properties', 'Num Int Params', 'LinearScaling', 'Geometry']
    atomtypeindices = kf_orig.read('ActiveFrag', 'atomtypeIndices')
    if not kf_orig._sections:
        kf_orig._create_index()
    for section in kf_orig._sections:
        for i in atomtypeindices:
            if section.startswith('Atyp%3i' % i):
                keepsections.append(section)
            elif section.startswith('Atyp%4i' % i):
                keepsections.append(section)
            elif section.startswith('Atyp%5i' % i):
                keepsections.append(section)
        if section.startswith('ActiveFrag'):
            keepsections.append(section)
    subprocess.Popen(['cpkf', fn_orig, fn] + keepsections,
                     stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
    os.remove(fn_orig)
