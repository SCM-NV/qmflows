__all__ = ['mfcc', 'MFCC_Result', 'Fragment', 'adf_fragmentsjob']

from qmworks import Settings, templates
from qmworks.packages import Result
from qmworks.packages.SCM import adf
from noodles import gather, schedule, Storable
from plams import Molecule

import numpy as np

@schedule
class MFCC_Result(Result):
    def __init__(self, frags, caps):
        self.frags = frags
        self.caps = caps

    @property
    def dipole(self):
        dipole = np.zeros(3)
        # loop over fragments
        for frag_result in self.frags:
            dipole += frag_result.dipole
        # loop over caps
        for cap_result in self.caps:
            dipole -= cap_result.dipole
        return dipole


def mfcc(package, frags, caps, settings=None):
    mfcc_settings = templates.singlepoint
    if settings:
        mfcc_settings = mfcc_settings.overlay(settings)
    frag_jobs = [package(mfcc_settings, frag, job_name="mfcc_frag_" + str(i))
                 for i, frag in enumerate(frags)]
    cap_jobs = [package(mfcc_settings, cap, job_name="mfcc_cap_" + str(i))
                for i, cap in enumerate(caps)]
    return MFCC_Result(gather(*frag_jobs), gather(*cap_jobs))


class Fragment(Storable):
    def __init__(self, result, mol_list):
        self.result = result
        self.mol_list = gather(*mol_list)

@schedule
def adf_fragmentsjob(settings, mol, *frozen_frags):
    mol_tot = Molecule()
    frag_settings = Settings()
    for i, frag in enumerate(frozen_frags):
        frag_id = 'frag' + str(i + 1)
        for m in frag.mol_list:
            for a in m:
                a.fragment = frag_id
            mol_tot += m
        path = frag.result.kf.path + ' type=FDE'
        frag_settings.specific.adf.fragments[frag_id] = path
        frag_settings.specific.adf.fde.PW91k = ""
    mol_tot += mol
    return adf(settings.overlay(frag_settings), mol_tot, job_name="fde")
