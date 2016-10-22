__all__ = ['mfcc', 'MFCC_Result']

from qmworks import templates, molkit
from qmworks.settings import Settings
from qmworks.packages import Result, run
from noodles import gather, schedule

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


