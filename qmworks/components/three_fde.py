__all__ = ['mfcc', 'MFCC_Result']

from qmworks import templates, rdkitTools
from qmworks.settings import Settings
from qmworks.packages import Result, run
from noodles import gather, schedule

@schedule
def mfcc(package, frags, caps, settings):
    mfcc_result = MFCC_Result()
    frag_jobs = [package(templates.singlepoint, rdkitTools.rdkit2plams(frag), job_name="mfcc_frag_" + str(i))
                 for i, frag in enumerate(frags)]
    frag_results = gather(*frag_jobs)
    cap_jobs = [package(templates.singlepoint, rdkitTools.rdkit2plams(cap), job_name="mfcc_cap_" + str(i))
                for i, cap in enumerate(caps)]
    cap_results = gather(*cap_jobs)
    mfcc_result.results = run(gather(frag_results, cap_results))
    return mfcc_result

class MFCC_Result(Result):
    def __init__(self):
        self.results = [[], []]

    def get_dipole_vector(self):
        # pylint: disable-msg=E1101
        import numpy
        dipole = numpy.zeros(3)
        # loop over fragments
        for frag_result in self.results[0]:
            dipole += frag_result.dipole
        for cap_result in self.results[1]:
            dipole -= cap_result.dipole
        return dipole

