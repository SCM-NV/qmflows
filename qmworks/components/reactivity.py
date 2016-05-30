
__all__ = ['PES_scan']

from qmworks import templates
from qmworks.settings import Settings
from noodles import gather


def PES_scan(package, settings, molecule, scan, job_name="PESscan"):
    """
    This function enables multidimensional potential energy surface scans.
    The argument 'scan' should be a dictionary with keys 'constraint' and 'surface'
    and optionally 'scan'. The latter allows nested (multidimensional) scans.
    
    Example
    .. code-block:: python

        scan = {'constraint': "dist 1 5",
                'surface': {'nsteps':6, 'start': 2.3, 'stepsize': 0.1},
                'scan': {'constraint': "dist 3 4",
                         'surface': {'nsteps':6, 'start': 2.3, 'stepsize': 0.1}
                         }
                }   
    
    It is also possible to scan over two coordinates in a concerted way, for example:
    .. code-block:: python

        scan = {'constraint': ["dist 1 5", "dist 3 4"],
                   'surface': {'nsteps':6, 'start':[2.3,2.3], 'stepsize': [0.1,0.1]} }
    
    In this example multiple listed values for 'start' and 'stepsize' correspond to
    the same number of listed 'constraint's
    
    :returns:
        A list of promised result objects
    
    """   
    lt_jobs=[]

    def lt_job(settings, job_name):
        if isinstance(package, list):
            optimized_mol = package[0](templates.geometry.overlay(settings),molecule, job_name=job_name+"_opt").molecule
            result = package[1](templates.singlepoint.overlay(settings),optimized_mol, job_name=job_name)
        else:
            result = package(templates.geometry.overlay(settings),molecule, job_name=job_name)
        return result

    def get_constraint_settings(constraint, start, stepsize, step):
        s = Settings()
        if isinstance(constraint, list):
            for c in range(len(constraint)):
                s.constraint[constraint[c]] = start[c] + stepsize[c] * step
        else:
            s.constraint[constraint] = start + stepsize * step
        return s 
        
    def perform_scan(scan, settings, job_name):
        for step in range(scan['surface']['nsteps']):
            constraint = scan['constraint']
            start = scan['surface']['start']
            stepsize = scan['surface']['stepsize']
            new_settings = settings.overlay(get_constraint_settings(constraint, start, stepsize, step))
            if 'scan' in scan:
                perform_scan(scan['scan'],new_settings, job_name+"_"+str(step))
            else:
                lt_jobs.append(lt_job(new_settings, job_name+"_"+str(step)))

    perform_scan(scan,settings,job_name)
    return gather(*lt_jobs)
    #return max(lt_jobs, key = lambda item: item.energy)
