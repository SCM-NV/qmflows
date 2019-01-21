
__all__ = ['find_first_job', 'select_max', 'select_min']

from noodles import schedule, find_first


@schedule
def find_first_job(pred, packagelist, settings, molecule, job_name, **kwargs):
    joblist = [package(
        settings, molecule, job_name=package.pkg_name + "_" + job_name, **kwargs)
        for package in packagelist]
    return find_first(pred, joblist)


@schedule
def select_max(results, prop='energy'):
    """
    Scheduled function to select a result with the maximum value for property from
    a list of results. If the property is not available from a result (e.g. because
    the job failed) the result is ignored.
    :param results:
    :param prop:
    :return:
    """
    filtered_results = [result for result in results if result and result.__getattr__(prop)]
    if len(filtered_results) > 0:
        selected = max(filtered_results, key=lambda item: getattr(item, prop))
        return selected
    else:
        return None


@schedule
def select_min(results, prop='energy'):
    """
    Scheduled function to select a result with the minimum value for property from
    a list of results. If the property is not available from a result (e.g. because
    the job failed) the result is ignored.
    :param results:
    :param prop:
    :return:
    """
    filtered_results = [result for result in results if result and result.__getattr__(prop)]
    if len(filtered_results) > 0:
        selected = min(filtered_results, key=lambda item: getattr(item, prop))
        return selected
    else:
        return None
