
__all__ = ['find_first_job', 'select_max', 'select_min', 'select_max_dict', 'select_min_dict']

from noodles import schedule_hint, find_first


@schedule_hint()
def find_first_job(pred, packagelist, settings, molecule, job_name, **kwargs):
    joblist = [package(
        settings, molecule, job_name=package.pkg_name + "_" + job_name, **kwargs)
        for package in packagelist]
    return find_first(pred, joblist)


@schedule_hint()
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
        selected = max(filtered_results, key=lambda item: item.__getattr__(prop))
        print("Selected {:s}: {:s}".format(selected.job_name, str(selected.__getattr__(prop))))
        return selected
    else:
        return None


@schedule_hint()
def select_max_dict(dictionary, prop):
    def key(k):
        ret = getattr(dictionary[k], prop)
        if isinstance(ret, (float, int)):
            return ret
        else:
            return float("-inf")
    selected_key = max(dictionary, key=key)
    print("Selected " + str(selected_key) + ": " + str(getattr(dictionary[selected_key], prop)))
    return dictionary[selected_key]

@schedule_hint()
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
        selected = min(filtered_results, key=lambda item: item.__getattr__(prop))
        print("Selected {:s}: {:s}".format(selected.job_name, str(selected.__getattr__(prop))))
        return selected
    else:
        return None

@schedule_hint()
def select_min_dict(dictionary, prop):
    def key(k):
        ret = getattr(dictionary[k], prop)
        if isinstance(ret, (float, int)):
            return ret
        else:
            return float("inf")
    selected_key = min(dictionary, key=key)
    print("Selected " + str(selected_key) + ": " + str(getattr(dictionary[selected_key], prop)))
    return dictionary[selected_key]

