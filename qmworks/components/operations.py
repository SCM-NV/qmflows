
__all__ = ['select_max', 'select_min']

from noodles import schedule_hint


@schedule_hint()
def select_max(results, prop='energy'):
    """
    Scheduled function to select a result with the maximum value for property
    from a list or list of lists
    :param results:
    :param prop:
    :return:
    """
    max_res = sel_max(results, prop)
    print("Selected " + str(max_res.job_name) + ": " + str(max_res.__getattr__(prop)))
    return max_res


def sel_max(results, prop):
    line = "From " + prop + " values: "
    for i in range(len(results)):
        if isinstance(results[i], list):
            n = sel_max(results[i], prop)
            results[i] = n
        else:
            line += "{:12.6f}".format(results[i].__getattr__(prop), end="")
    print(line)
    selected_result = max(results, key=lambda item: item.__getattr__(prop))
    return selected_result


@schedule_hint()
def select_min(results, prop='energy'):
    """
    Scheduled function to select a result with the minimum value for property from
    a list or list of lists
    :param results:
    :param prop:
    :return:
    """
    min_res = sel_min(results, prop)
    print("Selected " + str(max_res.job_name) + ": "+ str(min_res.__getattr__(prop)))
    return min_res


def sel_min(results, prop):
    line = "From " + prop + " values: "
    for i in range(len(results)):
        if isinstance(results[i], list):
            n = sel_min(results[i], prop)
            results[i] = n
        else:
            line += "{:12.6f}".format(results[i].__getattr__(prop), end="")
    print(line)
    selected_result = min(results, key=lambda item: item.__getattr__(prop))
    return selected_result
