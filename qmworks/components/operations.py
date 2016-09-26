
__all__ = ['select_max', 'select_min']

from noodles import schedule_hint
from qmworks.packages import Result


@schedule_hint()
def select_max(results, property='energy'):
    """
    Scheduled function to select a result with the maximum value for property
    from a list or list of lists
    :param results:
    :param property:
    :return:
    """
    max_res = sel_max(results, property)
    print("Appr TS energy " + str(max_res.job_name) + ": " + str(max_res.energy))
    return max_res

def sel_max(results, property):
    line = ""
    for i in range(len(results)):
        if isinstance(results[i], list):
            n = sel_max(results[i], property)
            results[i] = n
        else:
            line += "{:12.6f}".format(results[i].energy, end="")
    print(line)
    selected_result = max(results, key=lambda item: item.__getattr__(property))
    return selected_result

@schedule_hint()
def select_min(results, property='energy'):
    """
    Scheduled function to select a result with the minimum value for property from
    a list or list of lists
    :param results:
    :param property:
    :return:
    """
    min_res = sel_min(results, property)
    print("Appr TS energy: " + str(min_res.energy))
    return min_res

def sel_min(results, property):
    line = ""
    for i in range(len(results)):
        if isinstance(results[i], list):
            n = sel_min(results[i], property)
            results[i] = n
        else:
            line += "{:12.6f}".format(results[i].energy, end="")
    print(line)
    selected_result = min(results, key=lambda item: item.__getattr__(property))
    return selected_result