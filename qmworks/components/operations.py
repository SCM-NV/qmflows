
__all__ = ['select_max']

from noodles import schedule_hint


@schedule_hint()
def select_max(results, property='energy'):
    print([a.energy for a in results])
    selected_result = max(results, key=lambda item: item.__getattr__(property))
    print("max " + property + ":", selected_result.__getattr__(property))
    return selected_result
