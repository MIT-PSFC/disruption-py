from typing import Callable, List

def instantiate_classes(l : List):
    """
    Instantiate all classes in a list of classes and objects.
    
    Parameters
    ----------
    l : List
        List to instantiate classes from.

    Returns
    -------
    List
        The list with all classes instantiated.
    """
    return [x() for x in l if isinstance(x, type)]
    

def without_duplicates(l : List):
    """
    Get list without duplicates maintaining order.

    Parameters
    ----------
    l : List
        List to get without duplicates.

    Returns
    -------
    List
        The list l with duplicates removed.
    """
    seen = set()
    return [x for x in l if not (x in seen or seen.add(x))]