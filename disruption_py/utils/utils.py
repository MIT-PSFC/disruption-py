def without_duplicates(l):
    """
    Get list without duplicates maintaining order.

    Parameters
    ----------
    l : list
        List to get without duplicates.

    Returns
    -------
    _ : list
        The list l with duplicates removed.
    """
    seen = set()
    return [x for x in l if not (x in seen or seen.add(x))]