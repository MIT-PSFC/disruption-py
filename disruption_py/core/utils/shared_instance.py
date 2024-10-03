#!/usr/bin/env python3

"""
This module provides a singleton class that ensures only one instance of a 
given class is created for a specific set of arguments, allowing for 
shared instances across the same process.
"""

import os


class SharedInstance:
    """
    Singleton class for creating shared instances of a specified class.

    Attributes
    ----------
    cls_arg : type
        The class for which shared instances will be created.
    _instances : dict
        A dictionary to store shared instances keyed by a unique identifier.
    """

    _instances = {}

    def __init__(self, cls_arg):
        """
        Initialize with a class.

        Parameters
        ----------
        cls_arg : type
            The class for which instances will be shared.
        """
        self.cls_arg = cls_arg

    def get_instance(self, *args, **kwargs):
        """
        Get a shared instance of the specified class, creating it if necessary.

        Parameters
        ----------
        *args : tuple
            Positional arguments to pass to the class constructor.
        **kwargs : dict
            Keyword arguments to pass to the class constructor.

        Returns
        -------
        object
            The shared instance of the specified class.
        """

        # Convert any dictionary in args or kwargs to a hashable form
        def make_hashable(obj):
            if isinstance(obj, dict):
                return tuple(sorted((k, make_hashable(v)) for k, v in obj.items()))
            if isinstance(obj, (list, set)):
                return tuple(sorted(make_hashable(e) for e in obj))
            return obj

        pid = os.getpid()
        hashable_args = tuple(make_hashable(arg) for arg in args)
        hashable_kwargs = tuple(
            sorted((k, make_hashable(v)) for k, v in kwargs.items())
        )

        key = (pid, self.cls_arg, hashable_args, hashable_kwargs)
        if key not in SharedInstance._instances:
            instance = self.cls_arg(*args, **kwargs)
            SharedInstance._instances[key] = instance
        return SharedInstance._instances[key]
