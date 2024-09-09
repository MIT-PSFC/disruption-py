#!/usr/bin/env python3

import os


class SharedInstanceFactory:
    _instances = {}

    def __init__(self, cls_arg):
        self.cls_arg = cls_arg

    def get_instance(self, *args, **kwargs):
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
        if key not in SharedInstanceFactory._instances:
            instance = self.cls_arg(*args, **kwargs)
            SharedInstanceFactory._instances[key] = instance
        return SharedInstanceFactory._instances[key]
