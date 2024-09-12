#!/usr/bin/env python3
import os

from disruption_py.config import config


def change_directory(test, tmpdir="/tmp"):
    original_dir = os.getcwd()
    os.chdir(tmpdir)

    def wrapper():
        test()
        os.chdir(original_dir)

    return wrapper


@change_directory
def test_settings_file():
    """
    Temporarily change the current working directory to test if the config
    file is reachable.
    """
    assert config().TIME_CONST is not None
