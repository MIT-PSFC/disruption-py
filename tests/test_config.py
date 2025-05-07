#!/usr/bin/env python3

"""This module contains tests to ensure all of the config settings load properly."""

import os

from disruption_py.config import config


def change_directory(test, tmpdir="/tmp"):
    """
    Change the current working directory before a test and revert back to the
    original directory after the test completes.
    """

    def wrapper():
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        test()
        os.chdir(original_dir)

    return wrapper


@change_directory
def test_settings_file():
    """
    Temporarily change the current working directory to test if the config settings
    file is reachable.
    """
    assert config().time.time_const is not None


def test_access_tokamak_settings(tokamak):
    """
    Test each tokamak's unique settings are accessible.
    """
    assert config(tokamak)
