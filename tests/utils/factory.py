#!/usr/bin/env python3

"""This module provides utility functions for retrieving testing settings like the
shots and columns based on the provided tokamak."""

import os

from disruption_py.config import config
from disruption_py.machine.tokamak import Tokamak


def get_tokamak_test_expected_failure_columns(tokamak: Tokamak):
    """Return the columns expected to fail from the config"""
    return config(tokamak).testing.EXPECTED_FAILURE_COLUMNS


def get_tokamak_test_shotlist(tokamak: Tokamak) -> list[int]:
    """Return the shot ids used for testing from the config. Return a smaller
    shotlist when running as a Github action."""
    shot_id_dict = config(tokamak).testing.TEST_SHOTS

    if "GITHUB_ACTIONS" in os.environ:
        shot_id_dict = {
            key: value for key, value in shot_id_dict.items() if "_fast" in key
        }

    return list(shot_id_dict.values())


def get_tokamak_test_columns(tokamak: Tokamak):
    """Return the columns used for testing from the config."""
    return config(tokamak).testing.TEST_COLUMNS
