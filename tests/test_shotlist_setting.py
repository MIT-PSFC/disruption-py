#!/usr/bin/env python3

"""Unit tests for the shotlist_setting"""

import numpy as np
import pytest

from disruption_py.settings.shotlist_setting import (
    ShotlistSetting,
    ShotlistSettingParams,
    shotlist_setting_runner,
)

# Use a few integers instead of reasonable shot ids because there's nothing
# grabbing actual shot data here
REFERENCE_SHOTLIST = list(range(5))


def get_shotlists():
    """
    Generate custom shot requests, along with 1, 2, and 3 dimensional Python lists
    and Numpy arrays.
    """

    class CustomShotlistSetting(ShotlistSetting):
        """A custom shotlist setting for testing"""

        def _get_shotlist(self, params: ShotlistSettingParams):
            """Return Numpy array of shots"""
            return np.array([REFERENCE_SHOTLIST])

    return [
        CustomShotlistSetting(),
        np.array(REFERENCE_SHOTLIST),
        np.array([REFERENCE_SHOTLIST]),
        np.array([[[i] for i in REFERENCE_SHOTLIST]]),
        REFERENCE_SHOTLIST,
        [REFERENCE_SHOTLIST],
        [[[i] for i in REFERENCE_SHOTLIST]],
    ]


@pytest.mark.parametrize("shotlist", get_shotlists())
def test_shotlist_setting_runner(shotlist):
    """
    Ensure all variations of shotlist_settings get flattened to a single dimensional
    list.
    """
    shot_ids_request_params = ShotlistSettingParams(database=None, tokamak=None)
    result = shotlist_setting_runner(shotlist, shot_ids_request_params)
    assert REFERENCE_SHOTLIST == result
