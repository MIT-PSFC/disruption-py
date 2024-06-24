#!/usr/bin/env python3

import numpy as np

from disruption_py.settings.shotlist_setting import ShotlistSetting, ShotlistSettingParams, shotlist_setting_runner


def test_shotlist_setting_runner():
    """
    Test custom shotlist, along with 1, 2, and 3 dimensional Python lists
    and Numpy arrays to ensure they all get flattened to a single dimensional
    list.
    """

    class CustomShotlistSetting(ShotlistSetting):
        def _get_shot_ids(self, params):
            return np.array([[1160405002, 1140523021, 1140523026, 1160620011]])

    shot_lists = [
        CustomShotlistSetting(),
        np.array([1160405002, 1140523021, 1140523026, 1160620011]),
        np.array([[1160405002, 1140523021, 1140523026, 1160620011]]),
        np.array([[[1160405002], [1140523021], [1140523026], [1160620011]]]),
        [1160405002, 1140523021, 1140523026, 1160620011],
        [[1160405002, 1140523021, 1140523026, 1160620011]],
        [[[1160405002], [1140523021], [1140523026], [1160620011]]],
    ]
    expected = [1160405002, 1140523021, 1140523026, 1160620011]

    shot_ids_request_params = ShotlistSettingParams(
        database=None, tokamak=None, logger=None
    )
    for shot_ids_request in shot_lists:
        result = shotlist_setting_runner(shot_ids_request, shot_ids_request_params)
        assert expected == result