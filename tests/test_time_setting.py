"""
Unit tests for the time_setting

Implemented tests:
- 'ip_efit': tests SharedTimeSetting, IpTimeSetting, and EfitTimeSetting
- SignalTimeSetting

Time settings that are tested elsewhere:
 - EfitTimeSetting: tested in test_quick.py
 - DisruptionTimeSetting: tested in test_against_cache.py
"""

import os

import numpy as np
import pytest

from disruption_py.machine.tokamak import Tokamak
from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.settings.time_setting import SignalTimeSetting, TimeSetting
from disruption_py.workflow import get_shots_data


def run_test_time_setting(
    tokamak: Tokamak,
    time_setting: TimeSetting,
    efit_tree: str,
    shot: int,
    t_start: float,
    t_stop: float,
    length: float,
    test_folder: str,
):
    """
    Retrieve data, then the check time array against the specified targets.
    """
    # Retrieve data
    retrieval_settings = RetrievalSettings(
        efit_nickname_setting=efit_tree,
        run_columns=["kappa_area"],
        time_setting=time_setting,
    )
    shot_data = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=[shot],
        retrieval_settings=retrieval_settings,
        output_setting=os.path.join(test_folder, "output.nc"),
        log_settings=LogSettings(
            console_level="WARNING",
            file_path=os.path.join(test_folder, "output.log"),
        ),
    )
    times = shot_data["time"].to_numpy()
    # Check start, end, and length of time array
    assert times[0] == pytest.approx(t_start, rel=1e-4)
    assert times[-1] == pytest.approx(t_stop, rel=1e-4)
    assert len(times) == length
    # Check for duplicated time point
    assert len(times) == len(np.unique(times))


def test_shared_time_setting(tokamak: Tokamak, test_folder_f: str):
    """
    Test SharedTimeSetting by using the 'ip_efit' shortcut.
    """
    test_setup = {
        Tokamak.CMOD: ["analysis", 1150805012, 0.0601, 1.2799, 6100],
        Tokamak.D3D: ["efit01", 161228, 0.1, 5.0395, 9880],
        Tokamak.EAST: ["efit_east", 55012, 0.301, 5.7, 5401],
    }
    run_test_time_setting(tokamak, "ip_efit", *test_setup[tokamak], test_folder_f)


def test_signal_time_setting(tokamak: Tokamak, test_folder_f: str):
    """
    Test SignalTimeSetting using a signal that is not Ip or a EFIT signal.
    """
    test_setup = {
        Tokamak.CMOD: [
            SignalTimeSetting("spectroscopy", r"\twopi_diode"),
            "analysis",
            1150805012,
            -1.4997,
            3.9559,
            16384,
        ],
        Tokamak.D3D: [
            SignalTimeSetting("rf", r"\top.ech.total:echpwrc"),
            "efit01",
            161228,
            -0.05,
            10,
            201000,
        ],
        Tokamak.EAST: [
            SignalTimeSetting("pcs_east", r"\pcvloop"),
            "efit_east",
            55012,
            -5.5,
            9.199,
            14702,
        ],
    }
    run_test_time_setting(tokamak, *test_setup[tokamak], test_folder_f)
