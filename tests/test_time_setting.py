"""
Unit tests for the time_setting

- CMOD: 'ip_efit', SignalTimeSetting
- D3D: 'disruption', 'ip_efit', SignalTimeSetting
- EAST: 'disruption', 'ip_efit', SignalTimeSetting

(EfitTimeSetting is tested in test_quick.py)
('ip_efit' tests both IpTimeSetting and MixedTimeSetting)
"""

import numpy as np
import pytest

from disruption_py.machine.tokamak import Tokamak
from disruption_py.settings import RetrievalSettings
from disruption_py.settings.time_setting import SignalTimeSetting
from disruption_py.workflow import get_shots_data


def run_test_time_setting(time_setting, efit_tree, shotno, t_start, t_stop, length):
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
        shotlist_setting=[shotno],
        retrieval_settings=retrieval_settings,
    )
    times = shot_data["time"].to_numpy()
    # Check start, end, and length of time array
    assert times[0] == pytest.approx(t_start, rel=1e-4)
    assert times[-1] == pytest.approx(t_stop, rel=1e-4)
    assert len(times) == length
    # Check for duplicated time point
    assert len(times) == len(np.unique(times))


def test_shared_time_setting(tokamak: Tokamak):
    """
    Test SharedTimeSetting by using the 'ip_efit' shortcut.
    """
    test_setup = {
        Tokamak.CMOD: ["analysis", 1150805012, 0.0601, 1.2799, 6100],
        Tokamak.D3D: ["efit01", 161228, 0.1, 5.0395, 9880],
        Tokamak.EAST: ["efit_east", 55012, 0.301, 5.7, 5401],
    }
    run_test_time_setting("ip_efit", *test_setup[tokamak])


def test_signal_time_setting(tokamak: Tokamak):
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
    run_test_time_setting(*test_setup[tokamak])


def test_disruption_time_setting(tokamak: Tokamak):
    """
    Test DisruptionTimeSetting for DIII-D and EAST, skip for C-MOD.
    """
    test_setup = {
        Tokamak.D3D: [161228, 0.1, 5.0935, 247],
        Tokamak.EAST: [55012, 0.2, 5.7113, 79],
    }
    if tokamak not in test_setup:
        pytest.skip(f"DisruptionTimeSetting is not implemented for {tokamak.name}.")
    run_test_time_setting("disruption", "disruption_warning", *test_setup)
