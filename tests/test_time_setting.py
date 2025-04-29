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


def run_test_time_setting(efit_tree, time_setting, shotno, t_start, t_stop, length):
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
    if tokamak is Tokamak.CMOD:
        test_setup = ["analysis", "ip_efit", 1150805012, 0.0601, 1.2799, 6100]
    elif tokamak is Tokamak.D3D:
        test_setup = ["efit01", "ip_efit", 161228, 0.1, 5.0395, 9880]
    elif tokamak is Tokamak.EAST:
        test_setup = ["efit_east", "ip_efit", 55012, 0.301, 5.7, 5401]
    else:
        raise ValueError(f"Unsupported tokamak: {tokamak}.")

    run_test_time_setting(*test_setup)


def test_signal_time_setting(tokamak: Tokamak):
    """
    Test SignalTimeSetting using a signal that is not Ip or a EFIT signal.
    """
    if tokamak is Tokamak.CMOD:
        time_setting = SignalTimeSetting("spectroscopy", r"\twopi_diode")
        test_setup = ["analysis", time_setting, 1150805012, -1.4997, 3.9559, 16384]
    elif tokamak is Tokamak.D3D:
        time_setting = SignalTimeSetting("rf", r"\top.ech.total:echpwrc")
        test_setup = ["efit01", time_setting, 161228, -0.05, 10, 201000]
    elif tokamak is Tokamak.EAST:
        time_setting = SignalTimeSetting("pcs_east", r"\pcvloop")
        test_setup = ["efit_east", time_setting, 55012, -5.5, 9.199, 14702]
    else:
        raise ValueError(f"Unsupported tokamak: {tokamak}.")

    run_test_time_setting(*test_setup)


def test_disruption_time_setting(tokamak: Tokamak):
    """
    Test DisruptionTimeSetting for D3D and EAST.

    DisruptionTimeSetting is not implemented for CMOD and therefore we skip this test.
    """
    if tokamak is Tokamak.CMOD:
        return
    if tokamak is Tokamak.D3D:
        test_setup = ["disruption_warning", "disruption", 161228, 0.1, 5.0935, 247]
    elif tokamak is Tokamak.EAST:
        test_setup = ["disruption_warning", "disruption", 55012, 0.2, 5.7113, 79]
    else:
        raise ValueError(f"Unsupported tokamak: {tokamak}.")

    run_test_time_setting(*test_setup)
