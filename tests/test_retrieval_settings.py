#!/usr/bin/env python3

"""
Unit tests for the retrieval settings which covers data retrieval from various
sources and the time domain of the data retrieved.
"""

import os

import pytest
import xarray as xr

from disruption_py.machine.tokamak import Tokamak
from disruption_py.settings import RetrievalSettings
from disruption_py.settings.log_settings import LogSettings
from disruption_py.workflow import get_shots_data
from tests.conftest import skip_on_fast_execution


@pytest.fixture(scope="module", name="full_domain_data")
def full_domain_data_fixture(tokamak, shotlist, test_folder_m) -> xr.Dataset:
    """
    Get data for the full time domain.
    """
    return get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        output_setting=os.path.join(test_folder_m, "output.nc"),
        log_settings=LogSettings(
            console_level="WARNING",
            log_file_path=os.path.join(test_folder_m, "output.log"),
        ),
        num_processes=2,
    )


@skip_on_fast_execution
def test_only_requested_columns(tokamak, shotlist, test_folder_f):
    """
    Ensure `only_requested_columns` works. `ip` is returned by
    `get_ip_parameters`, so we should not see any of the other quantities like
    `dip_dt` returned. `q95` is from efit, so none of the other efit quantities
    should be returned.
    """
    run_columns = ["ip", "q95"]
    retrieval_settings = RetrievalSettings(
        run_columns=run_columns,
        only_requested_columns=True,
    )
    results = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        output_setting=os.path.join(test_folder_f, "output.nc"),
        log_settings=LogSettings(
            console_level="WARNING",
            log_file_path=os.path.join(test_folder_f, "output.log"),
        ),
        num_processes=2,
    )
    assert set(run_columns) == set(results.data_vars)


@skip_on_fast_execution
@pytest.mark.parametrize("domain_setting", ["flattop", "rampup_and_flattop"])
def test_domain_setting(
    tokamak, shotlist, domain_setting, full_domain_data, test_folder_f
):
    """
    Test the two partial domain settings by comparing their start and end times
    with the full domain.
    """
    if (
        tokamak in [Tokamak.D3D, Tokamak.EAST]
        and domain_setting == "rampup_and_flattop"
    ):
        pytest.skip(f"{domain_setting} domain setting not defined for {tokamak.value}")
    retrieval_settings = RetrievalSettings(
        efit_nickname_setting="default", domain_setting=domain_setting
    )
    part_domain_data = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        output_setting=os.path.join(test_folder_f, "output.nc"),
        log_settings=LogSettings(
            console_level="WARNING",
            log_file_path=os.path.join(test_folder_f, "output.log"),
        ),
        num_processes=2,
    )
    assert (
        set(shotlist)
        == set(part_domain_data.shot.values)
        == set(full_domain_data.shot.values)
    )
    for shot in shotlist:
        part_domain = part_domain_data.sel(idx=part_domain_data.shot == shot)
        full_domain = full_domain_data.sel(idx=full_domain_data.shot == shot)
        p_start, p_end = part_domain.time.values[0], part_domain.time.values[-1]
        f_start, f_end = full_domain.time.values[0], full_domain.time.values[-1]
        if domain_setting == "flattop":
            # Use <= because a shot may end during the flattop,
            assert f_start < p_start < p_end <= f_end
        else:
            assert f_start == p_start < p_end < f_end


@skip_on_fast_execution
@pytest.mark.parametrize(
    "run_methods, run_columns, expected_cols, forbidden_cols",
    [
        # Test run_methods with run_columns=None
        (None, None, None, []),
        ([], None, [], []),
        (["~get_kappa_area"], None, None, ["kappa_area"]),
        (["get_kappa_area"], None, {"kappa_area"}, []),
        # Test run_columns with run_methods=None
        (None, [], [], []),
        (None, ["kappa_area"], {"kappa_area"}, []),
        # Test run_methods and run_columns combo
        ([], [], [], []),
        (["get_kappa_area"], [], {"kappa_area"}, []),
        (
            ["get_kappa_area"],
            ["greenwald_fraction"],
            {"kappa_area", "n_e", "dn_dt", "greenwald_fraction"},
            [],
        ),
        (
            ["~get_kappa_area"],
            ["greenwald_fraction"],
            {"n_e", "dn_dt", "greenwald_fraction"},
            [],
        ),
        (["~get_kappa_area"], ["kappa_area"], [], []),
    ],
)
def test_run_methods_and_columns(
    tokamak,
    shotlist,
    run_methods,
    run_columns,
    expected_cols,
    forbidden_cols,
    full_domain_data,
    test_folder_f,
):
    """
    Test the `run_methods` and `run_columns` parameters of RetrievalSettings.

    - If both are None, all methods are run
    - If one is None or [], the methods/columns specified by the other are run
    - If both are specified, the combined methods/columns are run
    - If `run_methods` excludes a method returning a column specified in `run_columns`,
      the method is not run
    """
    retrieval_settings = RetrievalSettings(
        run_methods=run_methods,
        run_columns=run_columns,
    )
    results = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        output_setting=os.path.join(test_folder_f, "output.nc"),
        log_settings=LogSettings(
            console_level="CRITICAL",
            log_file_path=os.path.join(test_folder_f, "output.log"),
        ),
        num_processes=2,
    )
    # Expected columns None means all columns (except forbidden cols) are returned
    if expected_cols is None:
        assert set(results.data_vars) == set(
            k for k in full_domain_data.data_vars if k not in forbidden_cols
        )
    else:
        assert set(results.data_vars) == set(expected_cols)
    assert all(col not in results.data_vars for col in forbidden_cols)
