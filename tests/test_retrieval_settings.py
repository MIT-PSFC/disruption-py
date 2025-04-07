#!/usr/bin/env python3

"""
Unit tests for the retrieval settings which covers data retrieval from various
sources and the time domain of the data retrieved.
"""

import pytest

from disruption_py.core.physics_method.runner import REQUIRED_COLS
from disruption_py.inout.mds import ProcessMDSConnection
from disruption_py.machine.tokamak import Tokamak
from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.workflow import get_shots_data
from tests.conftest import skip_on_fast_execution


def dummy_mds_initializer():
    """Return a dummy MDS connection to ensure no fresh data is retrieved."""
    return ProcessMDSConnection(None)


@pytest.fixture(scope="module", name="full_time_domain_data")
def full_time_domain_data_fixture(tokamak, shotlist):
    """Get data for the full time domain"""
    retrieval_settings = RetrievalSettings(
        efit_nickname_setting="default", domain_setting="full"
    )
    results = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        output_setting="dict",
        num_processes=2,
    )
    # Output data in the order shots were given
    results = [results[shot] for shot in shotlist]
    return results


@skip_on_fast_execution
def test_only_requested_columns(tokamak, shotlist):
    """
    Ensure `only_requested_columns` works. `ip` is returned by
    `get_ip_parameters`, so we should not see any of the other quantities like
    `dip_dt` returned. `q95` is from efit, so none of the other efit quantities
    should be returned.
    """
    retrieval_settings = RetrievalSettings(
        run_columns=["ip", "q95"],
        only_requested_columns=True,
    )
    results = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        num_processes=2,
        output_setting="dataframe",
        log_settings="WARNING",
    )
    assert {"ip", "q95", "shot", "time"} == set(results.columns)


@skip_on_fast_execution
@pytest.mark.parametrize("domain_setting", ["flattop", "rampup_and_flattop"])
def test_domain_setting(tokamak, shotlist, domain_setting, full_time_domain_data):
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
    results = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        output_setting="dict",
        num_processes=2,
        log_settings="WARNING",
    )
    results = [results[shot] for shot in shotlist]

    assert len(shotlist) == len(results) == len(full_time_domain_data)
    for part_domain, full_domain in zip(results, full_time_domain_data):
        p_start, p_end = part_domain["time"].values[0], part_domain["time"].values[-1]
        f_start, f_end = full_domain["time"].values[0], full_domain["time"].values[-1]
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
        ([], None, REQUIRED_COLS, []),
        (["~get_kappa_area"], None, None, ["kappa_area"]),
        (["get_kappa_area"], None, REQUIRED_COLS | {"kappa_area"}, []),
        # Test run_columns with run_methods=None
        (None, [], REQUIRED_COLS, []),
        (None, ["kappa_area"], REQUIRED_COLS | {"kappa_area"}, []),
        # Test run_methods and run_columns combo
        ([], [], REQUIRED_COLS, []),
        (["get_kappa_area"], [], REQUIRED_COLS | {"kappa_area"}, []),
        (
            ["get_kappa_area"],
            ["greenwald_fraction"],
            REQUIRED_COLS | {"kappa_area", "n_e", "dn_dt", "greenwald_fraction"},
            [],
        ),
        (
            ["~get_kappa_area"],
            ["greenwald_fraction"],
            REQUIRED_COLS | {"n_e", "dn_dt", "greenwald_fraction"},
            [],
        ),
        (["~get_kappa_area"], ["kappa_area"], REQUIRED_COLS, []),
    ],
)
def test_run_methods_and_columns(
    tokamak,
    shotlist,
    run_methods,
    run_columns,
    expected_cols,
    forbidden_cols,
    full_time_domain_data,
):
    """
    Test the `run_methods` and `run_columns` parameters of RetrievalSettings.

    - If both are None, all methods are run
    - If one is None or [], the methods/columns specified by the other are run
    - If both are specified, the combined methods/columns are run
    - If `run_methods` excludes a method returning a column specified in `run_columns`,
      the method is not run
    """
    num_all_cols = len(full_time_domain_data[0].columns)
    retrieval_settings = RetrievalSettings(
        run_methods=run_methods,
        run_columns=run_columns,
    )
    results = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        num_processes=2,
        log_settings="WARNING",
    )
    # Expected columns None means all columns (except forbidden cols) are returned
    if expected_cols is None:
        assert len(results.columns) == num_all_cols - len(forbidden_cols)
    else:
        assert set(results.columns) == set(expected_cols)
    assert all(col not in results.columns for col in forbidden_cols)
