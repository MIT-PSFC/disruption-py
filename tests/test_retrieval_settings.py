#!/usr/bin/env python3

"""
Unit tests for the retrieval settings which covers data retrieval from various
sources and the time domain of the data retrieved.
"""

import pytest
import xarray as xr

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
        output_setting="dataset",
        num_processes=2,
    )
    return results


@pytest.mark.parametrize("num_processes", [1, 2])
def test_cache_setting_sql(tokamak, shotlist, num_processes):
    """
    Use `time_until_disrupt` to test retrieving cached data from SQL.
    `time_until_disrupt` exists in SQL and it is the only parameter returned from
    its physics method, so the physics method will not run. This test uses a dummy
    MDSconnection to ensure we don't call MDSplus.
    """
    retrieval_settings = RetrievalSettings(
        cache_setting="sql",
        use_cache_setting_timebase=True,
        run_columns=["time_until_disrupt"],
        only_requested_columns=True,
        efit_nickname_setting="default",
    )

    results = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        num_processes=num_processes,
        mds_connection_initializer=dummy_mds_initializer,
        output_setting="dataset",
    )

    # Ensure there is a connection to SQL -- if there is a dummy MDS connection,
    # but no cache data is retrieved from SQL, then results=[]
    assert len(shotlist) == len(results.shot)

    # Verify the correct data vars were retrieved from SQL
    assert {"time_until_disrupt"} == set(results.data_vars)
    assert {"shot", "time"} == set(results.dims)
    assert "commit_hash" in results.attrs


@skip_on_fast_execution
@pytest.mark.parametrize("output_format", [".nc", ".hdf5"])
def test_cache_setting_prev_output(tokamak, shotlist, test_file_path_f, output_format):
    """
    Use the file output from an initial call to `get_shots_data` as the cache for
    a subsequent call to `get_shots_data` and make sure the data remains the same.

    Only request `time_until_disrupt` rather than all data vars to speed up the test.
    """
    # Save data to file
    get_shots_data(
        tokamak=tokamak,
        retrieval_settings=RetrievalSettings(
            run_columns=["time_until_disrupt"], only_requested_columns=True
        ),
        shotlist_setting=shotlist,
        num_processes=2,
        output_setting=test_file_path_f(output_format),
    )

    cache_data = xr.open_dataset(test_file_path_f(output_format))

    # Use saved data as cache
    retrieval_settings = RetrievalSettings(
        cache_setting=cache_data,
        use_cache_setting_timebase=True,
        run_columns=["time_until_disrupt"],
        only_requested_columns=True,
    )

    results = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        num_processes=2,
        output_setting="dataset",
        mds_connection_initializer=dummy_mds_initializer,
    )
    xr.testing.assert_equal(cache_data, results)


@skip_on_fast_execution
def test_only_requested_columns(tokamak, shotlist):
    """
    Ensure `only_requested_columns` works. `ip` is returned by
    `get_ip_parameters`, so we should not see any of the other quantities like
    `dip_dt` returned. `q95` is from efit, so none of the other efit quantities
    should be returned.
    """
    columns = {"ip", "q95"}
    retrieval_settings = RetrievalSettings(
        run_columns=columns,
        only_requested_columns=True,
    )
    results = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        num_processes=2,
        output_setting="dataset",
        log_settings="WARNING",
    )
    assert columns == set(results.data_vars)


@skip_on_fast_execution
@pytest.mark.parametrize("domain_setting", ["flattop", "rampup_and_flattop"])
def test_domain_setting(tokamak, shotlist, domain_setting, full_time_domain_data):
    """
    Test the two partial domain settings by comparing their start and end times
    with the full domain.
    """
    if tokamak == Tokamak.D3D and domain_setting == "rampup_and_flattop":
        pytest.skip("rampup_and_flattop domain setting not defined for DIII-D")
    retrieval_settings = RetrievalSettings(
        efit_nickname_setting="default",
        domain_setting=domain_setting,
        run_columns=["ip"],
        only_requested_columns=True,
    )
    results = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        output_setting="dataset",
        num_processes=2,
        log_settings="WARNING",
    )

    assert len(shotlist) == len(results.shot) == len(full_time_domain_data.shot)
    for shot in shotlist:
        part_domain = results.sel(shot=shot).ip.dropna("time")
        full_domain = full_time_domain_data.sel(shot=shot).ip.dropna("time")
        p_start, p_end = part_domain["time"].values[0], part_domain["time"].values[-1]
        f_start, f_end = full_domain["time"].values[0], full_domain["time"].values[-1]
        if domain_setting == "flattop":
            # Use <= because a shot may end during the flattop,
            assert f_start < p_start < p_end <= f_end
        else:
            assert f_start == p_start < p_end <= f_end


@skip_on_fast_execution
@pytest.mark.parametrize(
    "run_methods, run_columns, expected_cols, forbidden_cols",
    [
        # Test run_methods with run_columns=None
        (None, None, None, []),
        ([], None, {}, []),
        (["~get_kappa_area"], None, None, ["kappa_area"]),
        (["get_kappa_area"], None, {"kappa_area"}, []),
        # Test run_columns with run_methods=None
        (None, [], {}, []),
        (None, ["kappa_area"], {"kappa_area"}, []),
        # Test run_methods and run_columns combo
        ([], [], {}, []),
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
        (["~get_kappa_area"], ["kappa_area"], {}, []),
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
    num_all_cols = len(full_time_domain_data.data_vars)
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
        output_setting="dataset",
    )
    # Expected columns None means all columns (except forbidden cols) are returned
    if expected_cols is None:
        assert len(results.data_vars) == num_all_cols - len(forbidden_cols)
    else:
        assert set(results.data_vars) == set(expected_cols)
    assert not any(col in results.data_vars for col in forbidden_cols)
