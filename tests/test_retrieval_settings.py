#!/usr/bin/env python3

"""
Unit tests for the retrieval settings which covers data retrieval from various
sources and the time domain of the data retrieved.
"""

import pandas as pd
import pytest

from disruption_py.core.physics_method.runner import REQUIRED_COLS
from disruption_py.core.utils.misc import safe_df_concat
from disruption_py.inout.mds import ProcessMDSConnection
from disruption_py.machine.tokamak import Tokamak
from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.workflow import get_shots_data
from tests.conftest import skip_on_fast_execution
from tests.utils.data_difference import assert_frame_equal_unordered


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
        output_setting="list",
    )

    # Ensure there is a connection to SQL -- if there is a dummy MDS connection,
    # but no cache data is retrieved from SQL, then results=[]
    assert len(shotlist) == len(results)

    # Verify the correct columns were retrieved from SQL
    for res in results:
        assert {"time_until_disrupt", "shot", "time", "commit_hash"} == set(res.columns)


@skip_on_fast_execution
@pytest.mark.parametrize("output_format", [".csv", ".hdf5"])
def test_cache_setting_prev_output(tokamak, shotlist, test_file_path_f, output_format):
    """
    Use the file output from an initial call to `get_shots_data` as the cache for
    a subsequent call to `get_shots_data` and make sure the data remains the same.
    """
    # Save data to file
    get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        num_processes=2,
        output_setting=test_file_path_f(output_format),
    )

    if output_format == ".csv":
        cache_data = pd.read_csv(test_file_path_f(".csv"), dtype={"commit_hash": str})
    else:
        cache_data_list = []
        for shot in shotlist:
            cache_data_list.append(
                pd.read_hdf(test_file_path_f(".hdf5"), key=f"df_{shot}")
            )
        cache_data = safe_df_concat(pd.DataFrame(), cache_data_list)

    # Use saved data as cache
    retrieval_settings = RetrievalSettings(
        cache_setting=cache_data,
        use_cache_setting_timebase=True,
    )

    results = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        num_processes=2,
        output_setting="dataframe",
        mds_connection_initializer=dummy_mds_initializer,
    )
    assert len(shotlist) == len(set(cache_data["shot"])) == len(set(results["shot"]))
    assert_frame_equal_unordered(cache_data, results)


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
        output_setting="list",
        log_settings="WARNING",
    )
    for res in results:
        assert {"ip", "q95", "shot", "time", "commit_hash"} == set(res.columns)


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
    assert not any(col in results.columns for col in forbidden_cols)
