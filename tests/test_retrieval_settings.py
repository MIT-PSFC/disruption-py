#!/usr/bin/env python3

import os

import pytest

from disruption_py.inout.mds import ProcessMDSConnection
from disruption_py.machine.tokamak import Tokamak
from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.workflow import get_shots_data


def skip_on_fast_execution(method):
    if "GITHUB_ACTIONS" in os.environ:

        @pytest.mark.skip("fast execution")
        def wrapper(method):
            return method

        return wrapper
    return method


@pytest.fixture(scope="module")
def full_time_domain_data(tokamak, shotlist):
    retrieval_settings = RetrievalSettings(
        efit_nickname_setting="default", domain_setting="full"
    )
    results = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        num_processes=2,
    )
    return results


@pytest.mark.parametrize("num_processes", [1, 2])
def test_sql_cache(tokamak, shotlist, num_processes):
    """
    Use `i_efc` to test retrieving cached data from SQL. `i_efc` exists in SQL and
    it is the only parameter returned from its physics method, so the physics
    method will not to run. This test uses a dummy MDSconnection to ensure we don't
    call MDSplus.
    """
    # `i_efc` does not exist on DIII-D
    if tokamak != Tokamak.CMOD:
        pytest.skip()

    retrieval_settings = RetrievalSettings(
        cache_setting="sql",
        use_cache_setting_timebase=True,
        run_columns=["i_efc"],
        run_tags=[],
        only_requested_columns=True,
        efit_nickname_setting="default",
    )

    def dummy_mds_initializer():
        return ProcessMDSConnection(ProcessMDSConnection.DUMMY_CONNECTION_STRING)

    results = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        num_processes=num_processes,
        mds_connection_initializer=dummy_mds_initializer,
    )
    for res in results:
        assert set(res.columns) == {"i_efc", "shot", "time", "commit_hash"}


@skip_on_fast_execution
@pytest.mark.parametrize("num_processes", [1, 2])
def test_only_requested_columns(tokamak, shotlist, num_processes):
    """
    Ensure `only_requested_columns` works. `v_loop` is returned by
    `get_ohmic_parameters`, so we should not see `p_oh` returned. `q95` is from
    efit, so none of the other efit quantities should be returned.
    """
    retrieval_settings = RetrievalSettings(
        run_columns=["v_loop", "q95"],
        run_tags=[],
        run_methods=[],
        only_requested_columns=True,
    )
    results = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        num_processes=num_processes,
    )
    for res in results:
        assert set(res.columns) == {"v_loop", "q95", "shot", "time", "commit_hash"}


@skip_on_fast_execution
@pytest.mark.parametrize("domain_setting", ["flattop", "rampup_and_flattop"])
def test_domain_setting(tokamak, shotlist, domain_setting, full_time_domain_data):
    """
    Test the two partial domain settings by comparing their start and end times
    with the full domain.
    """
    retrieval_settings = RetrievalSettings(
        efit_nickname_setting="default", domain_setting=domain_setting
    )
    results = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        num_processes=2,
    )
    for part_domain, full_domain in zip(results, full_time_domain_data):
        # Both partial domains should end before the full domain, but only flattop
        # starts after
        if domain_setting == "flattop":
            assert "start", part_domain["time"][0] > full_domain["time"][0]
        assert "end", part_domain["time"].values[-1] < full_domain["time"].values[-1]
