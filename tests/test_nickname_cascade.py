#!/usr/bin/env python3

"""
Integration test for the efit nickname cascade (``NicknameSettingList``, PR #559).

CMOD-specific: exercises the ``[efit18, analysis, efit21]`` cascade against shots
that do and do not have an ``efit18`` tree, and asserts on the retrieved output.
Requires MDSplus + SQL on the MFE workstation; skipped on fast (GITHUB_ACTIONS)
execution and on non-CMOD tokamaks.

Builds on the snippet @yumouwei posted in the PR thread.
"""

import numpy as np
import pytest
import xarray as xr

from disruption_py.machine.tokamak import Tokamak
from disruption_py.settings import RetrievalSettings
from disruption_py.workflow import get_shots_data
from tests.conftest import skip_on_fast_execution

# CMOD shots chosen for their efit tree availability:
DISRUPTIVE = 1150805012  # has efit18, analysis, and efit21
NON_DISRUPTIVE = 1150805013  # has analysis and efit21, but NO efit18
SHOTLIST = [DISRUPTIVE, NON_DISRUPTIVE]


def _fetch(efit_nickname_setting) -> xr.Dataset:
    """Retrieve ``li`` on the efit timebase for one efit nickname setting."""
    retrieval_settings = RetrievalSettings(
        efit_nickname_setting=efit_nickname_setting,
        run_columns=["li"],
        time_setting="efit",
    )
    return get_shots_data(
        tokamak=Tokamak.CMOD,
        shotlist_setting=SHOTLIST,
        retrieval_settings=retrieval_settings,
    )


def _rows_for(ds: xr.Dataset, shot: int) -> xr.Dataset:
    """Slice of the flat (idx,) dataset belonging to one shot.

    The default 'dataset' output concatenates every shot along ``idx`` with
    ``shot``/``time`` as non-dimension coords, so we mask on shot rather than
    use ``.sel(shot=...)`` (shot is not a dimension coordinate).
    """
    return ds.isel(idx=np.where(ds["shot"].values == shot)[0])


def _length(ds: xr.Dataset, shot: int) -> int:
    """Number of timeslices retrieved for a shot (0 if the shot is absent)."""
    return int((ds["shot"].values == shot).sum())


@pytest.fixture(scope="module", autouse=True)
def _require_cmod(tokamak):
    """Cascade tree availability is CMOD-specific; skip elsewhere."""
    if tokamak is not Tokamak.CMOD:
        pytest.skip(f"efit cascade test is CMOD-specific (got {tokamak.name})")


@pytest.fixture(scope="module", name="cascade")
def cascade_fixture() -> xr.Dataset:
    """Cascade preferring efit18, falling back to analysis, then efit21."""
    return _fetch(["efit18", "analysis", "efit21"])


@pytest.fixture(scope="module", name="efit18")
def efit18_fixture() -> xr.Dataset:
    """Single-tree efit18 baseline."""
    return _fetch("efit18")


@pytest.fixture(scope="module", name="analysis")
def analysis_fixture() -> xr.Dataset:
    """Single-tree analysis baseline."""
    return _fetch("analysis")


@skip_on_fast_execution
def test_cascade_first_item_wins_when_present(cascade, efit18, analysis):
    """When efit18 exists, the cascade selects it and ignores later items."""
    # Timeslice count is the primary difference between efit trees: the cascade
    # matches efit18's length, not analysis's.
    assert _length(efit18, DISRUPTIVE) != _length(analysis, DISRUPTIVE)
    assert _length(cascade, DISRUPTIVE) == _length(efit18, DISRUPTIVE)
    # Disruptive shot has efit18 -> cascade == efit18-only, bit for bit.
    assert _rows_for(cascade, DISRUPTIVE).equals(_rows_for(efit18, DISRUPTIVE))
    # ...and efit18 is a genuinely distinct equilibrium from analysis.
    assert not _rows_for(cascade, DISRUPTIVE).equals(_rows_for(analysis, DISRUPTIVE))


@skip_on_fast_execution
def test_cascade_falls_through_when_missing(cascade, efit18, analysis):
    """When efit18 is absent, the cascade falls through to the next item."""
    # Non-disruptive shot has no efit18 tree at all.
    assert _length(efit18, NON_DISRUPTIVE) == 0
    # ...so the cascade resolves to analysis: same length and same data.
    assert _length(cascade, NON_DISRUPTIVE) == _length(analysis, NON_DISRUPTIVE)
    assert _rows_for(cascade, NON_DISRUPTIVE).equals(
        _rows_for(analysis, NON_DISRUPTIVE)
    )
