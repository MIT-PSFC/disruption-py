#!/usr/bin/env python3

"""
Unit tests for nickname_setting.

These tests are pure-Python with a fake data connection — no MDSplus or
network access required, so they run under both the fast and full test suites.
"""

import pytest

from disruption_py.core.physics_method.errors import FetchDataError
from disruption_py.inout.mds import mdsExceptions
from disruption_py.machine.tokamak import Tokamak
from disruption_py.settings.nickname_setting import (
    DefaultNicknameSetting,
    DisruptionNicknameSetting,
    NicknameSetting,
    NicknameSettingDict,
    NicknameSettingList,
    NicknameSettingParams,
    StaticNicknameSetting,
    resolve_nickname_setting,
)


class FakeMDSConn:
    """Minimal stand-in for MDSConnection with controllable open_tree behavior."""

    def __init__(self, openable=None):
        # set of tree names that open_tree() should accept; all others raise
        self.openable = set(openable) if openable else set()
        self.opened = []  # ordered record of every open_tree() call

    def open_tree(self, tree_name):
        """Record the attempt; raise TreeFOPENR if the tree is not in the openable set."""
        self.opened.append(tree_name)
        if tree_name not in self.openable:
            # mdsthin's TreeFOPENR.__init__ takes no message arg; instantiate bare.
            raise mdsExceptions.TreeFOPENR()


def _params(data_conn, *, shot_id=1234567, disruption_time=None, tokamak=Tokamak.CMOD):
    """Build NicknameSettingParams with a fake data_conn and no database."""
    return NicknameSettingParams(
        shot_id=shot_id,
        data_conn=data_conn,
        database=None,
        disruption_time=disruption_time,
        tokamak=tokamak,
    )


# ---------------------------------------------------------------------------
# resolve_nickname_setting dispatch
# ---------------------------------------------------------------------------


def test_resolve_passthrough_instance():
    """A NicknameSetting instance is returned unchanged."""
    instance = StaticNicknameSetting("efit21")
    assert resolve_nickname_setting(instance) is instance


def test_resolve_string_registered_key():
    """Registered string keys map to their canonical class."""
    assert isinstance(resolve_nickname_setting("disruption"), DisruptionNicknameSetting)
    assert isinstance(resolve_nickname_setting("default"), DefaultNicknameSetting)


def test_resolve_string_static():
    """An unregistered single string becomes a StaticNicknameSetting."""
    out = resolve_nickname_setting("efit21")
    assert isinstance(out, StaticNicknameSetting)
    assert out.tree_name == "efit21"


def test_resolve_string_comma_cascade():
    """A comma-separated string of literal tree names becomes a NicknameSettingList."""
    out = resolve_nickname_setting("efit21,efit18,efit01")
    assert isinstance(out, NicknameSettingList)
    assert out.resolved_items == ["efit21", "efit18", "efit01"]


def test_resolve_string_comma_strips_whitespace():
    """Whitespace and empty entries in the comma string are tolerated."""
    out = resolve_nickname_setting(" efit21 , efit18 ,, efit01 ")
    assert out.resolved_items == ["efit21", "efit18", "efit01"]


def test_resolve_list_cascade():
    """A list dispatches to NicknameSettingList."""
    out = resolve_nickname_setting(["efit21", "analysis"])
    assert isinstance(out, NicknameSettingList)


def test_resolve_dict_dispatch():
    """A dict dispatches to NicknameSettingDict."""
    out = resolve_nickname_setting({Tokamak.CMOD: "efit21"})
    assert isinstance(out, NicknameSettingDict)


def test_resolve_dict_value_can_be_cascade():
    """A NicknameSettingDict value may itself be a cascade list."""
    out = resolve_nickname_setting({Tokamak.CMOD: ["efit21", "efit18", "analysis"]})
    assert isinstance(out, NicknameSettingDict)
    inner = out.resolved_nickname_setting_dict[Tokamak.CMOD]
    assert isinstance(inner, NicknameSettingList)
    assert inner.resolved_items == ["efit21", "efit18", "analysis"]


def test_resolve_invalid_type():
    """Unsupported types raise ValueError."""
    with pytest.raises(ValueError, match="Invalid nickname setting type"):
        resolve_nickname_setting(42)


# ---------------------------------------------------------------------------
# NicknameSettingList item validation
# ---------------------------------------------------------------------------


def test_list_empty_rejected():
    """An empty cascade is rejected at construction."""
    with pytest.raises(ValueError, match="at least one item"):
        NicknameSettingList([])


@pytest.mark.parametrize("bad_item", [42, 3.14, None, {"a": 1}, ["nested"]])
def test_list_bad_item_type(bad_item):
    """Cascade items must be str or NicknameSetting."""
    with pytest.raises(ValueError, match="must be str or NicknameSetting"):
        NicknameSettingList([bad_item])


def test_list_accepts_nicknamesetting_instances():
    """Cascade items may be NicknameSetting instances."""
    cascade = NicknameSettingList([StaticNicknameSetting("efit21"), "analysis"])
    assert isinstance(cascade.resolved_items[0], StaticNicknameSetting)
    assert cascade.resolved_items[1] == "analysis"


# ---------------------------------------------------------------------------
# NicknameSettingList cascade behavior
# ---------------------------------------------------------------------------


def test_cascade_first_succeeds():
    """First candidate opens cleanly → no fallback, no extra open_tree calls."""
    conn = FakeMDSConn(openable={"efit21"})
    cascade = NicknameSettingList(["efit21", "efit18", "analysis"])
    assert cascade.get_tree_name(_params(conn)) == "efit21"
    assert conn.opened == ["efit21"]


def test_cascade_falls_through_to_second():
    """First fails with TreeFOPENR → cascade tries the next."""
    conn = FakeMDSConn(openable={"efit18"})
    cascade = NicknameSettingList(["efit21", "efit18", "analysis"])
    assert cascade.get_tree_name(_params(conn)) == "efit18"
    assert conn.opened == ["efit21", "efit18"]


def test_cascade_falls_through_to_last():
    """Cascade walks all the way to the final candidate when needed."""
    conn = FakeMDSConn(openable={"analysis"})
    cascade = NicknameSettingList(["efit21", "efit18", "analysis"])
    assert cascade.get_tree_name(_params(conn)) == "analysis"
    assert conn.opened == ["efit21", "efit18", "analysis"]


def test_cascade_all_fail_raises():
    """If no candidate opens, FetchDataError is raised listing every attempt."""
    conn = FakeMDSConn(openable=set())
    cascade = NicknameSettingList(["efit21", "efit18", "analysis"])
    with pytest.raises(FetchDataError, match="No tree in nickname cascade"):
        cascade.get_tree_name(_params(conn))
    assert conn.opened == ["efit21", "efit18", "analysis"]


def test_cascade_with_nested_nicknamesetting():
    """A NicknameSetting instance inside the cascade resolves at access time."""
    conn = FakeMDSConn(openable={"analysis"})
    cascade = NicknameSettingList([StaticNicknameSetting("efit21"), "analysis"])
    assert cascade.get_tree_name(_params(conn)) == "analysis"
    assert conn.opened == ["efit21", "analysis"]


# ---------------------------------------------------------------------------
# DisruptionNicknameSetting CMOD path (post-refactor)
# ---------------------------------------------------------------------------


def test_disruption_cmod_picks_configured_tree_when_available():
    """When the configured tree (efit18 under pytest) opens, it wins."""
    conn = FakeMDSConn(openable={"efit18"})
    setting = DisruptionNicknameSetting()
    # pytest is in sys.modules during the test run, so tree is forced to efit18
    assert setting.get_tree_name(_params(conn)) == "efit18"


def test_disruption_cmod_falls_back_to_analysis():
    """When efit18 isn't openable, the CMOD path falls back to 'analysis'."""
    conn = FakeMDSConn(openable={"analysis"})
    setting = DisruptionNicknameSetting()
    assert setting.get_tree_name(_params(conn)) == "analysis"
    # confirm the cascade actually attempted efit18 first
    assert conn.opened == ["efit18", "analysis"]


# ---------------------------------------------------------------------------
# Comma-string tokens that match registered keys resolve to their classes
# ---------------------------------------------------------------------------


def test_resolve_comma_string_resolves_registered_keys():
    """Tokens in a comma-string that match registered keys map to instances."""
    out = resolve_nickname_setting("disruption,doesnotexist,analysis")
    assert isinstance(out, NicknameSettingList)
    assert isinstance(out.resolved_items[0], DisruptionNicknameSetting)
    assert out.resolved_items[1] == "doesnotexist"
    # "analysis" is a registered (deprecated) alias for DefaultNicknameSetting.
    assert isinstance(out.resolved_items[2], DefaultNicknameSetting)


# ---------------------------------------------------------------------------
# Cascade catches nested-cascade exhaustion (FetchDataError)
# ---------------------------------------------------------------------------


class _AlwaysFailsNickname(NicknameSetting):
    """Test helper: raises FetchDataError from get_tree_name regardless of params."""

    def _get_tree_name(self, params):
        raise FetchDataError("test: nested cascade exhausted")


def test_cascade_falls_through_on_nested_fetch_data_error():
    """A NicknameSetting item raising FetchDataError counts as a failed candidate."""
    conn = FakeMDSConn(openable={"analysis"})
    cascade = NicknameSettingList([_AlwaysFailsNickname(), "analysis"])
    assert cascade.get_tree_name(_params(conn)) == "analysis"
    # The failing nickname never reached open_tree, so only "analysis" was attempted.
    assert conn.opened == ["analysis"]
