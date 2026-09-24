#!/usr/bin/env python3

"""This module contains tests to ensure all of the config settings load properly."""

import os

from disruption_py.config import config, configs
from disruption_py.settings import LogSettings
from disruption_py.settings.log_settings import level_no, resolve_file_level


def change_directory(test, tmpdir="/tmp"):
    """
    Change the current working directory before a test and revert back to the
    original directory after the test completes.
    """

    def wrapper():
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        test()
        os.chdir(original_dir)

    return wrapper


@change_directory
def test_settings_file():
    """
    Temporarily change the current working directory to test if the config settings
    file is reachable.
    """
    assert config().time.time_const is not None


def test_access_tokamak_settings(tokamak):
    """
    Test each tokamak's unique settings are accessible.
    """
    assert config(tokamak)


def test_log_settings_defaults_from_config():
    """
    Unset LogSettings fields are resolved from the [log] config section.
    """
    settings = LogSettings()
    assert settings.console_level == config().log.console_level
    assert settings.file_level == resolve_file_level(
        config().log.file_level, config().log.console_level
    )


def test_log_settings_env_override(monkeypatch):
    """
    DISPY_LOG__* environment variables override the config file defaults.
    """
    monkeypatch.setenv("DISPY_LOG__CONSOLE_LEVEL", "WARNING")
    # drop the cached config so the environment variables are re-read
    configs.pop("default", None)
    try:
        settings = LogSettings()
        assert settings.console_level == "WARNING"
        # untouched fields still come from the layered configuration
        assert settings.file_level == config().log.file_level
    finally:
        # rebuild the cache without the overrides for subsequent tests
        configs.pop("default", None)


def test_log_settings_argument_wins():
    """
    Explicitly passed values take precedence over the configuration.
    """
    settings = LogSettings(console_level="warning", file_level="trace")
    assert settings.console_level == "WARNING"
    assert settings.file_level == "TRACE"
    # unset fields still resolve from the configuration
    settings = LogSettings(console_level="WARNING")
    assert settings.file_level == config().log.file_level


def test_log_settings_to_config():
    """
    The effective settings round-trip into the shape of the [log] section.
    """
    settings = LogSettings(console_level="warning")
    log_config = settings.to_config()
    # level names are normalized to the spelling loguru expects
    assert log_config["console_level"] == "WARNING"
    assert log_config["file_level"] == config().log.file_level
    assert set(log_config) == {"file_level", "console_level"}
    # written back, the configuration reflects the runtime value
    config("cmod").update({"log": log_config})
    try:
        assert config("cmod").log.console_level == "WARNING"
        assert config("cmod").log.file_level == config().log.file_level
    finally:
        configs.pop("cmod", None)


def test_file_level_floor():
    """
    The file log level is at least as verbose as the console level.
    """
    assert level_no("debug") == level_no("DEBUG") == 10
    # a more verbose console lowers the file level to match
    assert resolve_file_level("DEBUG", "TRACE") == "TRACE"
    # a quieter console never raises it
    assert resolve_file_level("DEBUG", "WARNING") == "DEBUG"
    # numeric levels are accepted as is
    assert resolve_file_level("DEBUG", 5) == 5
    # the floor is applied at construction, so the dump carries it too
    settings = LogSettings(file_level="DEBUG", console_level="trace")
    assert settings.file_level == "TRACE"
    assert settings.to_config()["file_level"] == "TRACE"
    settings = LogSettings(file_level="DEBUG", console_level="WARNING")
    assert settings.file_level == "DEBUG"
