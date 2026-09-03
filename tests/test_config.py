#!/usr/bin/env python3

"""This module contains tests to ensure all of the config settings load properly."""

import os

from disruption_py.config import config, configs
from disruption_py.settings import LogSettings


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
    assert settings.file_level == config().log.file_level
    assert settings.warning_threshold == config().log.warning_threshold
    assert settings.success_threshold == config().log.success_threshold
    assert settings.info_threshold == config().log.info_threshold


def test_log_settings_env_override(monkeypatch):
    """
    DISPY_LOG__* environment variables override the config file defaults.
    """
    monkeypatch.setenv("DISPY_LOG__WARNING_THRESHOLD", "2000")
    monkeypatch.setenv("DISPY_LOG__FILE_LEVEL", "INFO")
    # drop the cached config so the environment variables are re-read
    configs.pop("default", None)
    try:
        settings = LogSettings()
        assert settings.warning_threshold == 2000
        assert settings.file_level == "INFO"
        # untouched fields still come from the config file
        assert settings.success_threshold == 500
    finally:
        # rebuild the cache without the overrides for subsequent tests
        configs.pop("default", None)


def test_log_settings_argument_wins():
    """
    Explicitly passed values take precedence over the configuration.
    """
    settings = LogSettings(warning_threshold=7, file_level="ERROR")
    assert settings.warning_threshold == 7
    assert settings.file_level == "ERROR"
    # unset fields still resolve from the configuration
    assert settings.info_threshold == config().log.info_threshold
