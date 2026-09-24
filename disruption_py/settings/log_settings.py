#!/usr/bin/env python3

"""
This module defines the LogSettings class, which provides settings and setup for
logging in both files and console with customizable levels and formats.
"""

import multiprocessing
import os
import sys
from dataclasses import dataclass
from functools import partialmethod
from typing import Union

from loguru import logger
from tqdm.auto import tqdm

from disruption_py.config import config
from disruption_py.core.utils.misc import get_metadata, get_temporary_folder

# Register logger state at module import time so it is available in every
# process, including multiprocessing workers that use `spawn` (macOS/Windows)
# or `forkserver` (Linux, Python 3.14+) start methods. These do not inherit
# the main process's loguru state, so anything the workers need must live
# here — not in setup_logging(), which only runs in the main process.
if not hasattr(logger.__class__, "verbose"):
    try:
        logger.level("VERBOSE", no=15, color="<dim>")
    except TypeError:
        # Level already exists (re-import), just update the color.
        logger.level("VERBOSE", color="<dim>")
    logger.__class__.verbose = partialmethod(logger.__class__.log, "VERBOSE")

# Custom colors for built-in levels. Must run at module level so workers
# inherit them via re-import under spawn/forkserver.
logger.level("TRACE", color="<cyan><dim>")
logger.level("DEBUG", color="<blue>")
logger.level("VERBOSE", color="<dim>")
logger.level("INFO", color="")
logger.level("SUCCESS", color="<green>")
logger.level("WARNING", color="<yellow>")
logger.level("ERROR", color="<red>")

LogSettingsType = Union["LogSettings", str, int]


@dataclass
class LogSettings:
    """
    Settings for configuring logging.

    Unset levels are read from the `[log]` section of the layered configuration.
    In order of precedence: values passed explicitly, then `DISPY_LOG__*`
    environment variables (e.g. `DISPY_LOG__CONSOLE_LEVEL=WARNING`), then
    `~/.config/disruption-py/user.toml`, then the repo `config.toml`.

    The file log level is never less verbose than the console: if the console
    level resolves below `file_level` (e.g. "TRACE" against the default
    "DEBUG"), `file_level` is lowered to match. The floor is applied once at
    construction, so the object, its handlers and the configuration dump all
    carry the same effective values.

    Attributes
    ----------
    file_path : str, optional
        Path to the log file. If None, no log file will be created.
        By default, a log file will be created in a temporary folder.
    file_level : str or int, optional
        Logging level for the log file. Default is None, so the level is read
        from the configuration (`log.file_level`, "DEBUG" out of the box).
        A value less verbose than the console level is lowered to match it.
        Possible values are:
        "TRACE", "DEBUG", "VERBOSE" (custom), "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL".
        See: https://loguru.readthedocs.io/en/stable/api/logger.html#levels
    console_level : str or int, optional
        The log level for the console. Default is None, so the level is read
        from the configuration (`log.console_level`, "INFO" out of the box).
        Possible values are:
        "TRACE", "DEBUG", "VERBOSE" (custom), "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL".
        See: https://loguru.readthedocs.io/en/stable/api/logger.html#levels
    _logging_has_been_setup : bool
        Internal flag to prevent multiple setups (default is False).
    """

    file_path: str = os.path.join(get_temporary_folder(), "output.log")
    file_level: str | int = None
    console_level: str | int = None

    _logging_has_been_setup: bool = False

    def __post_init__(self):
        # Resolve unset levels from the layered config, where DISPY_LOG__*
        # environment variables override user.toml, which overrides the repo
        # config.toml. Explicit arguments always win.
        log_config = config().log
        if self.file_level is None:
            self.file_level = log_config.file_level
        if self.console_level is None:
            self.console_level = log_config.console_level
        # normalize level names once, so that the object and the config dump
        # carry the same spelling that loguru expects
        if isinstance(self.file_level, str):
            self.file_level = self.file_level.upper()
        if isinstance(self.console_level, str):
            self.console_level = self.console_level.upper()
        # the file is at least as verbose as the console
        self.file_level = resolve_file_level(self.file_level, self.console_level)

    def to_config(self) -> dict:
        """
        Return the effective settings in the shape of the `[log]` config section,
        so that they can be written back into the configuration at runtime.

        Returns
        -------
        dict
            The resolved log configuration.
        """
        return {
            "file_level": self.file_level,
            "console_level": self.console_level,
        }

    def reset_handlers(self):
        """
        Remove default logger and set up custom handlers.
        """
        # Remove default logger
        logger.remove()

        # formats
        message_format = "<level>[{level:^7s}] {message}</level>"
        console_format = "{time:HH:mm:ss.SSS} " + message_format
        file_format = "{time:YYYY-MM-DD HH:mm:ss.SSS} " + message_format

        # Add console handler
        logger.add(
            lambda msg: tqdm.write(msg, end=""),
            level=self.console_level,
            format=console_format,
            colorize=True,
            enqueue=True,
            backtrace=False,
            diagnose=True,
        )

        # Add file handler if log file path is provided. Main process truncates
        # to start a fresh log; workers append so they don't wipe the main
        # process's output under spawn/forkserver.
        if self.file_path is not None:
            is_main = multiprocessing.current_process().name == "MainProcess"
            logger.add(
                self.file_path,
                level=self.file_level,
                format=file_format,
                mode="w" if is_main else "a",
                enqueue=True,
                backtrace=False,
                diagnose=True,
            )

    def setup_logging(self):
        """
        Set up logging with custom styles and levels.
        """
        if self._logging_has_been_setup:
            return

        self.reset_handlers()

        # header
        metadata = get_metadata()
        commit = metadata.get("commit")
        logger.info(
            "Starting: {package} ~ v{version}{strcommit} / {user}@{host}",
            strcommit=f" # {commit}" if commit else "",
            **metadata,
        )
        if self.file_path is not None:
            logger.info("Logging: {l}", l=self.file_path)
        logger.debug(
            "Source: {source}",
            source=metadata["source"],
        )
        logger.debug("Executable: {e}", e=sys.executable)

        self._logging_has_been_setup = True


def level_no(level: str | int) -> int:
    """
    Convert a log level name or number to its loguru severity number.

    Parameters
    ----------
    level : str | int
        Level name (e.g. "DEBUG") or number.

    Returns
    -------
    int
        The severity number.
    """
    if isinstance(level, int):
        return level
    return logger.level(level.upper()).no


def resolve_file_level(file_level: str | int, console_level: str | int) -> str | int:
    """
    Resolve the file log level so that it is at least as verbose as the console.

    Parameters
    ----------
    file_level : str | int
        The configured file log level.
    console_level : str | int
        The effective console log level.

    Returns
    -------
    str | int
        The more verbose of the two levels, as given.
    """
    if level_no(console_level) < level_no(file_level):
        return console_level
    return file_level


def resolve_log_settings(
    log_settings: LogSettingsType,
) -> LogSettings:
    """
    Resolve the log settings to a LogSettings instance.

    Parameters
    ----------
    log_settings : LogSettingsType
        The log setting to resolve, which can be an instance of LogSettings, or
        a string or int representing the console log level

    Returns
    -------
    LogSettings
        The resolved LogSettings instance.
    """
    if isinstance(log_settings, LogSettings):
        return log_settings

    if isinstance(log_settings, (str, int)):
        return LogSettings(console_level=log_settings)

    if isinstance(log_settings, dict):
        return LogSettings(**log_settings)

    if log_settings is None:
        return LogSettings()

    raise ValueError(f"Invalid log settings {log_settings}")
