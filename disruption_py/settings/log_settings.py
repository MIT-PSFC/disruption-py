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

    Attributes
    ----------
    file_path : str, optional
        Path to the log file. If None, no log file will be created.
        By default, a log file will be created in a temporary folder.
    file_level : str
        Logging level for the log file (default is "DEBUG").
        Possible values are:
        "TRACE", "DEBUG", "VERBOSE" (custom), "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL".
        See: https://loguru.readthedocs.io/en/stable/api/logger.html#levels
    console_level : str or int, optional
        The log level for the console. Default is None, so log level will be determined
        dynamically based on the number of shots.
        Possible values are:
        "TRACE", "DEBUG", "VERBOSE" (custom), "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL".
        See: https://loguru.readthedocs.io/en/stable/api/logger.html#levels
    warning_threshold : int
        If number of shots is greater than this threshold, the console log level will
        be "WARNING". Default is 1000.
    success_threshold : int
        If number of shots is greater than this threshold and less than the warning_threshold,
        the console log level will be "SUCCESS". Default is 500.
    info_threshold : int
        If number of shots is greater than this threshold and less than the success_threshold,
        the console log level will be "INFO". Default is 50.
    _logging_has_been_setup : bool
        Internal flag to prevent multiple setups (default is False).
    """

    file_path: str = os.path.join(get_temporary_folder(), "output.log")
    file_level: str = "DEBUG"
    console_level: str = None

    warning_threshold: int = 1000
    success_threshold: int = 500
    info_threshold: int = 50

    _logging_has_been_setup: bool = False

    def reset_handlers(self, num_shots: int = None):
        """
        Remove default logger and set up custom handlers.

        Parameters
        ----------
        num_shots : int, optional
            Number of shots to determine the console log level dynamically.
        """
        # Remove default logger
        logger.remove()

        # formats
        message_format = "<level>[{level:^7s}] {message}</level>"
        console_format = "{time:HH:mm:ss.SSS} " + message_format
        file_format = "{time:YYYY-MM-DD HH:mm:ss.SSS} " + message_format

        if self.console_level is None:
            # Determine console log level dynamically based on the number of shots
            console_level = "VERBOSE"
            if num_shots and num_shots > self.warning_threshold:
                console_level = "WARNING"
            elif num_shots and num_shots > self.success_threshold:
                console_level = "SUCCESS"
            elif num_shots and num_shots > self.info_threshold:
                console_level = "INFO"
        elif isinstance(self.console_level, str):
            console_level = self.console_level.upper()
        else:
            console_level = self.console_level

        # Add console handler
        logger.add(
            lambda msg: tqdm.write(msg, end=""),
            level=console_level,
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

        self.reset_handlers(num_shots=None)

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
