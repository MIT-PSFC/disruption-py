#!/usr/bin/env python3

"""
This module defines the LogSettings class, which provides settings and setup for
logging in both files and console with customizable levels and formats.
"""

import importlib.metadata
import os
import sys
from dataclasses import dataclass
from functools import partialmethod
from typing import Union

from loguru import logger
from tqdm.auto import tqdm

from disruption_py.core.utils.misc import get_commit_hash, get_temporary_folder

LogSettingsType = Union[
    "LogSettings",
    str,
    int,
]


@dataclass
class LogSettings:
    """
    Settings for configuring logging.

    Attributes
    ----------
    log_file_path : str, optional
        Path to the log file. If None, no log file will be created.
        By default, a log file will be created in a temporary folder.
    file_log_level : str
        Logging level for the log file (default is "DEBUG").
        Possible values are:
        "TRACE", "DEBUG", "VERBOSE" (custom), "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL".
        See: https://loguru.readthedocs.io/en/stable/api/logger.html#levels
    log_file_write_mode : str, optional
        The write mode for the log file. Default is "w".
    log_to_console : bool
        Whether to log messages to the console (default is True).
    console_log_level : str or int, optional
        The log level for the console. Default is None, so log level will be determined
        dynamically based on the number of shots.
        Possible values are:
        "TRACE", "DEBUG", "VERBOSE" (custom), "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL".
        See: https://loguru.readthedocs.io/en/stable/api/logger.html#levels
    use_custom_logging : bool
        Whether to use custom logging. If set to true, no logging setup will be done.
        Default is False.
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

    log_file_path: str = os.path.join(get_temporary_folder(), "output.log")
    file_log_level: str = "DEBUG"
    log_file_write_mode: str = "w"

    log_to_console: bool = True
    console_log_level: str = None

    use_custom_logging: bool = False

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

        if self.console_log_level is None:
            # Determine console log level dynamically based on the number of shots
            console_level = "VERBOSE"
            if num_shots and num_shots > self.warning_threshold:
                console_level = "WARNING"
            elif num_shots and num_shots > self.success_threshold:
                console_level = "SUCCESS"
            elif num_shots and num_shots > self.info_threshold:
                console_level = "INFO"
        elif isinstance(self.console_log_level, str):
            console_level = self.console_log_level.upper()
        else:
            console_level = self.console_log_level

        # Add console handler
        if self.log_to_console:
            logger.add(
                lambda msg: tqdm.write(msg, end=""),
                level=console_level,
                format=console_format,
                colorize=True,
                enqueue=True,
                backtrace=False,
                diagnose=True,
            )

        # Add file handler if log file path is provided
        if self.log_file_path is not None:
            logger.add(
                self.log_file_path,
                level=self.file_log_level,
                format=file_format,
                mode=self.log_file_write_mode,
                enqueue=True,
                backtrace=False,
                diagnose=True,
            )

    def setup_logging(self):
        """
        Set up logging with custom styles and levels.
        """
        if self.use_custom_logging or self._logging_has_been_setup:
            return

        # Set custom style and add a VERBOSE level. This only needs to be done
        # once, so there is no need to add it to the reset_handlers method.
        logger.level("TRACE", color="<cyan><dim>")
        logger.level("DEBUG", color="<blue>")
        # Ensure the level does not already exist because the level no can only
        # be added once. This might happen if the logger is re-initialized.
        try:
            logger.level("VERBOSE", color="<dim>")
        except ValueError:
            logger.level("VERBOSE", color="<dim>", no=15)
        logger.level("INFO", color="")
        logger.level("SUCCESS", color="<green>")
        logger.level("WARNING", color="<yellow>")
        logger.level("ERROR", color="<red>")
        # Bind the verbose level to the class so it can be used in any file even
        # after changing the logger instance
        logger.__class__.verbose = partialmethod(logger.__class__.log, "VERBOSE")

        self.reset_handlers(num_shots=None)

        # header
        package, *_ = __name__.split(".")
        commit = get_commit_hash()
        logger.info(
            "Starting: {p} ~ v{v}{t}{c} / {u}@{h}",
            p=package,
            v=importlib.metadata.version(package),
            t=" # " if commit else "",
            c=commit,
            u=os.getenv("USER"),
            h=os.uname().nodename,
        )
        if self.log_file_path is not None:
            logger.info("Logging: {l}", l=self.log_file_path)
        logger.debug(
            "Repository: {url}{append}{commit}",
            url="https://github.com/MIT-PSFC/disruption-py",
            append="/commit/" if commit else "",
            commit=commit,
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
        return LogSettings(console_log_level=log_settings)

    if log_settings is None:
        return LogSettings()

    raise ValueError(f"Invalid log settings {log_settings}")
