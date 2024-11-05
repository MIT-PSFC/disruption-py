#!/usr/bin/env python3

"""
This module defines the LogSettings class, which provides settings and setup for
logging in both files and console with customizable levels and formats.
"""

import importlib.metadata
import os
import sys
from dataclasses import dataclass

from loguru import logger
from tqdm.auto import tqdm

from disruption_py.core.utils.misc import get_commit_hash, get_temporary_folder


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
        Possible values are: "TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR",
        "CRITICAL". See https://loguru.readthedocs.io/en/stable/api/logger.html#levels
    log_file_write_mode : str, optional
        The write mode for the log file. Default is "w".
    log_to_console : bool
        Whether to log messages to the console (default is True).
    console_log_level : str
        The log level for the console. Default is None, so log level will be determined
        dynamically based on the number of shots.
        Possible values are: "TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR",
        "CRITICAL". See https://loguru.readthedocs.io/en/stable/api/logger.html#levels
    use_custom_logging : bool
        Whether to use custom logging. If set to true, no logging setup will be done.
        Default is False.
    warning_threshold : int
        If number of shots is greater than this threshold, the console log level will
        be "WARNING". Default is 500.
    success_threshold : int
        If number of shots is greater than this threshold and less than the warning_threshold,
        the console log level will be "SUCCESS". Default is 100.
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

        # Set colors without any bolding for each level (levels are bold by default)
        logger.level("DEBUG", color="<dim>")
        logger.level("INFO", color="")
        logger.level("SUCCESS", color="<green>")
        logger.level("WARNING", color="<yellow>")
        logger.level("ERROR", color="<red>")

        # formats
        message_format = "<level>[{level:^7s}] {message}</level>"
        console_format = "{time:HH:mm:ss.SSS} " + message_format
        file_format = "{time:YYYY-MM-DD HH:mm:ss.SSS} " + message_format

        if self.console_log_level is None:
            # Determine console log level dynamically based on the number of shots
            console_level = "INFO"
            if num_shots and num_shots > self.warning_threshold:
                console_level = "WARNING"
            elif num_shots and num_shots > self.success_threshold:
                console_level = "SUCCESS"
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
        Set up logging based on the provided settings.

        Parameters
        ----------
        logger_name : str, optional
            Name of the logger (default is "disruption_py").

        Returns
        -------
        logging.Logger
            Configured logger instance.
        """
        if self.use_custom_logging or self._logging_has_been_setup:
            return

        self.reset_handlers(num_shots=None)

        # header
        package = "disruption_py"
        commit = get_commit_hash()
        logger.info(
            "Starting: {p} ~ v{v} # {c} / {u}@{h}",
            p=package,
            v=importlib.metadata.version(package),
            c=commit,
            u=os.getenv("USER"),
            h=os.uname().nodename,
        )
        if self.log_file_path is not None:
            logger.info("Logging: {l}", l=self.log_file_path)
        logger.debug(
            "Repository: {r}/commit/{c}",
            r="https://github.com/MIT-PSFC/disruption-py",
            c=commit,
        )
        logger.debug("Executable: {e}", e=sys.executable)

        self._logging_has_been_setup = True
