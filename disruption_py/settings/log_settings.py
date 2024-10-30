#!/usr/bin/env python3

"""
This module defines the LogSettings class, which provides settings and setup for 
logging in both files and console with customizable levels and formats.
"""

import os
import sys
from dataclasses import dataclass

from loguru import logger


@dataclass
class LogSettings:
    """
    Settings for configuring logging.

    Attributes
    ----------
    log_file_path : str, optional
        Path to the log file. If None, no log file will be created (default is None).
    file_log_level : str
        Logging level for the log file (default is "DEBUG").
        Possible values are: "TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "SUMMARY", "ERROR",
        "CRITICAL". See https://loguru.readthedocs.io/en/stable/api/logger.html#levels
    log_file_write_mode : str, optional
        The write mode for the log file. Default is "w".
    log_to_console : bool
        Whether to log messages to the console (default is True).
    console_log_level : str
        The log level for the console. Default is "WARNING".
        Possible values are: "TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "SUMMARY", "ERROR",
        "CRITICAL". See https://loguru.readthedocs.io/en/stable/api/logger.html#levels
    use_custom_logging : bool
        Whether to use custom logging. If set to true, no logging setup will be done.
        Default is False.
    _logging_has_been_setup : bool
        Internal flag to prevent multiple setups (default is False).
    """

    log_file_path: str = None
    file_log_level: str = "DEBUG"
    log_file_write_mode: str = "w"

    log_to_console: bool = True
    console_log_level: str = "WARNING"

    use_custom_logging: bool = False

    @property
    def _logging_has_been_setup(self) -> bool:
        """
        Return True if the logger has been set up, False otherwise.
        """
        try:
            # The logger is set up when the SUMMARY level exists.
            logger.level("SUMMARY")
            return True
        except ValueError:
            return False

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

        # Remove default logger
        logger.remove()

        # Add summary level after warning (30) but before error (40)
        logger.level("SUMMARY", no=35, color="<cyan>")

        # formats
        message_format = "<level>[{level:^7s}] {message}</level>"
        console_format = "{time:HH:mm:ss.SSS} " + message_format
        file_format = "{time:YYYY-MM-DD HH:mm:ss.SSS} " + message_format

        # Add console handler
        if self.log_to_console:
            logger.add(
                sys.stderr,
                level=self.console_log_level,
                format=console_format,
                colorize=True,
                enqueue=True,
                backtrace=False,
                diagnose=True,
            )

        # Add file handler if log file path is provided
        if self.log_file_path is not None:
            logger.debug("Logging to file {file_path}", file_path=self.log_file_path)
            logger.add(
                self.log_file_path,
                level=self.file_log_level,
                format=file_format,
                mode=self.log_file_write_mode,
                enqueue=True,
                backtrace=False,
                diagnose=True,
            )

        logger.debug(os.getenv("VIRTUAL_ENV"))
