#!/usr/bin/env python3

"""
This module defines the LogSettings class, which provides settings and setup for 
logging in both files and console with customizable levels and formats.
"""

import logging
from dataclasses import dataclass


@dataclass
class LogSettings:
    """
    Settings for configuring logging.

    Attributes
    ----------
    log_file_path : str, optional
        Path to the log file. If None, no log file will be created (default is None).
    file_log_level : int, optional
        Logging level for the log file (default is logging.WARNING).
        Can be set to logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR,
        or logging.CRITICAL.
    log_file_write_mode : str, optional
        File write mode for the log file (default is 'w').
    log_to_console : bool, optional
        Whether to log messages to the console (default is True).
    console_log_level : int, optional
        Logging level for console output (default is logging.WARNING).
        Can be set to logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR,
        or logging.CRITICAL.
    use_custom_logging : bool, optional
        Whether to use custom logging setup (default is False).
        If True, no logging setup is done within this class.
    _logging_has_been_setup : bool, optional
        Internal flag to prevent multiple setups (default is False).
    """

    log_file_path: str = None
    file_log_level: int = logging.WARNING
    log_file_write_mode: str = "w"

    log_to_console: bool = True
    console_log_level: int = logging.WARNING

    use_custom_logging: bool = False

    _logging_has_been_setup: bool = False

    def logger(self, logger_name="disruption_py") -> logging.Logger:
        """
        Get or set up the logger with the specified settings.

        Parameters
        ----------
        logger_name : str, optional
            Name of the logger (default is "disruption_py").

        Returns
        -------
        logging.Logger
            The configured logger instance.
        """
        return self.setup_logging(logger_name)

    def setup_logging(self, logger_name="disruption_py") -> logging.Logger:
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
        logger = logging.getLogger(logger_name)
        if self.use_custom_logging or self._logging_has_been_setup:
            return logger
        self._logging_has_been_setup = True

        formatter = logging.Formatter(
            "%(asctime)s,%(msecs)d %(name)s %(levelname)s | %(message)s", "%H:%M:%S"
        )
        logger.propagate = False
        logger.handlers.clear()
        logger.setLevel(logging.DEBUG)
        if self.log_file_path is not None:
            fh = logging.FileHandler(self.log_file_path, mode=self.log_file_write_mode)
            fh.setFormatter(formatter)
            fh.setLevel(self.file_log_level)
            logger.addHandler(fh)
        if self.log_to_console:
            sh = logging.StreamHandler()
            sh.setLevel(self.console_log_level)
            sh.setFormatter(formatter)
            logger.addHandler(sh)
        return logger
