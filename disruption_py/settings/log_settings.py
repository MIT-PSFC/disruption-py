import logging
from dataclasses import dataclass


@dataclass
class LogSettings:
    """
    Settings for logging.

    Attributes
    ----------
    log_file_path : str
        The path of the log file. Set to None if no log file should be created. Default is None.
    file_log_level : int
        The log level for the log file. Default is logging.INFO.
        Possible values are: logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL.
        See https://docs.python.org/3/library/logging.html#levels for more information.
    log_file_write_mode : str
        The write mode for the log file. Default is "w".
    log_to_console : bool
        Whether to log to the console. Default is True.
    console_log_level : int
        The log level for the console. Default is logging.INFO.
        Possible values are: logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL.
        See https://docs.python.org/3/library/logging.html#levels for more information.
    use_custom_logging : bool
        Whether to use custom logging. If set to true, no logging setup will be done. Default is False.
    """
    log_file_path: str = None
    file_log_level: int = logging.INFO
    log_file_write_mode: str = "w"
    
    log_to_console: bool = True
    console_log_level: int = logging.INFO
    
    use_custom_logging: bool = False
    
    _logging_has_been_setup: bool = False

    def default():
        return LogSettings()

    def setup_logging(self, logger_name = 'disruption_py'):
        if self.use_custom_logging or self._logging_has_been_setup:
            return
        self._logging_has_been_setup = True

        formatter = logging.Formatter('%(asctime)s,%(msecs)d %(name)s %(levelname)s | %(message)s', '%H:%M:%S')
        
        logger = logging.getLogger(logger_name)
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
            logger.addHandler(sh)