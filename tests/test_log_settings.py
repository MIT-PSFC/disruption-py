#!/usr/bin/env python3

"""
This module contains tests for the logging setup across processes: worker
console lines are relayed through the main process, and every queued line
reaches the log file before the pool exits.
"""

import os
from multiprocessing import Pool

from loguru import logger

from disruption_py.settings.log_settings import ConsoleRelay, LogSettings

NUM_WORKERS = 2
NUM_LINES = 20


# Module-scope so spawn/forkserver workers can unpickle the function reference.
def _log_lines(worker):
    for i in range(NUM_LINES):
        logger.info("worker {} line {}", worker, i)
    # as in `_execute_retrieval`
    logger.complete()
    return worker


def _run_pool(log_settings):
    """
    Run a small pool wired like `get_shots_data`.
    """
    log_settings.reset_handlers()
    with (
        ConsoleRelay() as relay,
        Pool(
            processes=NUM_WORKERS,
            initializer=log_settings.reset_handlers,
            initargs=(relay.queue,),
        ) as pool,
    ):
        list(pool.imap(_log_lines, range(NUM_WORKERS)))
    logger.complete()
    logger.remove()


def _expected_lines():
    return [
        f"worker {worker} line {i}"
        for worker in range(NUM_WORKERS)
        for i in range(NUM_LINES)
    ]


def test_worker_console_lines_relayed(capsys, test_folder_f):
    """
    Every worker console line is written by the main process.
    """
    _run_pool(
        LogSettings(
            console_level="INFO",
            file_path=os.path.join(test_folder_f, "output.log"),
        )
    )
    out = capsys.readouterr().out
    missing = [line for line in _expected_lines() if line not in out]
    assert not missing


def test_worker_file_lines_complete(test_folder_f):
    """
    Every worker line reaches the log file, including the tail written just
    before the pool terminates its workers.
    """
    file_path = os.path.join(test_folder_f, "output.log")
    _run_pool(LogSettings(console_level="WARNING", file_path=file_path))
    with open(file_path, encoding="utf8") as f:
        content = f.read()
    missing = [line for line in _expected_lines() if line not in content]
    assert not missing
    assert content.count("worker ") == NUM_WORKERS * NUM_LINES
