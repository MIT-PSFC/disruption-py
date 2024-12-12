#!/usr/bin/env python3

"""
simple module with command-line argument parsing.
"""

import argparse

from disruption_py.machine.tokamak import resolve_tokamak_from_environment
from disruption_py.settings import RetrievalSettings
from disruption_py.workflow import get_shots_data
from tests.utils.factory import get_tokamak_test_shotlist


def main(tokamak, methods, shots, processes, log_level):
    """
    simple workflow.
    """

    if not tokamak:
        tokamak = resolve_tokamak_from_environment()
    if not shots:
        shots, *_ = get_tokamak_test_shotlist(tokamak)
    tags = [] if methods else ["all"]

    out = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shots,
        retrieval_settings=RetrievalSettings(run_methods=methods, run_tags=tags),
        num_processes=processes,
        log_settings=log_level,
    )
    return out


def cli():
    """
    simple argument parser.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--tokamak", type=str)
    parser.add_argument("-s", "--shots", type=int, action="append")
    parser.add_argument("-m", "--methods", type=str, action="append")
    parser.add_argument("-p", "--processes", type=int, default=1)
    parser.add_argument("-l", "--log-level")

    return main(**vars(parser.parse_args()))


if __name__ == "__main__":
    print(cli())
