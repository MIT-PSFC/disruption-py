#!/usr/bin/env python3

"""
simple module with command-line argument parsing.
"""

import argparse
import os

from loguru import logger

from disruption_py.core.utils.misc import get_temporary_folder
from disruption_py.machine.tokamak import resolve_tokamak_from_environment
from disruption_py.settings import RetrievalSettings
from disruption_py.workflow import get_shots_data
from tests.utils.factory import get_tokamak_test_shotlist


def main(
    tokamak, methods, shots, efit_tree, time_base, output_file, processes, log_level
):
    """
    simple workflow.
    """

    if not tokamak:
        tokamak = resolve_tokamak_from_environment()
    if not shots:
        shots, *_ = get_tokamak_test_shotlist(tokamak)
    methods = methods or None

    sett = RetrievalSettings(
        efit_nickname_setting=efit_tree, time_setting=time_base, run_methods=methods
    )
    out = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shots,
        retrieval_settings=sett,
        num_processes=processes,
        log_settings=log_level,
        output_setting="dataset",
    )

    if output_file:
        output_file = os.path.realpath(output_file)
    else:
        output_file = os.path.join(get_temporary_folder(), "output.nc")
    logger.info("Output: {output_file}", output_file=output_file)
    out.to_netcdf(output_file)

    return out


def cli():
    """
    simple argument parser.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("shots", nargs="*")
    parser.add_argument("-t", "--tokamak", type=str)
    parser.add_argument("-m", "--methods", type=str, action="append")
    parser.add_argument("-e", "--efit-tree", type=str, default="disruption")
    parser.add_argument("-b", "--time-base", default="disruption_warning")
    parser.add_argument("-o", "--output-file", type=str)
    parser.add_argument("-p", "--processes", type=int, default=1)
    parser.add_argument("-l", "--log-level", default="VERBOSE")

    return main(**vars(parser.parse_args()))


if __name__ == "__main__":
    ds = cli()
    print(ds)
