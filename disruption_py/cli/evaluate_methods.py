#!/usr/bin/env python3

import argparse
import logging

from disruption_py.settings.shotlist_setting import (
    ShotlistSettingParams,
    shotlist_setting_runner,
)
from disruption_py.utils.eval.eval_against_sql import (
    eval_against_sql,
    get_failure_statistics_string,
)
from disruption_py.core.utils.enums import map_string_to_enum
from disruption_py.machine.tokamak import Tokamak, get_tokamak_from_environment
from tests.utils.factory import (
    get_tokamak_test_shotlist,
)
from tests.utils.factory import get_tokamak_test_expected_failure_columns


def evaluate_accuracy(
    tokamak: Tokamak,
    shotlist: list[int],
    fail_quick: bool = False,
    data_columns: list[str] = None,
):
    if shotlist is None or len(shotlist) == 0:
        shotlist = get_tokamak_test_shotlist(tokamak)
    else:
        shotlist = [int(shot_id) for shot_id in shotlist]

    expected_failure_columns = get_tokamak_test_expected_failure_columns(tokamak)

    data_differences = eval_against_sql(
        tokamak=tokamak,
        shotlist=shotlist,
        expected_failure_columns=expected_failure_columns,
        fail_quick=fail_quick,
        test_columns=data_columns,
    )

    print(get_failure_statistics_string(data_differences))


def main(args):
    """
    Test the accuracy of DisruptionPy on a Tokamak

    Prints a short report on the methods that have suceeded and failed.
    Success criteria is having more that 95% of results within 1% of known results.
    """
    print(f"Welcome to the test accuracy script for disruption_py")

    tokamak_string: str = input(
        'Which tokamak would you like to evaluate methods for? ("cmod" for cmod, "d3d" for d3d, leave blank to auto-detect): '
    )
    while True:
        if tokamak_string == "":
            tokamak = get_tokamak_from_environment()
        else:
            tokamak = map_string_to_enum(tokamak_string, Tokamak, should_raise=False)

        if tokamak is None:
            tokamak_string = input(
                'Invalid input, please enter "cmod" for cmod, "d3d" for d3d, or leave blank to auto-detect: '
            )
        else:
            break

    # Get shot ids
    if args.shotlist is None:
        print()
        print(
            "What shot ids would you like to test? (Please enter shot ids as comma seperated list and/or on seperate lines)"
        )
        print(
            "To finish please leave a blank line (if no shots are entered a short default shot list will be used)"
        )
        all_shotlist = []
        while True:
            shotlist: str = input()
            if shotlist == "":
                break

            shotlist_to_add = [feature.strip() for feature in shotlist.split(",")]
            for shot_id in shotlist_to_add:
                try:
                    all_shotlist.append(int(shot_id))
                except ValueError:
                    print(
                        f"Shot id {shot_id} does not appear to be an integer, skipping"
                    )
    else:
        shotlist_setting_params = ShotlistSettingParams(
            database=None,
            tokamak=tokamak,
            logger=logging.getLogger("test_accuracy_logger"),
        )
        all_shotlist = shotlist_setting_runner(args.shotlist, shotlist_setting_params)

    data_columns = [args.data_column] if args.data_column else None

    print("Running evaluation...")
    evaluate_accuracy(
        tokamak=tokamak,
        shotlist=all_shotlist,
        fail_quick=args.fail_quick,
        data_columns=data_columns,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Evaluate the accuracy of DisruptionPy methods on a Tokamak."
    )
    parser.add_argument(
        "--shotlist",
        type=str,
        help="Path to file specifying a shotlist, leave blank for interactive mode",
        default=None,
    )
    parser.add_argument(
        "--fail_quick", action="store_true", help="Fail quickly", default=False
    )
    parser.add_argument(
        "--data_column", type=str, help="Data column to test", default=None
    )
    return parser
