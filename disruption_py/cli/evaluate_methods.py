#!/usr/bin/env python3

import argparse
import logging

from disruption_py.settings.shot_ids_request import (
    ShotIdsRequestParams,
    shot_ids_request_runner,
)
from disruption_py.utils.eval.eval_against_sql import (
    eval_against_sql,
    get_failure_statistics_string,
)
from disruption_py.utils.mappings.mappings_helpers import map_string_to_enum
from disruption_py.utils.mappings.tokamak import Tokamak, get_tokamak_from_environment
from disruption_py.utils.mappings.tokamak_helpers import (
    get_tokamak_test_expected_failure_columns,
    get_tokamak_test_shot_ids,
)


def evaluate_accuracy(
    tokamak: Tokamak,
    shot_ids: list[int],
    fail_quick: bool = False,
    data_columns: list[str] = None,
):
    if shot_ids is None or len(shot_ids) == 0:
        shot_ids = get_tokamak_test_shot_ids(tokamak)
    else:
        shot_ids = [int(shot_id) for shot_id in shot_ids]

    expected_failure_columns = get_tokamak_test_expected_failure_columns(tokamak)

    data_differences = eval_against_sql(
        tokamak=tokamak,
        shot_ids=shot_ids,
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
        all_shot_ids = []
        while True:
            shot_ids: str = input()
            if shot_ids == "":
                break

            shot_ids_to_add = [feature.strip() for feature in shot_ids.split(",")]
            for shot_id in shot_ids_to_add:
                try:
                    all_shot_ids.append(int(shot_id))
                except ValueError:
                    print(
                        f"Shot id {shot_id} does not appear to be an integer, skipping"
                    )
    else:
        shot_ids_request_params = ShotIdsRequestParams(
            database=None,
            tokamak=tokamak,
            logger=logging.getLogger("test_accuracy_logger"),
        )
        all_shot_ids = shot_ids_request_runner(args.shotlist, shot_ids_request_params)

    data_columns = [args.data_column] if args.data_column else None

    print("Running evaluation...")
    evaluate_accuracy(
        tokamak=tokamak,
        shot_ids=all_shot_ids,
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
