#!/usr/bin/env python3

"""
example module for Two-Color Interferometry on C-Mod
"""

import pytest

from disruption_py.machine.tokamak import Tokamak, resolve_tokamak_from_environment
from disruption_py.settings import RetrievalSettings
from disruption_py.workflow import get_shots_data


def main():
    """
    execute a simple workflow to fetch TCI data.
    """

    tokamak = resolve_tokamak_from_environment()

    run_methods = ["get_ts_tci_comparison"]
    if tokamak is Tokamak.CMOD:
        shotlist = [1150805012]
    else:
        raise ValueError(f"Unspecified or unsupported tokamak: {tokamak}.")

    print(f"Initialized for tokamak: {tokamak.value}")

    retrieval_settings = RetrievalSettings(
        run_methods=run_methods,
        efit_nickname_setting="default",
    )

    result = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        output_setting="dataset",
    )

    print(result)

    var_is_all_empty = [result[var].isnull().all() for var in result.data_vars]
    assert not all(var_is_all_empty), "All variables are null, expected some non-null values."


if __name__ == "__main__":
    main()
