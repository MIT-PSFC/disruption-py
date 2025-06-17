#!/usr/bin/env python3

"""
example module for Mirnov FFTs on C-Mod.
"""

from disruption_py.machine.tokamak import Tokamak, resolve_tokamak_from_environment
from disruption_py.settings import RetrievalSettings
from disruption_py.workflow import get_shots_data


def main():
    """
    execute a simple workflow to fetch EFIT parameters.
    """

    tokamak = resolve_tokamak_from_environment()

    run_methods = ["get_all_mirnov_ffts"]
    if tokamak in [Tokamak.D3D, Tokamak.EAST]:
        raise ValueError(
            f"Mirnov FFTs are not supported for {tokamak.value}. "
            "Please use a different tokamak or method."
        )
    elif tokamak is Tokamak.CMOD:
        shotlist = [1110316031, 1160714026]
    else:
        raise ValueError(f"Unspecified or unsupported tokamak: {tokamak}.")

    print(f"Initialized for tokamak: {tokamak.value}")

    retrieval_settings = RetrievalSettings(
        run_methods=run_methods,
        efit_nickname_setting="default",
        time_setting="mirnov",
    )

    result = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        output_setting="dataset",
    )

    print(result)


if __name__ == "__main__":
    main()
