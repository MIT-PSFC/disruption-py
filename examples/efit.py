#!/usr/bin/env python3

"""
example module for EFIT.
"""

from disruption_py.machine.tokamak import Tokamak, resolve_tokamak_from_environment
from disruption_py.settings import RetrievalSettings
from disruption_py.workflow import get_shots_data


def main():
    """
    execute a simple workflow to fetch EFIT parameters.
    """

    tokamak = resolve_tokamak_from_environment()

    run_methods = ["get_efit_parameters"]
    if tokamak is Tokamak.D3D:
        shotlist = [161228]
        shape = (247, 16)
    elif tokamak is Tokamak.CMOD:
        shotlist = [1150805012]
        shape = (62, 22)
    elif tokamak is Tokamak.EAST:
        shotlist = [55555]
        shape = (69, 16)
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

    assert len(result), len(result.time) == shape


if __name__ == "__main__":
    main()
