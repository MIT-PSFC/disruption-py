#!/usr/bin/env python3

"""
Simple example to retrieve CER CERAUTO data for a single shot and save to netCDF.
"""

from disruption_py.machine.tokamak import resolve_tokamak_from_environment, Tokamak
from disruption_py.settings import RetrievalSettings
from disruption_py.workflow import get_shots_data


def main():
    tokamak = resolve_tokamak_from_environment()

    run_methods = ["get_all_cer_data"]
    if tokamak is Tokamak.D3D:
        shotlist = [144948]
    elif tokamak is Tokamak.CMOD:
        shotlist = [1120906030]
    else:
        raise ValueError(f"Unsupported tokamak: {tokamak}")

    retrieval_settings = RetrievalSettings(
        run_methods=run_methods,
        efit_nickname_setting="efit02",
        time_setting="mirnov",
    )

    for shot in shotlist:
        print(f"Processing shot: {shot}")
        result = get_shots_data(
            tokamak=tokamak,
            shotlist_setting=[shot],
            retrieval_settings=retrieval_settings,
            output_setting="dataset",
        )
        if result is None:
            print(f"No CER data retrieved for shot {shot}")
            continue

        outfn = f"./cer_{shot}.nc"
        result.to_netcdf(outfn, format="NETCDF4", mode="w")
        print(f"Saved: {outfn}")


if __name__ == "__main__":
    main()
