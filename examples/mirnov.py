#!/usr/bin/env python3

"""
example module for Mirnov FFTs on C-Mod.
"""

from disruption_py.machine.tokamak import Tokamak, resolve_tokamak_from_environment
from disruption_py.settings import RetrievalSettings
from disruption_py.workflow import get_shots_data

from mirnov_freq_correction import apply_freq_correction


def main():
    """
    execute a simple workflow to fetch EFIT parameters.
    """

    tokamak = resolve_tokamak_from_environment()

    run_methods = ["get_all_mirnov_ffts", "get_geqdsk_parameters"]
    if tokamak in [Tokamak.D3D, Tokamak.EAST]:
        raise ValueError(
            f"Mirnov FFTs are not supported for {tokamak.value}. "
            "Please use a different tokamak or method."
        )
    elif tokamak is Tokamak.CMOD:
        shotlist = [1160714026   ] #1160714026,1160930034 1160826001 1051202011  1160714026 1110316031
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

    # Apply frequency response correction to the Mirnov FFT data
    # result = apply_freq_correction(result)
    
    
    result.to_netcdf(f'../data_archive/{shotlist[0]}.nc',format='NETCDF4')
    result.to_netcdf(f'../TARS/tars/scratch/input_data/{shotlist[0]}.nc',format='NETCDF4',mode='w')


    print(result)
    

if __name__ == "__main__":
    main()
