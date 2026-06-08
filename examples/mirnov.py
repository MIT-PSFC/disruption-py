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

    run_methods = ["get_all_mirnov_ffts", "get_geqdsk_parameters",\
                   "get_ip_parameters", "get_efit_parameters", "efit_check"]
    if tokamak in [Tokamak.D3D]:
        shotlist = [179118]#161228
    elif tokamak is Tokamak.CMOD:
        shotlist = [1120906030
    ] #1160714026,1160930034 1160826001 1051202011  1160714026 1110316031
    else:
        raise ValueError(f"Unspecified or unsupported tokamak: {tokamak}.")

    print(f"Initialized for tokamak: {tokamak.value}")

    retrieval_settings = RetrievalSettings(
        run_methods=run_methods,
        efit_nickname_setting="default",
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

        # Apply frequency response correction to the Mirnov FFT data
        result = apply_freq_correction(result)
        
        
        # result.to_netcdf(f'../data_archive/{shot}.nc',format='NETCDF4', mode='w')
        result.to_netcdf(f'../TARS/tars/scratch/input_data/{shot}.nc',format='NETCDF4',mode='w')


    print(result)
    

if __name__ == "__main__":
    main()

'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

vmin, vmax = 0.0, 6  # set your desired limits
levels = np.linspace(vmin, vmax, 101)
norm = mcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)

fig, ax = plt.subplots(2, 1, layout='constrained', sharex=True, figsize=(4, 5))
dat = np.sqrt( result.sel(probe=0).mirnov_fft_real**2 + result.sel(probe=0).mirnov_fft_real**2 ) 
im = ax[0].contourf(
	result.time,
	result.frequency * 1e-3,
	dat,
	levels=levels,      # critical: cap contour levels at vmax
	norm=norm,
	extend='neither'
)
ax[0].set_ylabel('Frequency [kHz]')

cbar = fig.colorbar(im, ax=ax[0], label='Mirnov Signal [T/s?]', extend='neither')
cbar.set_ticks(np.linspace(vmin, vmax, 6))  # optional, for clean tick control

ax[1].plot(result.time, result.ip * 1e-6, label=f'MPI66M020 Shot {shotlist[0]}')
ax[1].legend()
ax[1].set_ylabel('$I_p$ [MA]')
ax[1].set_xlabel('Time [s]')
ax[1].grid()

fig.savefig('/home/chandrar/Mirnov_out_1.png')
plt.show()

'''