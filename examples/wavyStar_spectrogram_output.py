#!/usr/bin/env python3

"""
example module for Mirnov FFTs on C-Mod.
"""

from disruption_py.machine.tokamak import Tokamak, resolve_tokamak_from_environment
from disruption_py.settings import RetrievalSettings
from disruption_py.workflow import get_shots_data
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt


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
        shotlist = [ 1160714026]
    else:
        raise ValueError(f"Unspecified or unsupported tokamak: {tokamak}.")

    print(f"Initialized for tokamak: {tokamak.vspect_xralue}")

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
    save_output(result)
                
#########################################################################333
def save_output(result, filename="wavyStar_mirnov_spectrograms.nc"):
    """
    Save the output to a Zarr file.
    """

    output_index = slice(23,25)
    out_SSX = result['mirnov_fft_real'][output_index,:]**2 + result['mirnov_fft_imag'][output_index,:]**2
    out_SSX = np.mean(out_SSX,axis=0).values
    out_SSX  = np.log10(out_SSX)

   
    spect_xr = xr.Dataset({  "Sxx": (['frequency','time'],out_SSX) },
            coords = {'frequency':result['frequency'],'time':result['time']},
            attrs = {'sampling_frequency':(1/(result['time'][1]-result['time'][0])).values}
        )
        # spect_xr = xr.DataArray(out_spect[f_inds].T,dims=['time','frequency'],\
        #                         coords={'time':time,'frequency':freq[f_inds]})
    spect_xr.to_netcdf('/home/rianc/Documents/Synthetic_Mirnov/data_output/Spectrogram_Xarrays/%s'%filename)
    print(f"Saved spectrogram to /home/rianc/Documents/Synthetic_Mirnov/data_output/Spectrogram_Xarrays/{filename}")


    # Plotting
    plt.close('all')
    fig,ax = plt.subplots(1,1,figsize=(10,3),tight_layout=True)

    ax.contourf(spect_xr['time'],spect_xr['frequency']*1e-3,spect_xr['Sxx'],vmin=4,vmax=9,levels=100,zorder=-5)
    ax.set_rasterization_zorder(-1)

    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Frequency [kHz]')
    
    # Add shot number in upper right corner
    ax.text(0.98, 0.98, f"Shot: {result.shot.values[0]}", 
            transform=ax.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))


    fig.savefig('/home/rianc/Documents/Synthetic_Mirnov/output_plots/Spectrogram_Xarrays_%d.pdf'%result.shot.values[0],\
                dpi=200, transparent=True)

    plt.show()
if __name__ == "__main__":
    main()
