#!/usr/bin/env python3

"""
Module for retrieving and processing EFIT parameters for CMOD.
"""

import numpy as np
from MDSplus import mdsExceptions
from scipy.signal import ShortTimeFFT
from scipy.signal.windows import kaiser
import xarray as xr

from disruption_py.core.physics_method.decorator import physics_method, cache_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
from disruption_py.machine.tokamak import Tokamak
from disruption_py.machine.cmod.mirnov import setup_fft

class D3DMirnovMethods:
    """
    Class for retrieving and processing Mirnov-related signals for D3D
    """

    @staticmethod
    def get_mirnov_names_and_locations(params: PhysicsMethodParams, debug=False):
        """Get the names and locations of the Mirnov coils in the ANALYSIS tree.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters for the physics method.

        Returns
        -------
        all_mirnov_names : list[str]
            The names of all the Mirnov coils.
        phi_all : numpy.ndarray
            The toroidal angles of the Mirnov coils
        theta_all : numpy.ndarray
            The poloidal angles of the Mirnov coils
        theta_pol_all : numpy.ndarray
            The poloidal offset angles of the Mirnov coils.
        """
        
    @staticmethod
    def get_mirnov_fft(params: PhysicsMethodParams, mirnov_name: str, freq_resolution: float = 100, max_freq: float = 80e3):
        """Get and the interpolated fft of a Mirnov coil.
        
        This will work best if the params timebase is uniform.
        """
        path = r"\magnetics::top.active_mhd.signals"
        try:
            mirnov_signal, mirnov_times = params.mds_conn.get_data_with_dims(
                        path=f"{path}.{mirnov_name}",
                        tree_name="magnetics",
                        astype="float32",
            )
            # Get the sampling frequency of the Mirnov signal
            f_mirnov = 1 / np.mean(np.diff(mirnov_times))
            params.logger.verbose(f"Using Mirnov frequency of {f_mirnov} Hz on {mirnov_name}")
            if not np.isclose(f_mirnov, 2.5e6, rtol=0.01) and not np.isclose(f_mirnov, 5e6, rtol=0.01):
                params.logger.warning(f"[Shot {params.shot_id}] Got Mirnov frequency of {f_mirnov} Hz on {mirnov_name}, expected 2.5 MHz or 5 MHz")

            # Set up the fft taking into account the params timebase
            f_timebase = 1 / np.mean(np.diff(params.times))
            sft = setup_fft(f_mirnov=f_mirnov, f_target=f_timebase, frequency_resolution=freq_resolution)
            f_indices = np.where(sft.f < max_freq)
            freqs = sft.f[f_indices]

            mirnov_fft_full = sft.stft(mirnov_signal)[f_indices]

            # Interpolate the STFFTs to the timebase
            fft_times = (sft.delta_t * np.arange(mirnov_fft_full.shape[1])) + mirnov_times[0]
            mirnov_fft_interp = interp1(fft_times, mirnov_fft_full, params.times)

            # Check if the average difference between frequencies is NOT close to the frequency resolution
            if not np.isclose(np.mean(np.diff(freqs)), freq_resolution, atol=1):
                params.logger.warning(f"[Shot {params.shot_id}] The Mirnov frequency resolution is not {freq_resolution} Hz")
            # Replace the frequencies with nice integers (for the sake of consistency)
            freqs = np.arange(0, max_freq, freq_resolution)

            return mirnov_fft_interp, freqs
        except Exception as e:
            return None, None

    @staticmethod
    @physics_method(
        tokamak=Tokamak.D3D,
    )
    def get_all_mirnov_ffts(params: PhysicsMethodParams):
        """Get all FFTs of the available Mirnov coils for this shot.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters for the physics method.

        Returns
        -------
        mirnov_ffts : xarray.DataSet
            The FFT of the Mirnov coils.
            Dimensions are probe, frequency, and time.
            Coordinates are probe, frequency, time, phi, theta, and theta_pol.
        """

        all_mirnov_names, phi_all, theta_all, theta_pol_all = D3DMirnovMethods.get_mirnov_names_and_locations(params, debug=True)

        valid_mirnov_ffts = []
        valid_mirnov_names = []
        valid_mirnov_locations = []
        saved_freqs = None

        for mirnov_name, mirnov_phi, mirnov_theta, mirnov_theta_pol in zip(all_mirnov_names, phi_all, theta_all, theta_pol_all):
            mirnov_fft, freqs = D3DMirnovMethods.get_mirnov_fft(params, mirnov_name)
            if mirnov_fft is not None:
                valid_mirnov_ffts.append(mirnov_fft)
                valid_mirnov_names.append(mirnov_name)
                valid_mirnov_locations.append((mirnov_phi, mirnov_theta, mirnov_theta_pol))

            if saved_freqs is None:
                saved_freqs = freqs

        valid_mirnov_ffts = np.expand_dims(valid_mirnov_ffts, axis=0)  # Add a new axis for the probe dimension

        mirnov_ffts_real = xr.DataArray(
            np.array(valid_mirnov_ffts).real,
            dims=("idx", "probe", "frequency", "time"),
            coords={
                "shot": ("idx", [params.shot_id]),
                "probe": list(range(len(valid_mirnov_locations))),
                "probe_name": ("probe", valid_mirnov_names),
                "frequency": saved_freqs,
                "time": params.times,
                "phi": ("probe", [loc[0] for loc in valid_mirnov_locations]),
                "theta": ("probe", [loc[1] for loc in valid_mirnov_locations]),
                "theta_pol": ("probe", [loc[2] for loc in valid_mirnov_locations]),
            },
        )
        mirnov_ffts_imag = xr.DataArray(
            np.array(valid_mirnov_ffts).imag,
            dims=("idx", "probe", "frequency", "time"),
            coords={
                "shot": ("idx", [params.shot_id]),
                "probe": list(range(len(valid_mirnov_locations))),
                "probe_name": ("probe", valid_mirnov_names),
                "frequency": saved_freqs,
                "time": params.times,
                "phi": ("probe", [loc[0] for loc in valid_mirnov_locations]),
                "theta": ("probe", [loc[1] for loc in valid_mirnov_locations]),
                "theta_pol": ("probe", [loc[2] for loc in valid_mirnov_locations]),
            },
        )
        mirnov_ds_real = mirnov_ffts_real.astype(np.float32).to_dataset(name="mirnov_fft_real")
        mirnov_ds_imag = mirnov_ffts_imag.astype(np.float32).to_dataset(name="mirnov_fft_imag")
        mirnov_ds = xr.merge([mirnov_ds_real, mirnov_ds_imag])

        return mirnov_ds
    
    @staticmethod
    @physics_method(
        tokamak=Tokamak.D3D,
    )
    def get_preferred_mirnov_sxx(params: PhysicsMethodParams):
        """Get the Sxx of a single Mirnov coil in the shot, in order of preference from a few consistently 'good' probes"""

        preferred_mirnov_names = ["BP02_GHK", "BP01_ABK"]

        for mirnov_name in preferred_mirnov_names:
            mirnov_fft, freqs = D3DMirnovMethods.get_mirnov_fft(params, mirnov_name)
            if mirnov_fft is not None:
                params.logger.info(f"Using {mirnov_name} for Sxx")
                mirnov_sxx = np.abs(mirnov_fft) ** 2
                return xr.DataArray(
                    mirnov_sxx,
                    dims=("frequency", "time"),
                    coords={
                        "shot": params.shot_id,
                        "frequency": freqs,
                        "time": params.times,
                    },
                    attrs={
                        "probe_name": mirnov_name,
                    },
                ).astype(np.float32).to_dataset(name="mirnov_sxx")