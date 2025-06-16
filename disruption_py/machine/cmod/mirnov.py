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

CMOD_R0 = 0.68  # Major radius of C-Mod [m]

def setup_stfft(f_mirnov=2.5e6, f_target=1e4, frequency_resolution=None) -> ShortTimeFFT:
    """Set up the Short Time Fast Fourier Transform object for the mirnov signals

    Parameters
    ----------
    f_mirnov : float, optional (default=2.5e6)
        The sampling frequency of the mirnov signals
    f_target : float, optional (default=1e4)
        The sampling frequency of the target timebase. 10 kHz is good for labeling.
    frequency_resolution : float, optional (default=None)
        The width of the frequency bins in the FFT. If None, will be set based on the Mirnov frequency and target frequency.

    Returns
    -------
    sft : scipy.signal.ShortTimeFFT
    """

    if frequency_resolution is None:
        frequency_resolution = int(f_mirnov / f_target)

    # Window size is the number of samples in the window
    # Think of how many samples you need to measure the minimum frequency you care about
    window_size = int(f_mirnov / frequency_resolution)
    # Hop is how much the window is slid between timesteps.
    # Slide the window based on the target frequency so interpolation is minimized
    hop = int(f_mirnov / f_target)
    # Something something kaiser window is optimal, just pick a beta that works
    window = kaiser(window_size, beta=14)
    sft = ShortTimeFFT(win=window, hop=hop, fs=f_mirnov)
    return sft

class CmodMirnovMethods:
    """
    Class for retrieving and processing Mirnov-related signals for CMOD.
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

        # These are the coils with positions that exist in MDSPlus
        mirnov_names_ab = [f"BP0{p}_ABK" for p in range(1, 10)] + [f"BP{p}_ABK" for p in range(10, 21)]
        mirnov_names_gh = [f"BP0{p}_GHK" for p in range(1, 10)] + [f"BP{p}_GHK" for p in range(10, 21)]
        mirnov_names_k = [f"BP0{p}_K" for p in range(1, 7)]

        # There's some tomfoolery in here because the probe locations are in a different tree from the signals, so we need to cut off the extra probes
        # that are NOT being digitized
        phi_ab, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.phi_AB", tree_name="magnetics")
        phi_gh, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.phi_GH", tree_name="magnetics")
        phi_k, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.phi_K", tree_name="magnetics")

        R_ab, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.R_AB", tree_name="magnetics")
        R_gh, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.R_GH", tree_name="magnetics")
        R_k, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.R_K", tree_name="magnetics")

        Z_ab, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.Z_AB", tree_name="magnetics")
        Z_gh, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.Z_GH", tree_name="magnetics")
        Z_k, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.Z_K", tree_name="magnetics")

        theta_pol_ab, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.theta_pol_AB", tree_name="magnetics")
        theta_pol_gh, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.theta_pol_GH", tree_name="magnetics")
        # theta_pol_k, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.theta_pol_K", tree_name="magnetics")  # noqa: ERA001
        theta_pol_k = np.empty_like(Z_k)  # Calibration data doesn't exist, so we'll just say NaN and deal with it later

        # For each of the above, we need to cut off the extra probes that are NOT being digitized
        # This is a bit of a kludge, but I'm expecting probes which aren't digitized to not show up in the analysis tree
        # TODO(ZanderKeith): Just hardcode all the analysis tree probes
        # The digitized probes are the first 20 for AB and GH, and the first 6 for K
        phi_ab = phi_ab[:len(mirnov_names_ab)]
        phi_gh = phi_gh[:len(mirnov_names_gh)]
        phi_k = phi_k[:7]

        R_ab = R_ab[:len(mirnov_names_ab)]
        R_gh = R_gh[:len(mirnov_names_gh)]
        R_k = R_k[:7]

        Z_ab = Z_ab[:len(mirnov_names_ab)]
        Z_gh = Z_gh[:len(mirnov_names_gh)]
        Z_k = Z_k[:7]

        theta_pol_ab = theta_pol_ab[:len(mirnov_names_ab)]
        theta_pol_gh = theta_pol_gh[:len(mirnov_names_gh)]
        theta_pol_k = theta_pol_k[:7]

        # These are the coils without positions in MDSPlus, need to hard-code values.
        # Taken from here: https://cmodwiki.psfc.mit.edu/index.php/FastMagneticsLocations#2010_Locations
        # Might need to include some logic for pre and post 2010... TODO(ZanderKeith)

        mirnov_names_tab = ["BP1T_ABK", "BP2T_ABK", "BP3T_ABK", "BP4T_ABK", "BP5T_ABK", "BP6T_ABK"]
        mirnov_names_tgh = ["BP1T_GHK", "BP2T_GHK", "BP3T_GHK", "BP4T_GHK", "BP5T_GHK", "BP6T_GHK"]
        mirnov_names_top = ["BP_KA_TOP", "BP_AB_TOP", "BP_BC_TOP", "BP_EF_TOP"]
        mirnov_names_bot = ["BP_KA_BOT", "BP_BC_BOT", "BP_EF_BOT"]

        phi_tab = [-23.10, -25.50, -27.90, -23.10, -25.50, -27.90]
        phi_tgh = [-224.40, -226.80, -229.20, -224.40, -226.80, -229.20]
        phi_top = [-344.80, -10.16, -59.87, -169.55]
        phi_bot = [-344.80, -59.87, -169.55]

        # Yes, the order of these is different from the c-mod wiki website. I'm putting it like phi, R, Z, theta_pol
        R_tab = [0.9045, 0.9045, 0.9045, 0.9045, 0.9045, 0.9045]
        R_tgh = [0.9042, 0.9042, 0.9042, 0.9042, 0.9042, 0.9042]
        R_top = [0.9126, 0.9126, 0.9146, 0.9126]
        R_bot = [0.9131, 0.9151, 0.9131]

        Z_tab = [0.1030, 0.1030, 0.1030, -0.1030, -0.1030, -0.1030]
        Z_tgh = [0.1000, 0.1000, 0.1000, -0.1000, -0.1000, -0.1000]
        Z_top = [0.0985, 0.0985, 0.0985, 0.1082]
        Z_bot = [-0.0985, -0.0985, -0.1092]

        theta_pol_tab = [23.7, 23.7, 23.7, -23.7, -23.7, -23.7]
        theta_pol_tgh = [23.1, 23.1, 23.1, -23.1, -23.1, -23.1]
        theta_pol_top = [18.0, 18.0, 18.0, 19.8]
        theta_pol_bot = [-18.0, -18.0, -20.0]

        all_mirnov_names = mirnov_names_ab + mirnov_names_gh + mirnov_names_k + mirnov_names_tab + mirnov_names_tgh + mirnov_names_top + mirnov_names_bot

        phi_all = np.concatenate((phi_ab, phi_gh, phi_k, phi_tab, phi_tgh, phi_top, phi_bot))
        phi_all = np.deg2rad(phi_all)
        R_all = np.concatenate((R_ab, R_gh, R_k, R_tab, R_tgh, R_top, R_bot))
        Z_all = np.concatenate((Z_ab, Z_gh, Z_k, Z_tab, Z_tgh, Z_top, Z_bot))

        # The angle we care about is when the hypotenuse is the minor radius
        Rd = R_all - CMOD_R0
        theta_all = np.arctan2(Z_all, Rd)

        theta_pol_all = np.concatenate((theta_pol_ab, theta_pol_gh, theta_pol_k, theta_pol_tab, theta_pol_tgh, theta_pol_top, theta_pol_bot))
        theta_pol_all = np.deg2rad(theta_pol_all)

        if debug:
            return all_mirnov_names[20:24], phi_all[20:24], theta_all[20:24], theta_pol_all[20:24]
        else:
            return all_mirnov_names, phi_all, theta_all, theta_pol_all

    @staticmethod
    def get_mirnov_stfft(params: PhysicsMethodParams, mirnov_name: str, freq_resolution: float = 250, max_freq: float = 80e3):
        """Get and the interpolated stfft of a Mirnov coil.
        
        This will work best if the params timebase is uniform.
        """
        path = r"\magnetics::top.active_mhd.signals"
        try:
            mirnov_signal, mirnov_times = params.mds_conn.get_data_with_dims(
                        path=f"{path}.{mirnov_name}",
                        tree_name="magnetics",
            )
            # Get the sampling frequency of the Mirnov signal
            f_mirnov = 1 / np.mean(np.diff(mirnov_times))
            params.logger.verbose(f"Using Mirnov frequency of {f_mirnov} Hz on {mirnov_name}")
            if not np.isclose(f_mirnov, 2.5e6, rtol=0.01) and not np.isclose(f_mirnov, 5e6, rtol=0.01):
                params.logger.warning(f"[Shot {params.shot_id}] Got Mirnov frequency of {f_mirnov} Hz on {mirnov_name}, expected 2.5 MHz or 5 MHz")

            # Set up the fft taking into account the params timebase
            f_timebase = 1 / np.mean(np.diff(params.times))
            sft = setup_stfft(f_mirnov=f_mirnov, f_target=f_timebase, frequency_resolution=freq_resolution)
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
        tokamak=Tokamak.CMOD,
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

        all_mirnov_names, phi_all, theta_all, theta_pol_all = CmodMirnovMethods.get_mirnov_names_and_locations(params)

        valid_mirnov_ffts = []
        valid_mirnov_names = []
        valid_mirnov_locations = []
        saved_freqs = None

        for mirnov_name, mirnov_phi, mirnov_theta, mirnov_theta_pol in zip(all_mirnov_names, phi_all, theta_all, theta_pol_all):
            mirnov_fft, freqs = CmodMirnovMethods.get_mirnov_stfft(params, mirnov_name)
            if mirnov_fft is not None:
                valid_mirnov_ffts.append(mirnov_fft)
                valid_mirnov_names.append(mirnov_name)
                valid_mirnov_locations.append((mirnov_phi, mirnov_theta, mirnov_theta_pol))

            if saved_freqs is None:
                saved_freqs = freqs

        mirnov_ffts_real = xr.DataArray(
            np.array(valid_mirnov_ffts).real,
            dims=("probe", "frequency", "time"),
            coords={
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
            dims=("probe", "frequency", "time"),
            coords={
                "probe": list(range(len(valid_mirnov_locations))),
                "probe_name": ("probe", valid_mirnov_names),
                "frequency": saved_freqs,
                "time": params.times,
                "phi": ("probe", [loc[0] for loc in valid_mirnov_locations]),
                "theta": ("probe", [loc[1] for loc in valid_mirnov_locations]),
                "theta_pol": ("probe", [loc[2] for loc in valid_mirnov_locations]),
            },
        )
        mirnov_ds_real = mirnov_ffts_real.to_dataset(name="mirnov_fft_real", promote_attrs=True)
        mirnov_ds_imag = mirnov_ffts_imag.to_dataset(name="mirnov_fft_imag", promote_attrs=True)
        mirnov_ds = xr.merge([mirnov_ds_real, mirnov_ds_imag])

        return mirnov_ds