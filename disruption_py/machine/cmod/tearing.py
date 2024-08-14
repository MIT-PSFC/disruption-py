#!/usr/bin/env python3

import traceback
import warnings
from importlib import resources

import numpy as np
import pandas as pd
import xarray as xr
from MDSplus import mdsExceptions
from scipy.signal import ShortTimeFFT
from scipy.signal.windows import kaiser

import disruption_py.data
from disruption_py.core.physics_method.caching import cache_method
from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import gaussian_fit, interp1, smooth
from disruption_py.machine.tokamak import Tokamak


def setup_stfft(f_mirnov=2.5e6, f_target=1e3, frequency_resolution=None) -> ShortTimeFFT:
    """Set up the Short Time Fast Fourier Transform object for the mirnov signals

    Parameters
    ----------
    f_mirnov : float, optional (default=2.5e6)
        The sampling frequency of the mirnov signals
    f_target : float, optional (default=1e3)
        The sampling frequency of the target timebase
    frequency_resolution : float, optional (default=None)
        The width of the frequency bins in the FFT. If None, will be set to 250 Hz

    Returns
    -------
    SFT : scipy.signal.ShortTimeFFT
        The ShortTimeFFT object
    """
    # Frequency resolution is the minimum frequency of interest, the size of the bins
    if frequency_resolution is None:
        frequency_resolution = 250

    # Window size is the number of samples in the window
    # Think of how many samples you need to measure the minimum frequency you care about
    window_size = int(f_mirnov / frequency_resolution)
    # Hop is how much the window is slid between timesteps.
    # Slide the window based on the target frequency so interpolation is minimized
    hop = int(f_mirnov / f_target)
    # Something something kaiser window is optimal, just pick a beta that works
    window = kaiser(window_size, beta=14)
    SFT = ShortTimeFFT(win=window, hop=hop, fs=f_mirnov)
    return SFT


def get_expected_cross_phase(mode: tuple[int, int], theta_1: float, phi_1: float, theta_2: float, phi_2: float):
    """Calculate the expected cross phase between two mirnov probes if there were an m/n mode present
    This cross phase will match if it is calculated like so:
    cross_spectrum = probe_1_fft * np.conj(probe_2_fft)
    cross_phase = np.angle(cross_spectrum) % (2 * np.pi)

    Parameters:
    ----------
    mode : tuple[int, int]
        The m/n mode that is expected to be present
    theta_1 : float
        The poloidal angle of the first mirnov probe [rad]
    phi_1 : float
        The toroidal angle of the first mirnov probe [rad]
    theta_2 : float
        The poloidal angle of the second mirnov probe [rad]
    phi_2 : float
        The toroidal angle of the second mirnov probe [rad]

    Returns:
    ----------
    expected_cross_phase : float
        The expected cross phase between the two probes.
    """

    poloidal_phase = mode[0] * (theta_2 - theta_1)
    toroidal_phase = mode[1] * (phi_2 - phi_1)
    expected_cross_phase = (poloidal_phase + toroidal_phase) % (2 * np.pi)

    return expected_cross_phase


def cross_phase_mask(original_signal, cross_phases, expected_cross_phases, tolerance=0.5):
    """Mask the signal based on the cross phases between multiple probes

    Parameters:
    ----------
    original_signal : numpy.ndarray
        The original signal to mask
    cross_phases : list[float]
        The cross phases between the probes
    expected_cross_phases : list[float]
        The expected cross phases between the probes
    tolerance : float, optional (default=0.5) [rad]
        The tolerance for the cross phase to be considered a match

    Returns:
    ----------
    masked_signal : numpy.ndarray
        The original signal with the masked values set to np.nan
    """

    masked_signal = original_signal
    for cross_phase, expected_cross_phase in zip(cross_phases, expected_cross_phases):
        # If the expected cross phase is near 0 or 2pi, we risk a 2 pi jump messing things up
        # Gotta bring everything into the center of the circle
        if expected_cross_phase < np.pi / 8:
            tested_expected_cross_phase = expected_cross_phase + np.pi
            tested_cross_phase = (cross_phase + np.pi) % (2 * np.pi)
        elif expected_cross_phase > 15 * np.pi / 8:
            tested_expected_cross_phase = expected_cross_phase - np.pi
            tested_cross_phase = (cross_phase - np.pi) % (2 * np.pi)
        else:
            tested_expected_cross_phase = expected_cross_phase
            tested_cross_phase = cross_phase
        masked_signal = np.where(np.abs(tested_cross_phase - tested_expected_cross_phase) > tolerance, np.nan, masked_signal)
    return masked_signal


class CmodTearingMethods:
    @staticmethod
    @cache_method
    def get_three_mirnov_signals(params: PhysicsMethodParams):  # noqa: PLR0915
        """Get measurements from up to 3 available Mirnov coils.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters for the physics method.

        Returns
        -------
        mirnov_times : numpy.ndarray
            The times of the Mirnov signals.
        mirnov_signals: list[numpy.ndarray]
            The signals from the Mirnov coils.
        probe_locations: list[tuple[float, float, float]]
            The (phi, theta, theta_offset) locations of the probes in radians.
        """

        # Get all names of the Mirnov coils TODO: how do you do this in disruption-py?
        # I want the equivalent of .getNode, or .getNodeWild, and things like that
        # mirnov_coil_names = params.mds_connection.get_coil_names()
        # TODO: expand for every Mirnov coil. For now, just using the ABK, GHK, and K coils
        mirnov_names_ab = [f"BP0{p}_ABK" for p in range(1, 10)] + [f"BP{p}_ABK" for p in range(10, 19)]
        mirnov_names_gh = [f"BP0{p}_GHK" for p in range(1, 10)] + [f"BP{p}_GHK" for p in range(10, 19)]
        mirnov_names_k = [f"BP0{p}_K" for p in range(1, 7)]
        path = r"\magnetics::top.active_mhd.signals"

        phi_ab, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.phi_AB", tree_name="magnetics")
        phi_gh, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.phi_GH", tree_name="magnetics")
        phi_k, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.phi_K", tree_name="magnetics")
        phi_all = np.concatenate((phi_ab, phi_gh, phi_k))
        phi_all = np.deg2rad(phi_all)

        R_ab, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.R_AB", tree_name="magnetics")
        R_gh, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.R_GH", tree_name="magnetics")
        R_k, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.R_K", tree_name="magnetics")
        R_all = np.concatenate((R_ab, R_gh, R_k))

        Z_ab, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.Z_AB", tree_name="magnetics")
        Z_gh, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.Z_GH", tree_name="magnetics")
        Z_k, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.Z_K", tree_name="magnetics")
        Z_all = np.concatenate((Z_ab, Z_gh, Z_k))

        theta_all = np.arctan2(Z_all, R_all)

        theta_pol_ab, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.theta_pol_AB", tree_name="magnetics")
        theta_pol_gh, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.theta_pol_GH", tree_name="magnetics")
        # theta_pol_k, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.theta_pol_K", tree_name="magnetics")  # noqa: ERA001
        theta_pol_k = np.empty_like(Z_k)  # Calibration data doesn't exist, so we'll just say NaN and deal with it later
        theta_pol_all = np.concatenate((theta_pol_ab, theta_pol_gh, theta_pol_k))
        theta_pol_all = np.deg2rad(theta_pol_all)

        original_names = mirnov_names_ab + mirnov_names_gh + mirnov_names_k
        unchecked_names = original_names.copy()

        mirnov_times = None
        mirnov_signals = []
        mirnov_locations = []
        mirnov_probes = []
        while len(mirnov_signals) < 3 and len(unchecked_names) > 0:
            if len(mirnov_locations) > 0:
                # Pick an unchecked probe with the largest average distance from the presently selected probes
                all_distances = np.zeros((len(unchecked_names), len(mirnov_locations)))
                for i, probe in enumerate(unchecked_names):
                    probe_index = original_names.index(probe)
                    for j, loc in enumerate(mirnov_locations):
                        phi, theta, theta_pol = loc
                        all_distances[i, j] = np.sqrt((phi - phi_all[probe_index]) ** 2 + (theta - theta_all[probe_index]) ** 2)
                distances = np.mean(all_distances, axis=1)
                mirnov_name = unchecked_names[np.argmax(distances)]
            else:
                # Haven't picked any yet, just find the first probe in the unchecked list
                mirnov_name = unchecked_names[0]

            selected_index = original_names.index(mirnov_name)

            try:
                mirnov_signal, mirnov_time = params.mds_conn.get_data_with_dims(
                    path=f"{path}.{mirnov_name}",
                    tree_name="magnetics",
                )

                mirnov_signals.append(mirnov_signal)
                phi = phi_all[selected_index]
                theta = theta_all[selected_index]
                theta_pol = theta_pol_all[selected_index]
                if theta_pol == np.nan:
                    theta_pol = theta  # No calibration data exists, so assume it's perfectly aligned
                mirnov_locations.append((phi, theta, theta_pol))
                mirnov_probes.append(mirnov_name)

                if mirnov_times is None:
                    mirnov_times = mirnov_time
            except Exception as e:
                # Problem with this probe, skip it
                pass

            unchecked_names.remove(mirnov_name)

        return mirnov_times, mirnov_signals, mirnov_locations

    @staticmethod
    @cache_method
    def get_three_mirnov_stffts(params: PhysicsMethodParams):
        """Get the short-time fast Fourier transforms of up to 3 available Mirnov coils.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters for the physics method.

        Returns
        -------
        mirnov_ffts: xarray.DataArray
            The STFFT of the Mirnov coils.
        """

        mirnov_times, mirnov_signals, mirnov_locations = CmodTearingMethods.get_three_mirnov_signals(params)

        # Test to make sure the sampling frequency is consistent
        actual_mirnov_sample_freq = 1 / (mirnov_times[1] - mirnov_times[0])  # Mirnov sampling rate
        if actual_mirnov_sample_freq < 2.4e6 or actual_mirnov_sample_freq > 2.6e6:
            print("Mirnov sample frequency is not 2.5 MHz. Be warned!!!")
            f_mirnov = actual_mirnov_sample_freq
        else:
            f_mirnov = 2.5e6

        f_timebase = 1 / (params.times[1] - params.times[0])  # however fast the timebase is

        SFT = setup_stfft(f_mirnov=f_mirnov, f_target=f_timebase)
        f_indices = np.where(SFT.f < 80e3)
        freqs = SFT.f[f_indices]

        signal_ffts_full = []
        for signal in mirnov_signals:
            signal_fft = SFT.stft(signal)[f_indices]
            signal_ffts_full.append(signal_fft)

        # Interpolate the STFFTs to the timebase
        fft_times = (SFT.delta_t * np.arange(signal_ffts_full[0].shape[1])) + mirnov_times[0]
        signal_ffts_interp = []
        for signal_fft_full in signal_ffts_full:
            signal_fft_interp = interp1(fft_times, signal_fft_full, params.times)
            signal_ffts_interp.append(signal_fft_interp)

        # Make an xarray DataArray for the STFFT's, with dimensions of probe location, frequency, and time
        mirnov_ffts = xr.DataArray(
            np.array(signal_ffts_interp),
            dims=("probe", "frequency", "time"),
            coords={
                "probe": list(range(len(mirnov_locations))),
                "frequency": freqs,
                "time": params.times,
                "phi": ("probe", [loc[0] for loc in mirnov_locations]),
                "theta": ("probe", [loc[1] for loc in mirnov_locations]),
                "theta_pol": ("probe", [loc[2] for loc in mirnov_locations]),
            },
        )
        return mirnov_ffts

    @staticmethod
    def get_complex_cross_spectra(params: PhysicsMethodParams, mirnov_signal_1, mirnov_signal_2):
        """Calculate the complex cross spectra between 2 Mirnov coils.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters for the physics method.
        mirnov_signal_1 : numpy.ndarray
            The signal from the first Mirnov coil.
        mirnov_signal_2 : numpy.ndarray
            The signal from the second Mirnov coil.
        """

        mirnov_signals = CmodTearingMethods.get_three_mirnov_signals(params)

    @staticmethod
    @physics_method(
        tags=["mirnov_spectrogram"],
        tokamak=Tokamak.CMOD,
    )
    def _get_mirnov_spectrogram(params: PhysicsMethodParams):
        """Calculate the spectrogram from a Mirnov coil."""

        mirnov_data = CmodTearingMethods.get_three_mirnov_stffts(params)
        # Pick the first Mirnov coil
        mirnov_data = mirnov_data[0]
        # Get the absolute value of the complex signal
        mirnov_spectrogram = np.abs(mirnov_data) ** 2

        # Convert to a dictionary where the keys are tuples of the signal name and the freuqency,
        # and the values are the spectrogram data for that signal and frequency over time

        spectrogram_dict = {}
        for i, freq in enumerate(mirnov_data.frequency.values):
            spectrogram_dict[("mirnov_spectrogram", freq)] = mirnov_spectrogram[i]
        return spectrogram_dict

