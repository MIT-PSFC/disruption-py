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

CHECKED_N_MODES = list(range(2, 6))  # Disregard n=1, go up to n=9. Anything higher is getting in the realm of "maybe that's just a really slow TAE"
#CHECKED_N_MODES = list(range(11, 12))
CHECKED_M_MODES = list(range(1, 11))  # TODO(ZanderKeith) This is bad. We have zero clue what the m number is.
PHASE_TOLERANCE = np.deg2rad(10)       # The tolerance for the phase to be considered a match

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

def get_expected_cross_phase(mode: tuple[int, int], theta_1: float, phi_1: float, theta_2: float, phi_2: float, cylindrical=True):
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

    if cylindrical:
        poloidal_phase = mode[0] * (theta_2 - theta_1)
        toroidal_phase = mode[1] * (phi_2 - phi_1)
        expected_cross_phase = (poloidal_phase + toroidal_phase) % (2 * np.pi)
    else:
        # The phase of the mode must take into account toroidal curvature
        # The full expression is m(theta - lambda sin theta) + n phi https://www.ias.ac.in/article/fulltext/pram/055/05-06/0727-0732
        # Where lambda = (a/R) * (poloidal beta + internal inducatance/2 + 1)
        # TODO(ZanderKeith): this implementation is a placeholder, need to get the actual values for lambda
        R = 0.68
        a = 0.22
        epsilon = a / R

        probe_1_phase = mode[0] * (theta_1 - epsilon * np.sin(theta_1)) + mode[1] * phi_1
        probe_2_phase = mode[0] * (theta_2 - epsilon * np.sin(theta_2)) + mode[1] * phi_2

        expected_cross_phase = (probe_2_phase - probe_1_phase) % (2 * np.pi)

    return float(expected_cross_phase)

def cross_phase_mask(original_signal, probe_cross_phase, mask_cross_phases, tolerance=PHASE_TOLERANCE, inverted=False):
    """Mask the signal based on the cross phases between multiple probes
    Things to look into: Change tolerance based on the mode number (higher q gets more lenience)
    somehow correct for the curvature of the torus
    only filter by n number for certain things
    only filter by odd vs even m and odd vs even n

    Parameters:
    ----------
    original_signal : numpy.ndarray
        The original signal to mask
    probe_cross_phase : np.ndarray
        The cross phase between two mirnov probes
    mask_cross_phases : list[float]
        The cross phases that are allowed to be kept
    tolerance : float, optional (default=PHASE_TOLERANCE) [rad]
        The tolerance for the cross phase to be considered a match
    inverted : bool, optional (default=False)
        If True, will invert the mask and remove the values in mask_cross_phases

    Returns:
    ----------
    masked_signal : numpy.ndarray
        The original signal with the masked values set to np.nan
    """

    if not inverted:
        valid_signal = np.zeros_like(original_signal)
    else:
        valid_signal = np.ones_like(original_signal)

    for cross_phase in mask_cross_phases:
        # We know the probe's cross phase should be sufficiently different from 0 or 2pi
        # But if the expected cross phase is near 0 or 2pi, we risk a 2 pi jump messing things up
        # Gotta bring everything into the center of the circle
        if cross_phase < np.pi / 8:
            tested_cross_phase = cross_phase + np.pi
            tested_probe_cross_phase = (probe_cross_phase + np.pi) % (2 * np.pi)
        elif cross_phase > 15 * np.pi / 8:
            tested_cross_phase = cross_phase - np.pi
            tested_probe_cross_phase = (probe_cross_phase - np.pi) % (2 * np.pi)
        else:
            tested_cross_phase = cross_phase
            tested_probe_cross_phase = probe_cross_phase
        cross_phase_difference = np.abs(tested_cross_phase - tested_probe_cross_phase)

        if not inverted:
            valid_signal = np.where(cross_phase_difference < tolerance, 1, valid_signal)
        else:
            valid_signal = np.where(cross_phase_difference < tolerance, 0, valid_signal)

    masked_signal = np.where(valid_signal == 1, original_signal, np.nan)

    return masked_signal

def aliasing_check(phi_1: float, phi_2: float, n_modes_of_interest: list[int], n_avoided_modes: list[int], tolerance=PHASE_TOLERANCE, self_aliasing=False):
    """Check if the toroidal separation between two probes causes aliasing between the modes of interest and the avoided modes.

    Parameters:
    ----------
    phi_1 : float
        The toroidal angle of the first mirnov probe [rad]
    phi_2 : float
        The toroidal angle of the second mirnov probe [rad]
    n_modes_of_interest : list[int]
        The modes that are of interest
    n_avoided_modes : list[int]
        The modes to avoid aliasing with
    tolerance : float, optional (default=PHASE_TOLERANCE) [rad]
        The tolerance for the phase to be considered a match
    self_aliasing : bool, optional (default=False)
        If True, will also check for aliasing between the modes of interest

    Returns:
    ----------
    aliasing : bool
        True if the toroidal separation causes aliasing between the modes of interest and the avoided modes, False otherwise
    """

    interest_n_cross_phases = [get_expected_cross_phase((0, n), 0, phi_1, 0, phi_2) for n in n_modes_of_interest]
    avoided_n_cross_phases = [get_expected_cross_phase((0, n), 0, phi_1, 0, phi_2) for n in n_avoided_modes]

    for interest_cross_phase in interest_n_cross_phases:
        for avoided_cross_phase in avoided_n_cross_phases:
            if np.abs(interest_cross_phase - avoided_cross_phase) < tolerance:
                return True

    if self_aliasing:
        for interest_cross_phase_1 in interest_n_cross_phases:
            for interest_cross_phase_2 in interest_n_cross_phases:
                if np.abs(interest_cross_phase_1 - interest_cross_phase_2) < tolerance:
                    return True

    return False

class CmodTearingMethods:
    @staticmethod
    @cache_method
    def get_mirnov_names_and_locations(params: PhysicsMethodParams):
        """Get the names and locations of the Mirnov coils.

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

        # How do you do this in disruption-py?
        # I want the equivalent of .getNode, or .getNodeWild, and things like that
        # something along these lines: mirnov_coil_names = params.mds_connection.get_coil_names()
        # TODO(ZanderKeith): expand for every Mirnov coil. For now, just using the ABK, GHK, and K coils
        mirnov_names_ab = [f"BP0{p}_ABK" for p in range(1, 10)] + [f"BP{p}_ABK" for p in range(10, 19)]
        mirnov_names_gh = [f"BP0{p}_GHK" for p in range(1, 10)] + [f"BP{p}_GHK" for p in range(10, 19)]
        mirnov_names_k = [f"BP0{p}_K" for p in range(1, 7)]
        all_mirnov_names = mirnov_names_ab + mirnov_names_gh + mirnov_names_k

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

        return all_mirnov_names, phi_all, theta_all, theta_pol_all

    @staticmethod
    @cache_method
    def get_mirnov_stfft(params: PhysicsMethodParams, mirnov_name: str):
        """Get and cache the interpolated stfft of a Mirnov coil.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters for the physics method.
        mirnov_name : str
            The name of the Mirnov coil to get the STFFT of.

        Returns
        -------
        mirnov_fft_interp : numpy.ndarray
            The interpolated STFFT of the Mirnov coil
        freqs : numpy.ndarray
            The frequencies of the STFFT

        """
        path = r"\magnetics::top.active_mhd.signals"
        try:
            mirnov_signal, mirnov_times = params.mds_conn.get_data_with_dims(
                        path=f"{path}.{mirnov_name}",
                        tree_name="magnetics",
            )
            # Get the sampling frequency of the Mirnov signal
            f_mirnov = 1 / np.mean(np.diff(mirnov_times))
            print(f"Using Mirnov frequency of {f_mirnov} Hz on {mirnov_name}")
            # Print a warning if the sample rate is very slow
            if f_mirnov < 2.4e6:
                params.logger.warning(f"[Shot {params.shot_id}] You're too slow! Got Mirnov frequency of {f_mirnov} Hz on {mirnov_name}, expected 2.5 MHz or 5 MHz")
            # Interpolate the signal to the actual correct Mirnov frequency
            #f_mirnov = 2.5e6
            #mirnov_signal = interp1(mirnov_times, mirnov_signal, np.arange(mirnov_times[0], mirnov_times[-1], 1/f_mirnov))

            f_timebase = 1 / (params.times[1] - params.times[0])  # however fast the timebase is

            SFT = setup_stfft(f_mirnov=f_mirnov, f_target=f_timebase)
            f_indices = np.where(SFT.f < 80e3)
            freqs = SFT.f[f_indices]

            mirnov_fft_full = SFT.stft(mirnov_signal)[f_indices]

            # Interpolate the STFFTs to the timebase
            fft_times = (SFT.delta_t * np.arange(mirnov_fft_full.shape[1])) + mirnov_times[0]
            mirnov_fft_interp = interp1(fft_times, mirnov_fft_full, params.times)

            # Check if the average difference between frequencies is NOT close to 250
            if not np.isclose(np.mean(np.diff(freqs)), 250, atol=1):
                params.logger.warning(f"[Shot {params.shot_id}] The Mirnov frequency resolution is not 250 Hz")
            # Replace the frequencies with nice integers (for the sake of consistency)
            freqs = np.arange(0, 80e3, 250)

            if (mirnov_fft_interp.shape[0] != 320):
                params.logger.warning(f"[Shot {params.shot_id}] The Mirnov signal is not 320 samples long. This may cause issues with the n mode filtering.")
            return mirnov_fft_interp, freqs
        except Exception as e:
            return None, None

    @staticmethod
    @cache_method
    def get_two_mirnov_ffts(params: PhysicsMethodParams):  # noqa: PLR0915
        """Get measurements from two available Mirnov coils to be used for n mode filtering.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters for the physics method.

        Returns
        -------
        mirnov_ffts: xarray.DataArray
            The STFFT of the Mirnov coils.
        """

        all_mirnov_names, phi_all, theta_all, theta_pol_all = CmodTearingMethods.get_mirnov_names_and_locations(params)

        # Get all pairs of probes where their toroidal separation doesn't cause aliasing with the n=1 mode
        # And sort them by the poloidal distance
        candidate_pairs = [(i, j) for i in range(len(all_mirnov_names)) for j in range(i + 1, len(all_mirnov_names)) if (aliasing_check(phi_all[i], phi_all[j], CHECKED_N_MODES, [1]) is False)]
        candidate_pairs = sorted(candidate_pairs, key=lambda x: np.abs(theta_all[x[0]] - theta_all[x[1]]))

        # Get the signals from the first pair that works
        # Eliminate pairs where one of the probes doesn't work
        while len(candidate_pairs) > 0:
            candidate_pair = candidate_pairs.pop(0)
            mirnov_1_index = candidate_pair[0]
            mirnov_2_index = candidate_pair[1]
            mirnov_1_name = all_mirnov_names[mirnov_1_index]
            mirnov_2_name = all_mirnov_names[mirnov_2_index]

            mirnov_1_fft, freqs = CmodTearingMethods.get_mirnov_stfft(params=params, mirnov_name=mirnov_1_name)
            if mirnov_1_fft is None:
                candidate_pairs = [pair for pair in candidate_pairs if mirnov_1_index not in pair]
                continue
            mirnov_2_fft, _ = CmodTearingMethods.get_mirnov_stfft(params=params, mirnov_name=mirnov_2_name)
            if mirnov_2_fft is None:
                candidate_pairs = [pair for pair in candidate_pairs if mirnov_2_index not in pair]
                continue

            mirnov_ffts = np.array([mirnov_1_fft, mirnov_2_fft])
            if np.abs(theta_pol_all[mirnov_1_index] - theta_pol_all[mirnov_2_index]) > PHASE_TOLERANCE:
                params.logger.warning(f"[Shot {params.shot_id}] The poloidal angles of the probes are not aligned. This may cause issues with the n mode filtering.")

            mirnov_1_location = (phi_all[mirnov_1_index], theta_all[mirnov_1_index], theta_pol_all[mirnov_1_index])
            mirnov_2_location = (phi_all[mirnov_2_index], theta_all[mirnov_2_index], theta_pol_all[mirnov_2_index])
            mirnov_locations = [mirnov_1_location, mirnov_2_location]

            # Make an xarray DataArray for the STFFT's, with dimensions of probe location, frequency, and time
            mirnov_ffts = xr.DataArray(
                mirnov_ffts,
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

        # If we get here, we couldn't find a pair of probes that worked
        return None

    @staticmethod
    @cache_method
    def get_mirnov_spectrogram(params: PhysicsMethodParams):
        all_mirnov_names, _, _, _ = CmodTearingMethods.get_mirnov_names_and_locations(params)
        for mirnov_name in all_mirnov_names:
            mirnov_data, freqs = CmodTearingMethods.get_mirnov_stfft(params=params, mirnov_name=mirnov_name)
            if mirnov_data is not None:
                break
        # Get the absolute value of the complex signal
        mirnov_spectrogram = np.abs(mirnov_data) ** 2

        return mirnov_spectrogram, freqs

    @staticmethod
    @physics_method(
        tags=["mirnov_spectrogram"],
        tokamak=Tokamak.CMOD,
    )
    def _get_mirnov_spectrogram(params: PhysicsMethodParams):
        """Calculate the spectrogram from a Mirnov coil."""

        mirnov_spectrogram, freqs = CmodTearingMethods.get_mirnov_spectrogram(params)

        # Convert to a dictionary where the keys are tuples of the signal name and the freuqency,
        # and the values are the spectrogram data for that signal and frequency over time
        spectrogram_dict = {}
        for i, freq in enumerate(freqs):
            spectrogram_dict[("mirnov_spectrogram", freq)] = mirnov_spectrogram[i]
        return spectrogram_dict

    @staticmethod
    @physics_method(
        tags=["n_mode_spectrogram"],
        tokamak=Tokamak.CMOD,
    )
    def _get_n_mode_spectrograms(params: PhysicsMethodParams):
        """Calculate the spectrogram from a Mirnov coil that is masked to only show specific n modes."""

        mirnov_data = CmodTearingMethods.get_two_mirnov_ffts(params)
        original_spectrogram, freqs = CmodTearingMethods.get_mirnov_spectrogram(params)

        # Find the cross phase between the probes
        cross_phase = np.angle(mirnov_data[0] * np.conj(mirnov_data[1])) % (2 * np.pi)

        # For each mode, mask the original spectrogram
        masked_spectrograms = {}
        for mode in CHECKED_N_MODES:
            expected_cross_phase = get_expected_cross_phase((0, mode), mirnov_data.theta[0], mirnov_data.phi[0], mirnov_data.theta[1], mirnov_data.phi[1])
            masked_spectrogram = cross_phase_mask(original_spectrogram, [cross_phase], [expected_cross_phase])
            masked_spectrograms[mode] = masked_spectrogram

        # Make a dictionary where the keys are ("mode_spectrogram", mode, frequency) and the values are the masked spectrogram data
        mode_spectrogram_dict = {}
        for i, freq in enumerate(freqs):
            for mode in CHECKED_N_MODES:
                mode_spectrogram_dict[("n_mode_spectrogram", mode, freq)] = masked_spectrograms[mode][i]

        return mode_spectrogram_dict

    @staticmethod
    @physics_method(
        tags=["n_greater_1_spectrogram"],
        tokamak=Tokamak.CMOD,
    )
    def _get_n_mode_spectrograms(params: PhysicsMethodParams):
        """Calculate the spectrogram from a Mirnov coil that is masked to show modes CHECKED_N_MODES."""

        mirnov_data = CmodTearingMethods.get_two_mirnov_ffts(params)
        original_spectrogram, freqs = CmodTearingMethods.get_mirnov_spectrogram(params)

        # Find the cross phase between the probes
        cross_phase = np.angle(mirnov_data[0] * np.conj(mirnov_data[1])) % (2 * np.pi)

        # For each mode, get the expected cross phase
        expected_cross_phases = [get_expected_cross_phase((0, mode), 0, mirnov_data.phi[0], 0, mirnov_data.phi[1]) for mode in CHECKED_N_MODES]

        masked_spectrogram = cross_phase_mask(original_spectrogram, cross_phase, expected_cross_phases)
        # Remove anything remotely close to n=1
        n_1_cross_phase = get_expected_cross_phase((0, 1), 0, mirnov_data.phi[0], 0, mirnov_data.phi[1])
        masked_spectrogram = cross_phase_mask(masked_spectrogram, cross_phase, [n_1_cross_phase], PHASE_TOLERANCE, inverted=True)

        # Make a dictionary where the keys are ("n_greater_1_spectrogram", frequency) and the values are the masked spectrogram data
        spectrogram_dict = {}
        for i, freq in enumerate(freqs):
            spectrogram_dict[("n_greater_1_spectrogram", freq)] = masked_spectrogram[i]

        return spectrogram_dict

    @staticmethod
    @physics_method(
        tags=["mn_mode_spectrogram"],
        tokamak=Tokamak.CMOD,
    )
    def _get_mn_mode_spectrograms(params: PhysicsMethodParams):
        """Calculate the spectrogram from a Mirnov coil that is masked to only show specific m/n modes."""
        # Everything up to 5/4
        searched_modes = [(m, n) for n in range(1, 10) for m in range(n, 11)]

        mirnov_data = CmodTearingMethods.get_three_mirnov_stffts(params)
        original_spectrogram = CmodTearingMethods.get_mirnov_spectrogram(params)

        # Find the cross phase between the first probe and the other two
        cross_phase_1_2 = np.angle(mirnov_data[0] * np.conj(mirnov_data[1])) % (2 * np.pi)
        cross_phase_1_3 = np.angle(mirnov_data[0] * np.conj(mirnov_data[2])) % (2 * np.pi)

        # For each mode, mask the original spectrogram
        masked_spectrograms = {}
        for mode in searched_modes:
            expected_cross_phase_1_2 = get_expected_cross_phase(mode, mirnov_data.theta[0], mirnov_data.phi[0], mirnov_data.theta[1], mirnov_data.phi[1], cylindrical=False)
            expected_cross_phase_1_3 = get_expected_cross_phase(mode, mirnov_data.theta[0], mirnov_data.phi[0], mirnov_data.theta[2], mirnov_data.phi[2], cylindrical=False)

            masked_spectrogram = cross_phase_mask(original_spectrogram.values, [cross_phase_1_2, cross_phase_1_3], [expected_cross_phase_1_2, expected_cross_phase_1_3])
            masked_spectrograms[mode] = masked_spectrogram

        # Make a dictionary where the keys are ("mode_spectrogram", mode, frequency) and the values are the masked spectrogram data
        mode_spectrogram_dict = {}
        for i, freq in enumerate(original_spectrogram.frequency.values):
            for mode in searched_modes:
                mode_spectrogram_dict[("mode_spectrogram", mode, freq)] = masked_spectrograms[mode][i]

        return mode_spectrogram_dict




