#!/usr/bin/env python3

import traceback
import warnings
from importlib import resources

import numpy as np
import pandas as pd
from MDSplus import mdsExceptions
from scipy.signal import ShortTimeFFT
from scipy.signal.windows import kaiser

import disruption_py.data
from disruption_py.core.physics_method.caching import cache_method
from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import gaussian_fit, interp1, smooth
from disruption_py.machine.tokamak import Tokamak


class CmodTearingMethods:

    @staticmethod
    @cache_method
    def get_three_mirnov_signals(params: PhysicsMethodParams):
        """Get measurements from up to n available Mirnov coils.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters for the physics method.
        n : int
            The number of Mirnov coils to retrieve.

        Returns
        -------
        mirnov_times : numpy.ndarray
            The times of the Mirnov signals.
        mirnov_signals: list[numpy.ndarray]
            The signals from the Mirnov coils.
        probe_locations: list[tuple]
            The (phi, theta, theta_offset) locations of the probes.
        """

        # Get all names of the Mirnov coils TODO: how do you do this in disruption-py?
        # I want the equivalent of .getNode, or .getNodeWild, and things like that
        # mirnov_coil_names = params.mds_connection.get_coil_names()
        # TODO: expand for every Mirnov coil. For now, just using the ABK and GHK coils
        mirnov_names_ab = [f"BP0{p}_ABK" for p in range(1, 10)] + [f"BP{p}_ABK" for p in range(10, 19)]
        mirnov_names_gh = [f"BP0{p}_GHK" for p in range(1, 10)] + [f"BP{p}_GHK" for p in range(10, 19)]
        mirnov_names_k = [f"BP0{p}_K" for p in range(1, 7)]
        path = r"\magnetics::top.active_mhd.signals"

        phi_ab, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.phi_AB", tree_name="magnetics")
        phi_gh, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.phi_GH", tree_name="magnetics")
        phi_k, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.phi_K", tree_name="magnetics")
        phi_all = np.concatenate((phi_ab, phi_gh, phi_k))

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
        # theta_pol_k, _ = params.mds_conn.get_data_with_dims(r"\magnetics::top.rf_lim_coils.theta_pol_K", tree_name="magnetics")
        theta_pol_k = np.empty_like(Z_k) # Calibration data doesn't exist, so we'll just say NaN and deal with it later
        theta_pol_all = np.concatenate((theta_pol_ab, theta_pol_gh, theta_pol_k))

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
                        all_distances[i, j] = np.sqrt(
                            (phi - phi_all[probe_index]) ** 2 + (theta - theta_all[probe_index]) ** 2
                        )
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

                mirnov_locations.append((phi, theta, theta_pol))
                mirnov_probes.append(mirnov_name)

                if mirnov_times is None:
                    mirnov_times = mirnov_time
            except Exception as e:
                pass

            unchecked_names.remove(mirnov_name)

        return mirnov_times, mirnov_signals, mirnov_locations

    @staticmethod
    def get_complex_cross_spectra(params: PhysicsMethodParams, mirnov_signal_1, mirnov_signal_2, mirnov_signal_3):
        """Calculate the complex cross spectra between 3 Mirnov coils.
        If any of the signals are not provided, the method will return NaN's in calculations which required them.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters for the physics method.
        mirnov_times : numpy.ndarray
            The times of the Mirnov signals.
        mirnov_signal_1 : numpy.ndarray
            The signal from the first Mirnov coil.
        mirnov_signal_2 : numpy.ndarray
            The signal from the second Mirnov coil.
        mirnov_signal_3 : numpy.ndarray
            The signal from the third Mirnov coil.
        """

        mirnov_signals = CmodTearingMethods.get_three_mirnov_signals(params)

    @staticmethod
    @physics_method(
        tags=["mirnov_spectrogram"],
        tokamak=Tokamak.CMOD,
    )
    def _get_mirnov_spectrogram(params: PhysicsMethodParams):
        """Calculate the spectrogram from a Mirnov coil."""

        CmodTearingMethods.get_complex_cross_spectra(params, 1, 2, 3)
        return 5


