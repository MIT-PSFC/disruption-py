#!/usr/bin/env python3

"""
Module for retrieving and calculating data for HBTEP physics methods.
"""

import numpy as np
from scipy.integrate import cumulative_trapezoid

from disruption_py.core.physics_method.caching import cache_method
from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import (
    butterworth_filter,
    gaussian_low_pass_filter,
    interp1,
)
from disruption_py.machine.hbtep.util import HbtepUtilMethods
from disruption_py.machine.tokamak import Tokamak


class HbtepPhysicsMethods:
    """
    This class provides methods to retrieve and calculate physics-related data
    for HBTEP.
    """

    @staticmethod
    @physics_method(columns=["ip"], tokamak=Tokamak.HBTEP)
    def get_ip(params: PhysicsMethodParams):
        """
        Get the plasma current
        """
        ip, t_ip = params.mds_conn.get_data_with_dims(
            r"\TOP.SENSORS.ROGOWSKIS:IP", tree_name="hbtep2"
        )  # [A], [s]
        ip = interp1(t_ip, ip, params.times, "linear")
        return {"ip": ip}

    @staticmethod
    @physics_method(columns=["i_vfc", "i_ohc"], tokamak=Tokamak.HBTEP)
    def get_cap_bank_parameters(params: PhysicsMethodParams):
        """
        Get the ohmic heating and vertical field coil currents
        """
        i_vfc, t_vfc = params.mds_conn.get_data_with_dims(
            r"\TOP.SENSORS.VF_CURRENT", tree_name="hbtep2"
        )  # [A], [s]
        i_ohc, t_ohc = params.mds_conn.get_data_with_dims(
            r"\TOP.SENSORS.OH_CURRENT", tree_name="hbtep2"
        )  # [A], [s]
        i_vfc = interp1(t_vfc, i_vfc, params.times, "linear")
        i_ohc = interp1(t_ohc, i_ohc, params.times, "linear")
        return {"i_vfc": i_vfc, "i_ohc": i_ohc}

    @staticmethod
    @physics_method(columns=["r", "aminor"], tokamak=Tokamak.HBTEP)
    def get_plasma_radii(params: PhysicsMethodParams):
        """
        Get the major & minor radii
        # TODO: get the updated limitor configuration post ~2023(?) from hbtplot.py
        """
        # Determined by Daisuke during copper plasma calibration
        a = 0.00643005
        b = -1.10423
        c = 48.2567

        # Calculated by Jeff, but still has errors
        vf_pickup = 0.0046315133 * -1e-3
        oh_pickup = 7.0723416e-08

        # get vf and oh data
        i_vfc, t_vfc = params.mds_conn.get_data_with_dims(
            r"\TOP.SENSORS.VF_CURRENT", tree_name="hbtep2"
        )  # [A], [s]
        i_ohc, t_ohc = params.mds_conn.get_data_with_dims(
            r"\TOP.SENSORS.OH_CURRENT", tree_name="hbtep2"
        )  # [A], [s]

        # get plasma current
        ip, t_ip = params.mds_conn.get_data_with_dims(
            r"\TOP.SENSORS.ROGOWSKIS:IP", tree_name="hbtep2"
        )  # [A], [s]
        ip *= 1212.3 * 1e-9  # ip gain

        # get cosine Rogowski data
        cos1_raw, t_cos1_raw = params.mds_conn.get_data_with_dims(
            r"\TOP.SENSORS.ROGOWSKIS:COS_1:RAW", tree_name="hbtep2"
        )  # TODO: units?, TODO: check if we need to specify time range
        # calculate and subtract offset
        (indices,) = np.where((-1e-3 < t_cos1_raw) & (t_cos1_raw < 0))
        cos1_raw_offset = np.mean(cos1_raw[indices])
        cos1_raw -= cos1_raw_offset

        # integrate cos1 raw -- hbteplib has the wrong slicing indices for cos1_raw ([:-1])
        cos1 = cumulative_trapezoid(cos1_raw, t_cos1_raw) + cos1_raw[1:] * 0.004571

        # Interpolate all signals to ip timebase
        i_vfc = interp1(t_vfc, i_vfc, t_ip, "linear")
        i_ohc = interp1(t_ohc, i_ohc, t_ip, "linear")
        cos1 = interp1(t_cos1_raw[1:], cos1, t_ip, "linear")

        # Calculate r major
        pickup = i_vfc * vf_pickup + i_ohc * oh_pickup
        ratio = ip / (cos1 - pickup)
        arg = b**2 - 4 * a * (c - ratio)
        arg[arg < 0] = 0
        r_major = (-b + np.sqrt(arg)) / (2 * a)
        r = r_major / 100  # Convert to meters

        # Calculate aminor
        aminor = np.ones(len(r)) * 0.15
        (outboard_limited_indices,) = np.where(r > 0.92)
        aminor[outboard_limited_indices] = (
            1.07 - r[outboard_limited_indices]
        )  # Outboard limited
        (inboard_limited_indices,) = np.where(r < (0.92 - 0.01704))
        aminor[inboard_limited_indices] = (
            r[inboard_limited_indices] - 0.75296
        )  # inward limited

        # Interpolate to requested timebase
        r = interp1(t_ip, r, params.times, "linear")
        aminor = interp1(t_ip, aminor, params.times, "linear")

        return {"r": r, "aminor": aminor}

    @staticmethod
    @physics_method(columns=["btor"], tokamak=Tokamak.HBTEP)
    def get_btor(params: PhysicsMethodParams):
        """
        Calculate B_tor from the TF probe data
        """
        btor, t_btor = params.mds_conn.get_data_with_dims(
            r"\TOP.SENSORS.TF_PROBE", tree_name="hbtep2"
        )  # [T], [s]
        btor = interp1(t_btor, btor, params.times, "linear")

        r = HbtepPhysicsMethods.get_plasma_radii(params)["r"]
        btor = btor * 1.23 / r
        return {"btor": btor}

    @staticmethod
    @physics_method(columns=["qstar"], tokamak=Tokamak.HBTEP)
    def get_qstar(params: PhysicsMethodParams):
        """
        Calculate the edge safety factor from magnetic measurements
        """
        ip = HbtepPhysicsMethods.get_ip(params)["ip"]
        radii = HbtepPhysicsMethods.get_plasma_radii(params)
        r, aminor = radii["r"], radii["aminor"]
        btor = HbtepPhysicsMethods.get_btor(params)["btor"]

        # Calculate qstar
        qstar = aminor**2 * btor / (2e-7 * ip * r)
        qstar *= 1.15  # 15% correction factor -- Jeff Levesque
        return {"qstar": qstar}

    @staticmethod
    @physics_method(columns=["v_loop"], tokamak=Tokamak.HBTEP)
    def get_v_loop(params: PhysicsMethodParams):
        """
        Get v_loop measurement
        """
        v_loop, t_v_loop = params.mds_conn.get_data_with_dims(
            r"\TOP.SENSORS.LOOP_VOLTAGE", tree_name="hbtep2"
        )  # [V], [s]
        v_loop = interp1(t_v_loop, v_loop, params.times, "linear")
        return {"v_loop": v_loop}

    @staticmethod
    @physics_method(columns=["h_alpha"], tokamak=Tokamak.HBTEP)
    def get_h_alpha(params: PhysicsMethodParams):
        """
        Get D alpha line emission from spectrometer
        """
        h_alpha, t_h_alpha = params.mds_conn.get_data_with_dims(
            r"\TOP.SENSORS.SPECTROMETER", tree_name="hbtep2"
        )  # [arb], [s]
        h_alpha = interp1(t_h_alpha, h_alpha, params.times, "linear")
        return {"h_alpha": h_alpha}

    @staticmethod
    @physics_method(columns=["sxr_midplane"], tokamak=Tokamak.HBTEP)
    def get_sxr_midplane(params: PhysicsMethodParams):
        """
        Get the soft x-ray midplane sensor data
        """
        sxr, t_sxr = params.mds_conn.get_data_with_dims(
            r"\TOP.DEVICES.NORTH_RACK:CPCI:INPUT_74", tree_name="hbtep2"
        )  # [arb], [s]
        # offset subtraction
        (offset_indices,) = np.where((0 < t_sxr) & (t_sxr < 0.5e-3))
        offset = sxr[offset_indices].mean()
        sxr = -(sxr - offset)
        # interpolate to requested timebase
        sxr = interp1(t_sxr, sxr, params.times, "linear")
        return {"sxr_midplane": sxr}

    @staticmethod
    @physics_method(
        columns=[
            "n_equal_1_mode",
            "n_equal_1_normalized",
            "n_equal_1_phase",
            "n_equal_1_frequency",
            "n_equal_2_mode",
            "n_equal_2_normalized",
            "n_equal_2_phase",
            "n_equal_2_frequency",
        ],
        tokamak=Tokamak.HBTEP,
    )
    def get_n_mode_data(params: PhysicsMethodParams):
        """
        Get the n=1 and 2 mode amplitdes, phases, and frequencies from the toroidal array (TA)
        # TODO: implement this function
        """
        # Get TA data
        ta_data = HbtepPhysicsMethods._get_ta_data(params)
        data = ta_data["ta_pol_data_filt"]
        time = ta_data["ta_pol_time"]
        phi = ta_data["ta_pol_phi"]
        n, m = len(data), len(
            data[0]
        )  # TODO: decide if we need to convert data to np.ndarray

        # Construct A matrix and calculate its inversion
        A = np.zeros((n, 5))
        A[:, 0] = np.ones(n)
        A[:, 1] = np.sin(phi)
        A[:, 2] = np.cos(phi)
        A[:, 3] = np.sin(2 * phi)
        A[:, 4] = np.cos(2 * phi)
        Ainv = np.linalg.pinv(A)

        # Use least squares to solve for coefficients
        x = np.matmul(Ainv, data)
        n1_amp = np.sqrt(x[1, :] ** 2 + x[2, :] ** 2)  # [T]
        n2_amp = np.sqrt(x[1, :] ** 2 + x[2, :] ** 2)
        n1_phase = -np.arctan2(x[1, :], x[2, :])  # *-1 to corrected the slope of phase
        n2_phase = -np.arctan2(x[3, :], x[4, :])

        # Calculate mode frequency. Note that gaussian_low_pass_filter is not a causal filter.
        n1_phase_filt = gaussian_low_pass_filter(
            HbtepUtilMethods.unwrap_phase(n1_phase), time, 1e-5
        )
        n2_phase_filt = gaussian_low_pass_filter(
            HbtepUtilMethods.unwrap_phase(n2_phase), time, 1e-5
        )
        n1_freq = np.gradient(n1_phase_filt) / np.gradient(time) / (2 * np.pi)
        n2_freq = np.gradient(n2_phase_filt) / np.gradient(time) / (2 * np.pi)

        # Interpolate output data to the requested timebase
        n1_amp = interp1(time, n1_amp, params.times, "linear")
        n2_amp = interp1(time, n2_amp, params.times, "linear")
        n1_phase = interp1(time, n1_phase, params.times, "linear")
        n2_phase = interp1(time, n2_phase, params.times, "linear")
        n1_freq = interp1(time, n1_freq, params.times, "linear")
        n2_freq = interp1(time, n2_freq, params.times, "linear")

        # Calculate normalized mode amplitudes
        btor = HbtepPhysicsMethods.get_btor(params)["btor"]
        n1_amp_norm = n1_amp / btor
        n2_amp_norm = n2_amp / btor

        return {
            "n_equal_1_mode": n1_amp,
            "n_equal_1_normalized": n1_amp_norm,
            "n_equal_1_phase": n1_phase,
            "n_equal_1_frequency": n1_freq,
            "n_equal_2_mode": n2_amp,
            "n_equal_2_normalized": n2_amp_norm,
            "n_equal_2_phase": n2_phase,
            "n_equal_2_frequency": n2_freq,
        }

    @staticmethod
    @cache_method
    def _get_ta_data(params: PhysicsMethodParams):
        r"""
        Get poloidal and radial field measurements from the toroidal array (TA) sensors
        TA sensors are stored in \HBTEP2::TOP.SENSORS.MAGNETIC:[SENSOR_NAME]
        Return TA pol and rad data in their original timebases
        """
        output = {}
        # Sensor names
        output["ta_pol_names"] = [
            "TA01_S1P",
            "TA01_S2P",
            "TA01_S3P",
            "TA02_S1P",
            "TA02_S2P",
            "TA02_S3P",
            "TA03_S1P",
            "TA03_S2P",
            "TA03_S3P",
            "TA04_S1P",
            "TA04_S2P",
            "TA04_S3P",
            "TA05_S1P",
            "TA05_S2P",
            "TA05_S3P",
            "TA06_S1P",
            "TA06_S2P",
            "TA06_S3P",
            "TA07_S1P",
            "TA07_S2P",
            "TA07_S3P",
            "TA08_S1P",
            "TA08_S2P",
            "TA08_S3P",
            "TA09_S1P",
            "TA09_S2P",
            "TA09_S3P",
            "TA10_S1P",
            "TA10_S2P",
            "TA10_S3P",
        ]
        output["ta_rad_names"] = [
            "TA01_S2R",
            "TA02_S2R",
            "TA03_S2R",
            "TA04_S2R",
            "TA05_S2R",
            "TA06_S2R",
            "TA07_S2R",
            "TA08_S2R",
            "TA09_S2R",
            "TA10_S2R",
        ]
        # Toroidal and poloidal positions of the poloidal sensors
        output["ta_pol_phi"] = (
            np.pi
            / 180
            * np.array(
                [
                    241.5,
                    250.5,
                    259.5,
                    277.5,
                    286.5,
                    295.5,
                    313.5,
                    322.5,
                    331.5,
                    349.5,
                    358.5,
                    7.5,
                    25.5,
                    34.5,
                    43.5,
                    61.5,
                    70.5,
                    79.5,
                    97.5,
                    106.5,
                    115.5,
                    133.5,
                    142.5,
                    151.5,
                    169.5,
                    178.5,
                    187.5,
                    205.5,
                    214.5,
                    223.5,
                ]
            )
        )
        output["ta_pol_theta"] = (
            np.ones(len(output["ta_pol_phi"])) * (189 - 360) * np.pi / 180
        )
        # Toroidal and poloidal positions of the radial sensors -- not implemented in hbteplib

        # Get TA time
        output["ta_pol_time"] = params.mds_conn.get_dims(
            r"\TOP.SENSORS.MAGNETIC:TA01_S1P", tree_name="hbtep2"
        )[
            0
        ]  # [s]
        output["ta_rad_time"] = params.mds_conn.get_dims(
            r"\TOP.SENSORS.MAGNETIC:TA01_S2R", tree_name="hbtep2"
        )[
            0
        ]  # [s]
        fs_ta_pol = (
            1 / (output["ta_pol_time"][1] - output["ta_pol_time"][0])
        ).round()  # should be 1/(2e-6)=500kHz
        fs_ta_rad = (1 / (output["ta_rad_time"][1] - output["ta_rad_time"][0])).round()

        # Get TA data
        ta_pol_data_raw, ta_pol_data_filt = [], []
        for sensor_name in output["ta_pol_names"]:
            sensor_data_raw = params.mds_conn.get_data(
                r"\TOP.SENSORS.MAGNETIC:" + sensor_name, tree_name="hbtep2"
            )  # TODO: unit
            sensor_data_filt = butterworth_filter(
                sensor_data_raw, fs=fs_ta_pol, cutoff=2e3, order=2, btype="high"
            )
            ta_pol_data_raw.append(sensor_data_raw)
            ta_pol_data_filt.append(sensor_data_filt)
        ta_rad_data_raw, ta_rad_data_filt = [], []
        for sensor_name in output["ta_rad_names"]:
            sensor_data_raw = params.mds_conn.get_data(
                r"\TOP.SENSORS.MAGNETIC:" + sensor_name, tree_name="hbtep2"
            )  # TODO: unit
            sensor_data_filt = butterworth_filter(
                sensor_data_raw, fs=fs_ta_rad, cutoff=2e3, order=2, btype="high"
            )
            ta_rad_data_raw.append(sensor_data_raw)
            ta_rad_data_filt.append(sensor_data_filt)
        output["ta_pol_data_raw"] = ta_pol_data_raw
        output["ta_pol_data_filt"] = ta_pol_data_filt
        output["ta_rad_data_raw"] = ta_rad_data_raw
        output["ta_rad_data_filt"] = ta_rad_data_filt

        return output
