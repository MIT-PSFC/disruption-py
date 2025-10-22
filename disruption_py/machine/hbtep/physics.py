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
    @cache_method
    @physics_method(columns=["ip"], tokamak=Tokamak.HBTEP)
    def get_ip(params: PhysicsMethodParams):
        """
        Get the plasma current
        """
        ip, t_ip = params.mds_conn.get_data_with_dims(
            r"\top.sensors.rogowskis:ip", tree_name="hbtep2"
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
            r"\top.sensors.vf_current", tree_name="hbtep2"
        )  # [A], [s]
        i_ohc, t_ohc = params.mds_conn.get_data_with_dims(
            r"\top.sensors.oh_current", tree_name="hbtep2"
        )  # [A], [s]
        i_vfc = interp1(t_vfc, i_vfc, params.times, "linear")
        i_ohc = interp1(t_ohc, i_ohc, params.times, "linear")
        return {"i_vfc": i_vfc, "i_ohc": i_ohc}

    @staticmethod
    @cache_method
    @physics_method(columns=["r", "aminor"], tokamak=Tokamak.HBTEP)
    def get_plasma_radii(params: PhysicsMethodParams):
        """
        Get the major & minor radii
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
            r"\top.sensors.vf_current", tree_name="hbtep2"
        )  # [A], [s]
        i_ohc, t_ohc = params.mds_conn.get_data_with_dims(
            r"\top.sensors.oh_current", tree_name="hbtep2"
        )  # [A], [s]

        # get plasma current
        ip, t_ip = params.mds_conn.get_data_with_dims(
            r"\top.sensors.rogowskis:ip", tree_name="hbtep2"
        )  # [A], [s]
        ip *= 1212.3 * 1e-9  # ip gain

        # get cosine Rogowski data
        cos1_raw, t_cos1_raw = params.mds_conn.get_data_with_dims(
            r"\top.sensors.rogowskis:cos_1:raw", tree_name="hbtep2"
        )  # [A/s], [s]
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
        # Default: top-down limited, aminor=0.15
        aminor = np.ones(len(r)) * 0.15
        # Get outboard-limited time points
        (outboard_limited_indices,) = np.where(r > 0.92)
        aminor[outboard_limited_indices] = 1.07 - r[outboard_limited_indices]
        # Get inboard-limited time points
        if params.shot_id >= 117122:
            # Use new limiting surface post 2023 HFS SOL tiles installation
            (inboard_limited_indices,) = np.where(r < 0.92)
            aminor[inboard_limited_indices] = np.sqrt(
                (r[inboard_limited_indices] - 0.77244) ** 2 + 0.03478**2
            )
            aminor[aminor > 0.15] = 0.15
        else:
            (inboard_limited_indices,) = np.where(r < (0.92 - 0.01704))
            aminor[inboard_limited_indices] = r[inboard_limited_indices] - 0.75296

        # Interpolate to requested timebase
        r = interp1(t_ip, r, params.times, "linear")
        aminor = interp1(t_ip, aminor, params.times, "linear")

        return {"r": r, "aminor": aminor}

    @staticmethod
    @cache_method
    @physics_method(columns=["btor"], tokamak=Tokamak.HBTEP)
    def get_btor(params: PhysicsMethodParams):
        """
        Calculate B_tor from the TF probe data
        """
        btor, t_btor = params.mds_conn.get_data_with_dims(
            r"\top.sensors.tf_probe", tree_name="hbtep2"
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
            r"\top.sensors.loop_voltage", tree_name="hbtep2"
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
            r"\top.sensors.spectrometer", tree_name="hbtep2"
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
            r"\top.devices.north_rack:cpci:input_74", tree_name="hbtep2"
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
        """
        # Get TA data
        ta_data = HbtepPhysicsMethods._get_ta_data(params)

        # Get the shortest time array and interpolate every signal to that time base
        time = ta_data["ta_pol_time"][0]
        for t_signal in ta_data["ta_pol_time"][1:]:
            if len(t_signal) < len(time):
                time = t_signal
        data = []
        for signal, t_signal in zip(
            ta_data["ta_pol_data_filt"], ta_data["ta_pol_time"]
        ):
            data.append(interp1(t_signal, signal, time))
        n = len(data)
        phi = ta_data["ta_pol_phi"]

        # Construct A matrix and calculate its inversion
        a_matrix = np.zeros((n, 5))
        a_matrix[:, 0] = np.ones(n)
        a_matrix[:, 1] = np.sin(phi)
        a_matrix[:, 2] = np.cos(phi)
        a_matrix[:, 3] = np.sin(2 * phi)
        a_matrix[:, 4] = np.cos(2 * phi)
        a_matrix_inv = np.linalg.pinv(a_matrix)

        # Use least squares to solve for coefficients
        x = np.matmul(a_matrix_inv, data)
        n1_amp = np.sqrt(x[1, :] ** 2 + x[2, :] ** 2)  # [T]
        n2_amp = np.sqrt(x[3, :] ** 2 + x[4, :] ** 2)
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
    @physics_method(
        columns=[
            "m_equal_2_mode",
            "m_equal_2_normalized",
            "m_equal_2_phase",
            "m_equal_2_frequency",
            "m_equal_3_mode",
            "m_equal_3_normalized",
            "m_equal_3_phase",
            "m_equal_3_frequency",
            "m_equal_4_mode",
            "m_equal_4_normalized",
            "m_equal_4_phase",
            "m_equal_4_frequency",
        ],
        tokamak=Tokamak.HBTEP,
    )
    def get_m_mode_data(params: PhysicsMethodParams):
        """
        Get the m=2,3,4 mode amplitdes, phases, and frequencies from the poloidal array (PA1)
        """
        # Get PA1 data
        pa_data = HbtepPhysicsMethods._get_pa_data(params)

        # Get the shortest time array and interpolate every signal to that time base
        time = pa_data["pa1_time"][0]
        for t_signal in pa_data["pa1_time"][1:]:
            if len(t_signal) < len(time):
                time = t_signal
        data = []
        for signal, t_signal in zip(pa_data["pa1_data_filt"], pa_data["pa1_time"]):
            data.append(interp1(t_signal, signal, time))
        n = len(data)
        theta = pa_data["pa1_theta"]
        phi = pa_data["pa1_phi"]

        # Construct A matrix and calculate its inversion
        a_matrix = np.zeros((n, 11))
        a_matrix[:, 0] = np.ones(n)
        for i in range(1, 6):
            a_matrix[:, i * 2 - 1] = np.sin(i * theta - phi)
            a_matrix[:, i * 2] = np.cos(i * theta - phi)
        a_matrix_inv = np.linalg.pinv(a_matrix)

        # Use least squares to solve for coefficients
        x = np.matmul(a_matrix_inv, data)
        m2_amp = np.sqrt(x[3, :] ** 2 + x[4, :] ** 2)  # [T]
        m3_amp = np.sqrt(x[5, :] ** 2 + x[6, :] ** 2)
        m4_amp = np.sqrt(x[7, :] ** 2 + x[8, :] ** 2)
        m2_phase = np.arctan2(x[3, :], x[4, :])
        m3_phase = np.arctan2(x[5, :], x[6, :])
        m4_phase = np.arctan2(x[7, :], x[8, :])

        # Calculate mode frequency. Note that gaussian_low_pass_filter is not a causal filter.
        m2_phase_filt = gaussian_low_pass_filter(
            HbtepUtilMethods.unwrap_phase(m2_phase), time, 1e-5
        )
        m3_phase_filt = gaussian_low_pass_filter(
            HbtepUtilMethods.unwrap_phase(m3_phase), time, 1e-5
        )
        m4_phase_filt = gaussian_low_pass_filter(
            HbtepUtilMethods.unwrap_phase(m4_phase), time, 1e-5
        )
        m2_freq = np.gradient(m2_phase_filt) / np.gradient(time) / (2 * np.pi)
        m3_freq = np.gradient(m3_phase_filt) / np.gradient(time) / (2 * np.pi)
        m4_freq = np.gradient(m4_phase_filt) / np.gradient(time) / (2 * np.pi)

        # Interpolate output data to the requested timebase
        m2_amp = interp1(time, m2_amp, params.times, "linear")
        m3_amp = interp1(time, m3_amp, params.times, "linear")
        m4_amp = interp1(time, m4_amp, params.times, "linear")
        m2_phase = interp1(time, m2_phase, params.times, "linear")
        m3_phase = interp1(time, m3_phase, params.times, "linear")
        m4_phase = interp1(time, m4_phase, params.times, "linear")
        m2_freq = interp1(time, m2_freq, params.times, "linear")
        m3_freq = interp1(time, m3_freq, params.times, "linear")
        m4_freq = interp1(time, m4_freq, params.times, "linear")

        # Calculate normalized mode amplitudes
        btor = HbtepPhysicsMethods.get_btor(params)["btor"]
        m2_amp_norm = m2_amp / btor
        m3_amp_norm = m3_amp / btor
        m4_amp_norm = m4_amp / btor

        return {
            "m_equal_2_mode": m2_amp,
            "m_equal_2_normalized": m2_amp_norm,
            "m_equal_2_phase": m2_phase,
            "m_equal_2_frequency": m2_freq,
            "m_equal_3_mode": m3_amp,
            "m_equal_3_normalized": m3_amp_norm,
            "m_equal_3_phase": m3_phase,
            "m_equal_3_frequency": m3_freq,
            "m_equal_4_mode": m4_amp,
            "m_equal_4_normalized": m4_amp_norm,
            "m_equal_4_phase": m4_phase,
            "m_equal_4_frequency": m4_freq,
        }

    @staticmethod
    @physics_method(columns=["sxr_total"], tokamak=Tokamak.HBTEP)
    def get_sxr_total(params: PhysicsMethodParams):
        """
        Get the top sxr fan array data and calculate the sum of the 10 good channels
        """
        sxr_data = HbtepPhysicsMethods._get_sxr_data(params)
        sxr_total = np.sum(sxr_data["data"], axis=0)
        sxr_total = interp1(sxr_data["time"], sxr_total, params.times, "linear")
        return {"sxr_total": sxr_total}

    @staticmethod
    @physics_method(
        columns=["euv_000_total", "euv_025_total", "euv_090_total", "euv_270_total"],
        tokamak=Tokamak.HBTEP,
    )
    def get_euv_total(params: PhysicsMethodParams):
        """
        Get the total emissions from the 4 EUV fan arrays
        - 000: outboard midplane array, total emission
        - 025: outboard upper midplane array, emission from top one-third of the plasma
        - 090: top array, emission from the inboard half of the plasma
        - 270: bottom array, emission from the outboard half of the plasma
        """
        euv_data = HbtepPhysicsMethods._get_euv_data(params)
        euv_000_total = interp1(
            euv_data["time"],
            np.sum(euv_data["euv_000_data"], axis=0),
            params.times,
            "linear",
        )
        euv_025_total = interp1(
            euv_data["time"],
            np.sum(euv_data["euv_025_data"], axis=0),
            params.times,
            "linear",
        )
        euv_090_total = interp1(
            euv_data["time"],
            np.sum(euv_data["euv_090_data"], axis=0),
            params.times,
            "linear",
        )
        euv_270_total = interp1(
            euv_data["time"],
            np.sum(euv_data["euv_270_data"], axis=0),
            params.times,
            "linear",
        )

        return {
            "euv_000_total": euv_000_total,
            "euv_025_total": euv_025_total,
            "euv_090_total": euv_090_total,
            "euv_270_total": euv_270_total,
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
            f"ta{i:02}_s{j}p" for i in range(1, 11) for j in range(1, 4)
        ]
        output["ta_rad_names"] = [f"ta{i:02}_s2r" for i in range(1, 11)]
        # Toroidal and poloidal positions of the poloidal sensors
        output["ta_pol_phi"] = np.deg2rad(
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
        output["ta_pol_theta"] = np.deg2rad(
            np.ones(len(output["ta_pol_phi"])) * (189 - 360)
        )
        # Toroidal and poloidal positions of the radial sensors -- not implemented in hbteplib

        # Get TA data and time
        ta_pol_data_raw = []
        ta_pol_data_filt = []
        ta_pol_time = []
        for sensor_name in output["ta_pol_names"]:
            sensor_data_raw, sensor_time = params.mds_conn.get_data_with_dims(
                r"\top.sensors.magnetic:" + sensor_name, tree_name="hbtep2"
            )  # [T], [s]
            fs_ta_pol = round(1 / (sensor_time[1] - sensor_time[0]))
            sensor_data_filt = butterworth_filter(
                sensor_data_raw, fs=fs_ta_pol, cutoff=2e3, order=2, btype="high"
            )
            ta_pol_data_raw.append(sensor_data_raw)
            ta_pol_data_filt.append(sensor_data_filt)
            ta_pol_time.append(sensor_time)
        ta_rad_data_raw = []
        ta_rad_data_filt = []
        ta_rad_time = []
        for sensor_name in output["ta_rad_names"]:
            sensor_data_raw, sensor_time = params.mds_conn.get_data_with_dims(
                r"\top.sensors.magnetic:" + sensor_name, tree_name="hbtep2"
            )  # [T], [s]
            fs_ta_rad = round(1 / (sensor_time[1] - sensor_time[0]))
            sensor_data_filt = butterworth_filter(
                sensor_data_raw, fs=fs_ta_rad, cutoff=2e3, order=2, btype="high"
            )
            ta_rad_data_raw.append(sensor_data_raw)
            ta_rad_data_filt.append(sensor_data_filt)
            ta_rad_time.append(sensor_time)
        output["ta_pol_data_raw"] = ta_pol_data_raw
        output["ta_pol_data_filt"] = ta_pol_data_filt
        output["ta_pol_time"] = ta_pol_time
        output["ta_rad_data_raw"] = ta_rad_data_raw
        output["ta_rad_data_filt"] = ta_rad_data_filt
        output["ta_rad_time"] = ta_rad_time

        return output

    @staticmethod
    @cache_method
    def _get_pa_data(params: PhysicsMethodParams):
        r"""
        Get poloidal field measurements from the 2 poloidal array (PA) sensors
        PA sensors are stored in \HBTEP2::TOP.SENSORS.MAGNETIC:[SENSOR_NAME]
        Return PA1 and PA2 data in their original timebases

        Known bad sensor(s): PA2_S14P, PA2_S27P
        """
        output = {}
        output["bad_sensors"] = ["pa2_s14p", "pa2_s27p"]
        # Sensor names
        output["pa1_names"] = [f"pa1_s{i:02}p" for i in range(1, 33)]
        output["pa2_names"] = [f"pa2_s{i:02}p" for i in range(1, 33)]
        # Toroidal and poloidal positions of the sensors
        output["pa1_theta"] = np.deg2rad(
            [
                -174.74778518,
                -164.23392461,
                -153.66901098,
                -143.01895411,
                -132.24974382,
                -121.3277924,
                -110.22067715,
                -98.93591492,
                -87.23999699,
                -75.60839722,
                -63.97679673,
                -52.34519359,
                -40.71359604,
                -29.08199717,
                -17.45039318,
                -5.81879416,
                5.81280487,
                17.44440438,
                29.07600466,
                40.70760263,
                52.33920936,
                63.97080017,
                75.60240749,
                87.23400093,
                98.93591492,
                110.22067715,
                121.3277924,
                132.24974382,
                143.01895411,
                153.66901098,
                164.23392461,
                174.74778518,
            ]
        )
        output["pa2_theta"] = np.deg2rad(
            [
                -174.74778518,
                -164.23392461,
                -153.66901098,
                -143.01895411,
                -132.24974382,
                -121.3277924,
                -110.22067715,
                -98.93591492,
                -87.23999699,
                -75.60839722,
                -63.97679673,
                -52.34519359,
                -40.71359604,
                -29.08199717,
                -17.45039318,
                -5.81879416,
                5.81280487,
                17.44440438,
                29.07600466,
                40.70760263,
                52.33920936,
                63.97080017,
                75.60240749,
                87.23400093,
                98.93591492,
                110.22067715,
                121.3277924,
                132.24974382,
                143.01895411,
                153.66901098,
                164.23392461,
                174.74778518,
            ]
        )
        output["pa1_phi"] = np.deg2rad(np.ones(len(output["pa1_theta"])) * 317.5)
        output["pa2_phi"] = np.deg2rad(np.ones(len(output["pa2_theta"])) * 317.5)

        # Get PA data and time
        pa1_data_raw = []
        pa1_data_filt = []
        pa1_time = []
        for sensor_name in output["pa1_names"]:
            sensor_data_raw, sensor_time = params.mds_conn.get_data_with_dims(
                r"\top.sensors.magnetic:" + sensor_name, tree_name="hbtep2"
            )
            fs_pa1 = round(1 / (sensor_time[1] - sensor_time[0]))
            sensor_data_filt = butterworth_filter(
                sensor_data_raw, fs=fs_pa1, cutoff=2e3, order=2, btype="high"
            )
            pa1_data_raw.append(sensor_data_raw)
            pa1_data_filt.append(sensor_data_filt)
            pa1_time.append(sensor_time)
        pa2_data_raw = []
        pa2_data_filt = []
        pa2_time = []
        for sensor_name in output["pa2_names"]:
            sensor_data_raw, sensor_time = params.mds_conn.get_data_with_dims(
                r"\top.sensors.magnetic:" + sensor_name, tree_name="hbtep2"
            )
            fs_pa2 = round(1 / (sensor_time[1] - sensor_time[0]))
            sensor_data_filt = butterworth_filter(
                sensor_data_raw, fs=fs_pa2, cutoff=2e3, order=2, btype="high"
            )
            pa2_data_raw.append(sensor_data_raw)
            pa2_data_filt.append(sensor_data_filt)
            pa2_time.append(sensor_time)
        output["pa1_data_raw"] = pa1_data_raw
        output["pa1_data_filt"] = pa1_data_filt
        output["pa1_time"] = pa1_time
        output["pa2_data_raw"] = pa2_data_raw
        output["pa2_data_filt"] = pa2_data_filt
        output["pa2_time"] = pa2_time

        return output

    @staticmethod
    @cache_method
    def _get_sxr_data(params: PhysicsMethodParams):
        """
        Get the top sxr fan array data

        Only 10 (of 16) of the SXR sensors are included in the data below.  Some of the
        missing sensors are broken and others include anomalous or attenuated
        results.
        """
        output = {}
        output["sensor_num"] = [1, 2, 3, 4, 6, 9, 11, 12, 14, 16]
        # channels = output['sensor_num'] + 75

        time = params.mds_conn.get_dims(
            r"\top.sensors.sxr_fan.channel_01:raw", tree_name="hbtep2"
        )
        output["time"] = time[0]

        # Fetch data
        data = []
        r = []
        z = []
        midplane = []
        for n in output["sensor_num"]:
            data.append(
                params.mds_conn.get_data(
                    rf"\top.sensors.sxr_fan.channel_{n:02d}:raw",
                    tree_name="hbtep2",
                )
            )
            r.append(
                params.mds_conn.get_data(
                    rf"\top.sensors.sxr_fan.channel_{n:02d}:r",
                    tree_name="hbtep2",
                )
            )
            z.append(
                params.mds_conn.get_data(
                    rf"\top.sensors.sxr_fan.channel_{n:02d}:z",
                    tree_name="hbtep2",
                )
            )
            midplane.append(
                params.mds_conn.get_data(
                    rf"\top.sensors.sxr_fan.channel_{n:02d}:midplane",
                    tree_name="hbtep2",
                )
            )
        output["data"] = data
        output["r"] = r
        output["z"] = z
        output["midplane"] = midplane

        # Aperture location
        output["det_ap_r"] = 0.9654660224914551
        output["det_ap_z"] = 0.2211250066757202
        output["det_ap_pol"] = 0.001270000007934868
        output["det_ap_tor"] = 0.02539999969303608

        return output

    @staticmethod
    @cache_method
    def _get_euv_data(params: PhysicsMethodParams):
        """
        Get data from the 4 EUV fan array
        """
        output = {}
        output["detectors"] = [0, 25, 90, 270]  # poloidal location of the 4 fan arrays

        time = params.mds_conn.get_dims(
            r"\top.sensors.euv.pol.det000.channel_01:raw", tree_name="hbtep2"
        )
        output["time"] = time[0]

        for detector in output["detectors"]:
            data, r, z, gain = [], [], [], []
            for channel in range(1, 17):
                address = (
                    rf"\top.sensors.euv.pol.det{detector:03d}.channel_{channel:02d}"
                )
                data.append(
                    params.mds_conn.get_data(address + ":raw", tree_name="hbtep2")
                )
                r.append(params.mds_conn.get_data(address + ":r", tree_name="hbtep2"))
                z.append(params.mds_conn.get_data(address + ":z", tree_name="hbtep2"))
                gain.append(
                    params.mds_conn.get_data(address + ":gain", tree_name="hbtep2")
                )
            output[f"euv_{detector:03d}_data"] = data
            output[f"euv_{detector:03d}_r"] = r
            output[f"euv_{detector:03d}_z"] = z
            output[f"euv_{detector:03d}_gain"] = gain

        # Aperture location
        output["det_ap_r"] = [1.160508, 1.090508, 0.907057, 0.929746]
        output["det_ap_z"] = [0.000000, 0.099975, 0.166773, -0.173440]
        output["det_ap_pol"] = [0.000635, 0.000635, 0.000635, 0.000635]
        output["det_ap_tor"] = [0.00635, 0.0254, 0.0254, 0.0244]
        output["orientation"] = [90, 90, 150, 145]  # direction of the aperture [deg]
        # Impact parameters in Tree ordering
        output["impact_parameters"] = [
            -15.20239282,
            -13.88379581,
            -12.34664158,
            -10.57283874,
            -8.55643016,
            -6.31046679,
            -3.87255431,
            -1.30599562,
            1.30599562,
            3.87255431,
            6.31046679,
            8.55643016,
            10.57283874,
            12.34664158,
            13.88379581,
            15.20239282,
            4.30322029,
            5.05782374,
            5.82121675,
            6.58965689,
            7.3590901,
            8.12525048,
            8.88376096,
            9.63024278,
            10.36044694,
            11.07036625,
            11.75634593,
            12.41518118,
            13.04417992,
            13.64121203,
            14.20472504,
            14.73373985,
            0.48160587,
            -0.45909578,
            -1.48806338,
            -2.60489231,
            -3.80409537,
            -5.0735163,
            -6.39312401,
            -7.73487081,
            -9.06416928,
            -10.34320912,
            -11.53564661,
            -12.61149041,
            -13.55074303,
            -14.34481908,
            -14.99565891,
            -15.51324615,
            14.72875124,
            14.0599001,
            13.28190271,
            12.39369218,
            11.39994292,
            10.31173957,
            9.14636266,
            7.92605538,
            6.67595335,
            5.42160618,
            4.18662977,
            2.99095252,
            1.84984336,
            0.7737403,
            -0.23139142,
            -1.16322417,
        ]
        return output
