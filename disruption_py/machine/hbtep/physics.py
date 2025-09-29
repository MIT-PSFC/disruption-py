#!/usr/bin/env python3

"""
Module for retrieving and calculating data for HBTEP physics methods.
"""

import numpy as np
from scipy.integrate import cumulative_trapezoid

from disruption_py.core.physics_method.caching import cache_method
from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
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
