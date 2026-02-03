#!/usr/bin/env python3

"""
Module for retrieving and processing EFIT parameters for MAST.
"""

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.inout.xr import XarrayConnection
from disruption_py.machine.mast.util import MastUtilMethods


class MastEfitMethods:
    """
    Class for retrieving and processing EFIT parameters for MAST.
    """

    efit_properties = {
        "beta_n": "beta_tor_normal",
        "beta_t": "beta_tor",
        "beta_p": "beta_pol",
        "kappa": "elongation",
        "rmagx": "magnetic_axis_r",
        "rmagz": "magnetic_axis_z",
        "tribot": "triangularity_lower",
        "tritop": "triangularity_upper",
        "a_minor": "minor_radius",
        "bvac_rmag": "bvac_rmag",
        "bphi_rmag": "bphi_rmag",
        "li": "li",
        "v_loop_static": "vloop_static",
        "v_loop_dynamic": "vloop_dynamic",
        "q95": "q95",
    }

    @staticmethod
    @physics_method(columns=list(efit_properties.keys()))
    def get_efit_parameters(params: PhysicsMethodParams):
        """
        Retrieve EFIT parameters for MAST.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing theconnection and shot information.

        Returns
        -------
        dict
            A dictionary containing the retrieved EFIT parameters.
        """
        conn: XarrayConnection = params.mds_conn
        eq_time = conn.get_data(params.shot_id, "equilibrium/time")
        times = params.times

        outputs = {}
        for key, prop in MastEfitMethods.efit_properties.items():
            signal = conn.get_data(params.shot_id, f"equilibrium/{prop}")
            item = MastUtilMethods.interpolate_1d(eq_time, signal, times)
            outputs[key] = item

        return outputs
