#!/usr/bin/env python3

"""
Module for retrieving and processing EFIT parameters for CMOD.
"""

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.inout.xarray import XarrayConnection
from disruption_py.machine.mast.physics import MastPhysicsMethods


class MastEfitMethods:
    """
    Class for retrieving and processing EFIT parameters for MAST.
    """

    efit_properties = [
        "beta_tor_normal",
        "beta_tor",
        "beta_pol",
        "elongation",
        "elongation_axis",
        "magnetic_axis_r",
        "magnetic_axis_z",
        "triangularity_lower",
        "triangularity_upper",
        "minor_radius",
        "bvac_rmag",
        "bphi_rmag",
        "li",
        "whmd",
        "vloop_static",
        "q95",
    ]

    @staticmethod
    @physics_method(columns=efit_properties)
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
        for prop in MastEfitMethods.efit_properties:
            signal = conn.get_data(params.shot_id, f"equilibrium/{prop}")
            item = MastPhysicsMethods.interpolate_1d(eq_time, signal, times)
            outputs[prop] = item

        if "whmd" in outputs:
            outputs["wmhd"] = outputs.pop("whmd")

        return outputs
