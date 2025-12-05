"""
Module for retrieving and processing EFIT parameters for CMOD.
"""

from disruption_py.core.utils.math import interp1
import xarray as xr

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.inout.xarray import XarrayConnection


class MastEfitMethods:
    """
    Class for retrieving and processing EFIT parameters for MAST.
    """

    efit_properties = [
        "beta_tor_normal",
        "elongation",
        "elongation_axis",
        "magnetic_axis_r",
        "magnetic_axis_z",
        "triangularity_lower",
        "triangularity_upper",
        "minor_radius",
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
            item = interp1(eq_time, signal, times)
            outputs[prop] = item

        return outputs
