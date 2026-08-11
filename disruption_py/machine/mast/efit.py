#!/usr/bin/env python3

"""
Module for retrieving and processing EFIT parameters for MAST.
"""

import numpy as np

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.errors import DataError
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.machine.mast.util import MastUtilMethods
from disruption_py.machine.tokamak import Tokamak


class MastEfitMethods:
    """
    Class for retrieving and processing EFIT parameters for MAST.
    """

    efit_properties = {
        "a_minor": "minor_radius",
        "beta_n": "beta_tor_normal",
        "beta_p": "beta_pol",
        "beta_t": "beta_tor",
        "bphi_rmag": "bphi_rmag",
        "bvac_rmag": "bvac_rmag",
        "kappa": "elongation",
        "li": "li",
        "q95": "q95",
        "rmagx": "magnetic_axis_r",
        "rmagz": "magnetic_axis_z",
        "tribot": "triangularity_lower",
        "tritop": "triangularity_upper",
        "volume": "volume",
        "v_loop_dynamic": "vloop_dynamic",
        "v_loop_static": "vloop_static",
        "wmhd": "wmhd",
    }

    @staticmethod
    @physics_method(columns=list(efit_properties.keys()), tokamak=Tokamak.MAST)
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
            A dictionary containing the retrieved EFIT parameters. Properties whose
            underlying signal is missing for this shot are returned as NaN rather
            than dropping the whole equilibrium reconstruction.
        """
        times = params.times
        try:
            eq_time = params.get_data("equilibrium/time", required=True)
        except DataError:
            return {key: [np.nan] for key in MastEfitMethods.efit_properties}

        outputs = {}
        for key, prop in MastEfitMethods.efit_properties.items():
            signal = params.get_data(f"equilibrium/{prop}")
            outputs[key] = MastUtilMethods.interpolate_1d(eq_time, signal, times)

        return outputs
