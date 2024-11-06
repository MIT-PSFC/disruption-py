#!/usr/bin/env python3

"""
Module for generic physics methods.
"""

from typing import Dict

import numpy as np

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
from disruption_py.machine.tokamak import Tokamak


class GenericPhysicsMethods:
    """
    Class to hold generic physics methods.
    """

    @staticmethod
    def _config_kappa_area(params: PhysicsMethodParams) -> Dict[str, str]:
        """
        Input configuration method for `kappa_area`.
        """

        defaults = {
            "aminor": r"\efit_a_eqdsk:aminor",
            "area": r"\efit_a_eqdsk:area",
            "atime": r"\efit_a_eqdsk:time",
        }

        overrides = {
            Tokamak.D3D: {
                "atime": r"\efit_a_eqdsk:atime/1000",
            },
        }

        config = defaults.copy()
        config.update(overrides.get(params.tokamak, {}))

        return config

    @staticmethod
    def _inputs_kappa_area(
        params: PhysicsMethodParams, config: Dict[str, str]
    ) -> Dict[str, np.ndarray]:
        """
        Data fetching method for `kappa_area`.
        """

        inputs = {
            k: params.mds_conn.get_data(v, tree_name="_efit_tree", astype="float64")
            for k, v in config.items()
        }

        return inputs

    @staticmethod
    def _output_kappa_area(
        times: np.ndarray, aminor: np.ndarray, area: np.ndarray, atime: np.ndarray
    ) -> Dict[str, np.ndarray]:
        """
        Computation method for `kappa_area`.
        """

        # filter
        aminor[aminor == 0] = np.nan
        area[area == 0] = np.nan

        # interpolate
        kappa_area = interp1(x=atime, y=area / (np.pi * aminor**2), new_x=times)

        # output
        outputs = {"kappa_area": kappa_area}
        return outputs

    @staticmethod
    @physics_method(columns=["kappa_area"])
    def get_kappa_area(params: PhysicsMethodParams) -> Dict[str, np.ndarray]:
        """
        # TODO: typeset mathematical formula
        Compute `kappa_area`, the elongation parameter, defined as:
        plasma area / (pi * aminor**2)

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/CMOD/matlab-core/get_kappa_area.m
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/get_kappa_area.m
        """

        config = GenericPhysicsMethods._config_kappa_area(params)
        inputs = GenericPhysicsMethods._inputs_kappa_area(params, config)
        output = GenericPhysicsMethods._output_kappa_area(times=params.times, **inputs)

        return output
