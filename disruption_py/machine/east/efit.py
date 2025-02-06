#!/usr/bin/env python3

"""
Module for retrieving and processing EFIT parameters for EAST.
"""

import numpy as np

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
from disruption_py.machine.tokamak import Tokamak


class EastEfitMethods:
    """
    Class for retrieving and processing EFIT parameters for EAST.

    Attributes
    ----------
    efit_cols : dict
        A dictionary mapping parameter names to their corresponding EFIT data paths.
    """

    # TODO: confirm units
    efit_cols = {
        "beta_n": r"\efit_aeqdsk:betan",
        "beta_p": r"\efit_aeqdsk:betap",
        "kappa": r"\efit_aeqdsk:kappa",
        "li": r"\efit_aeqdsk:li",
        "q0": r"\efit_aeqdsk:q0",
        "qstar": r"\efit_aeqdsk:qstar",
        "q95": r"\efit_aeqdsk:q95",
        "wmhd": r"\efit_aeqdsk:wmhd",
        "area": r"\efit_aeqdsk:area*1e-4",  # [cm^2] -> [m^2], from get_kappa_area
        "aminor": r"\efit_aeqdsk:aout",
        "chisq": r"\efit_aeqdsk:chisq",
    }
    efit_derivs = {"dbetap_dt": "beta_p", "dli_dt": "li", "dwmhd_dt": "wmhd"}

    pefit_cols = {
        "pbeta_n": r"\efit_aeqdsk:betan",
        "pbeta_p": r"\efit_aeqdsk:betap",
        "pkappa": r"\efit_aeqdsk:kappa",
        "pli": r"\efit_aeqdsk:li",
        "pq0": r"\efit_aeqdsk:q0",
        "pqstar": r"\efit_aeqdsk:qstar",
        "pq95": r"\efit_aeqdsk:q95",
        "pwmhd": r"\efit_aeqdsk:wmhd",
        "parea": r"\efit_aeqdsk:area*1e-4",  # [cm^2] -> [m^2]
        "paminor": r"\efit_aeqdsk:aout",
        "pchisq": r"\efit_aeqdsk:chisq",
        "pconvergence": r"\efit_aeqdsk:error",
    }

    @staticmethod
    @physics_method(
        columns=[*efit_cols.keys(), *efit_derivs.keys()],
        tokamak=Tokamak.EAST,
    )
    def get_efit_parameters(params: PhysicsMethodParams):
        """
        Retrieve EFIT parameters for EAST.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDS connection and shot information.

        Returns
        -------
        dict
            A dictionary containing the retrieved EFIT parameters.
        """
        efit_data = {
            k: params.mds_conn.get_data(v, tree_name="_efit_tree")
            for k, v in EastEfitMethods.efit_cols.items()
        }
        efit_time = params.mds_conn.get_data(
            r"\efit_aeqdsk:atime", tree_name="_efit_tree", astype="float64"
        )  # TODO: [unit?]

        # EFIT reconstructions are sometimes invalid, particularly when very close
        # to a disruption.  There are a number of EFIT parameters that can indicate
        # invalid reconstructions, such as 'terror' and 'chisq'.  Here we use
        # 'chisq' to determine which time slices should be excluded from our
        # disruption warning database.
        invalid_indices = np.where(efit_data["chisq"] > 50)

        for param in efit_data:
            efit_data[param][invalid_indices] = np.nan
        # pylint: disable=duplicate-code
        for deriv_param, param in EastEfitMethods.efit_derivs.items():
            efit_data[deriv_param] = np.gradient(efit_data[param], efit_time)
        if not np.array_equal(params.times, efit_time):
            for param in efit_data:
                efit_data[param] = interp1(efit_time, efit_data[param], params.times)
        return efit_data

    @staticmethod
    @physics_method(
        columns=[
            *pefit_cols.keys(),
        ],
        tokamak=Tokamak.EAST,
    )
    def get_pefit_parameters(params: PhysicsMethodParams):
        """
        Retrieve real-time P-EFIT (Parallel-EFIT) parameters from
        the 'pefit_east' tree.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDS connection and shot information.

        Returns
        -------
        dict
            A dictionary containing the retrieved P-EFIT parameters.
        """
        efit_time = params.mds_conn.get_data(
            r"\efit_a_eqdsk:atime", tree_name="pefit_east"
        )  # TODO: [unit?]

        # Deal with bug
        efit_time, unique_indices = np.unique(efit_time, return_index=True)

        efit_data = {
            k: params.mds_conn.get_data(v, tree_name="pefit_east")[unique_indices]
            for k, v in EastEfitMethods.pefit_cols.items()
        }

        # P-EFIT reconstructions are sometimes invalid, particularly when very close
        # to a disruption.  There are a number of P-EFIT parameters that can indicate
        # invalid reconstructions, such as 'error' and 'chisq'.  Here we use
        # 'chisq' to determine which time slices should be excluded from our
        # disruption warning database.
        # Yao Huang suggests to use data with:
        #   - chisq < 20
        #   - convergence error < 1 (data stored in MDS+ is multiplied by 1e3)
        #   - ip > 180 kA
        # For now, we only check the first two conditions.
        # If ever we want to extend analysis to ramp up or down we need to check ip.
        invalid_indices = np.where(
            (efit_data["pchisq"] > 50) & (efit_data["pconvergence"] < 1)
        )
        # pylint: disable=duplicate-code
        for param in efit_data:
            efit_data[param][invalid_indices] = np.nan
        if not np.array_equal(params.times, efit_time):
            for param in efit_data:
                efit_data[param] = interp1(efit_time, efit_data[param], params.times)
        return efit_data
