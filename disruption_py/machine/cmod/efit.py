#!/usr/bin/env python3

"""
Module for retrieving and processing EFIT parameters for CMOD.
"""

import numpy as np
from MDSplus import mdsExceptions

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
from disruption_py.machine.tokamak import Tokamak


class CmodEfitMethods:
    """
    Class for retrieving and processing EFIT parameters for CMOD.

    Attributes
    ----------
    efit_cols : dict
        A dictionary mapping parameter names to their corresponding EFIT data paths.
    efit_derivs : dict
        A dictionary mapping derivative parameter names to their corresponding base parameters.
    """

    efit_cols = {
        "beta_p": r"\efit_aeqdsk:betap",
        "kappa": r"\efit_aeqdsk:eout",
        "li": r"\efit_aeqdsk:ali",
        "upper_gap": r"\efit_aeqdsk:otop/100",
        "lower_gap": r"\efit_aeqdsk:obott/100",
        "q0": r"\efit_aeqdsk:qqmagx",
        "qstar": r"\efit_aeqdsk:qsta",
        "q95": r"\efit_aeqdsk:qpsib",
        "v_loop_efit": r"\efit_aeqdsk:vloopt",
        "wmhd": r"\efit_aeqdsk:wplasm",
        "ssep": r"\efit_aeqdsk:ssep/100",
        "n_over_ncrit": r"-\efit_aeqdsk:xnnc",
        "tritop": r"\efit_aeqdsk:doutu",
        "tribot": r"\efit_aeqdsk:doutl",
        "a_minor": r"\efit_aeqdsk:aout/100",
        "rmagx": r"\efit_aeqdsk:rmagx/100",
        "chisq": r"\efit_aeqdsk:tsaisq",
        # additional signals for defuse demo
        "area": r"\efit_aeqdsk:areao/1e4",
        "delta": r"(\efit_aeqdsk:doutu + \efit_aeqdsk:doutl)/2",
        "rout": r"\efit_aeqdsk:rout/100",
        "volume": r"\efit_aeqdsk:volume",
        "rcur": r"\efit_aeqdsk:rcur",
        "zmagx": r"\efit_aeqdsk:zmagx/100",
    }

    efit_derivs = {"dbetap_dt": "beta_p", "dli_dt": "li", "dwmhd_dt": "wmhd"}

    @staticmethod
    @physics_method(
        columns=[
            *efit_cols.keys(),
            *efit_derivs.keys(),
        ],
        tokamak=Tokamak.CMOD,
    )
    def get_efit_parameters(params: PhysicsMethodParams):
        """
        Retrieve EFIT parameters for CMOD.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDS connection and shot information.

        Returns
        -------
        dict
            A dictionary containing the retrieved EFIT parameters.
        """
        efit_time = params.mds_conn.get_data(
            r"\efit_aeqdsk:time", tree_name="_efit_tree"
        )  # [s]
        efit_data = {}

        # Get data from each of the columns in efit_cols one at a time
        for param, path in CmodEfitMethods.efit_cols.items():
            try:
                efit_data[param] = params.mds_conn.get_data(
                    path=path,
                    tree_name="_efit_tree",
                )
            except mdsExceptions.MdsException as e:
                params.logger.warning(repr(e))
                params.logger.opt(exception=True).debug(e)
                efit_data[param] = np.full(len(efit_time), np.nan)

        for deriv_param, param in CmodEfitMethods.efit_derivs.items():
            efit_data[deriv_param] = np.gradient(
                efit_data[param],
                efit_time,
                edge_order=1,
            )

        if not np.array_equal(params.times, efit_time):
            for param in efit_data:
                efit_data[param] = interp1(efit_time, efit_data[param], params.times)

        return efit_data

    @staticmethod
    def efit_check(params: PhysicsMethodParams):
        """
        Check the validity of EFIT parameters for the given shot.
        # TODO: Get description from Jinxiang

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDS connection and shot information.

        Returns
        -------
        tuple
            A tuple containing valid indices and corresponding times.
        """
        values = [
            params.mds_conn.get(expr, tree_name="analysis")
            for expr in [
                r"_lf=\efit_aeqdsk:lflag",
                r"_l0=((sum(_lf,1) - _lf[*,20] - _lf[*,1])==0)",
                r"_n=\efit_fitout:nitera,(_l0 and (_n>4))",
            ]
        ]
        _n = values[2].data()
        valid_indices = np.nonzero(_n)
        (times,) = params.mds_conn.get_dims(r"\efit_aeqdsk:lflag", tree_name="analysis")
        return valid_indices, times[valid_indices]
