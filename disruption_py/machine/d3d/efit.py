#!/usr/bin/env python3

import numpy as np

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
from disruption_py.machine.tokamak import Tokamak


class D3DEfitMethods:

    # EFIT Variables
    efit_cols = {
        "beta_n": r"\efit_a_eqdsk:betan",
        "beta_p": r"\efit_a_eqdsk:betap",
        "kappa": r"\efit_a_eqdsk:kappa",
        "li": r"\efit_a_eqdsk:li",
        "upper_gap": r"\efit_a_eqdsk:gaptop",
        "lower_gap": r"\efit_a_eqdsk:gapbot",
        "q0": r"\efit_a_eqdsk:q0",
        "qstar": r"\efit_a_eqdsk:qstar",
        "q95": r"\efit_a_eqdsk:q95",
        "wmhd": r"\efit_a_eqdsk:wmhd",
        "chisq": r"\efit_a_eqdsk:chisq",
    }
    # 'v_loop_efit': ,r'\efit_a_eqdsk:vsurf', 'bt0': r'\efit_a_eqdsk:bt0'
    efit_derivs = {"dbetap_dt": "beta_p", "dli_dt": "li", "dwmhd_dt": "wmhd"}
    rt_efit_cols = {
        "beta_p_rt": r"\efit_a_eqdsk:betap",
        "li_rt": r"\efit_a_eqdsk:li",
        "q95_rt": r"\efit_a_eqdsk:q95",
        "wmhd_rt": r"\efit_a_eqdsk:wmhd",
        "chisq": r"\efit_a_eqdsk:chisq",
    }
    # 'v_loop_efit_RT': r'\efit_a_eqdsk:vsurf',

    @staticmethod
    @physics_method(
        columns=[*efit_cols.keys(), *efit_derivs.keys()],
        tokamak=Tokamak.D3D,
    )
    def _get_efit_parameters(params: PhysicsMethodParams):
        efit_data = {
            k: params.mds_conn.get_data(v, tree_name="_efit_tree")
            for k, v in D3DEfitMethods.efit_cols.items()
        }
        efit_time = (
            params.mds_conn.get_data(r"\efit_a_eqdsk:atime", tree_name="_efit_tree")
            / 1.0e3
        )  # [ms] -> [s]

        # EFIT reconstructions are sometimes invalid, particularly when very close
        # to a disruption.  There are a number of EFIT parameters that can indicate
        # invalid reconstructions, such as 'terror' and 'chisq'.  Here we use
        # 'chisq' to determine which time slices should be excluded from our
        # disruption warning database.
        invalid_indices = np.where(efit_data["chisq"] > 50)
        del efit_data["chisq"]
        for param in efit_data:
            efit_data[param][invalid_indices] = np.nan
        for deriv_param in D3DEfitMethods.efit_derivs:
            efit_data[deriv_param] = np.gradient(
                efit_data[D3DEfitMethods.efit_derivs[deriv_param]], efit_time
            )
        if not np.array_equal(params.times, efit_time):
            for param in efit_data:
                efit_data[param] = interp1(efit_time, efit_data[param], params.times)
        return efit_data

    @staticmethod
    @physics_method(
        columns=[*rt_efit_cols.keys()],
        tokamak=Tokamak.D3D,
    )
    def _get_rt_efit_parameters(params: PhysicsMethodParams):
        efit_data = {
            k: params.mds_conn.get_data(v, tree_name="efitrt1")
            for k, v in D3DEfitMethods.rt_efit_cols.items()
        }
        efit_time = (
            params.mds_conn.get_data(r"\efit_a_eqdsk:atime", tree_name="efitrt1")
            / 1.0e3
        )  # [ms] -> [s]
        # EFIT reconstructions are sometimes invalid, particularly when very close
        # to a disruption.  There are a number of EFIT parameters that can indicate
        # invalid reconstructions, such as 'terror' and 'chisq'.  Here we use
        # 'chisq' to determine which time slices should be excluded from our
        # disruption warning database.
        invalid_indices = np.where(efit_data["chisq"] > 50)
        del efit_data["chisq"]
        for param in efit_data:
            efit_data[param][invalid_indices] = np.nan
        if not np.array_equal(params.times, efit_time):
            for param in efit_data:
                efit_data[param] = interp1(efit_time, efit_data[param], params.times)
        return efit_data
