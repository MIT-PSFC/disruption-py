#!/usr/bin/env python3

import traceback

import numpy as np

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
from disruption_py.machine.tokamak import Tokamak


class CmodEfitMethods:

    efit_cols = {
        "beta_n": r"\efit_aeqdsk:betan",
        "beta_p": r"\efit_aeqdsk:betap",
        "kappa": r"\efit_aeqdsk:eout",
        "li": r"\efit_aeqdsk:li",
        "upper_gap": r"\efit_aeqdsk:otop/100",
        "lower_gap": r"\efit_aeqdsk:obott/100",
        "q0": r"\efit_aeqdsk:q0",
        "qstar": r"\efit_aeqdsk:qstar",
        "q95": r"\efit_aeqdsk:q95",
        "v_loop_efit": r"\efit_aeqdsk:vloopt",
        "wmhd": r"\efit_aeqdsk:wplasm",
        "ssep": r"\efit_aeqdsk:ssep/100",
        "n_over_ncrit": r"-\efit_aeqdsk:xnnc",
        "tritop": r"\efit_aeqdsk:doutu",
        "tribot": r"\efit_aeqdsk:doutl",
        "a_minor": r"\efit_aeqdsk:aminor",
        "rmagx": r"\efit_aeqdsk:rmagx",  # TODO: change units to [m] (current [cm])
        "chisq": r"\efit_aeqdsk:chisq",
    }

    # EFIT column names for data before 2000 TODO: confirm with Bob that these are
    # the right back-ups and make sure that these are similar to standard EFIT columns
    efit_cols_pre_2000 = {
        "a_minor": r"\efit_aeqdsk:aout",
        "li": r"\efit_aeqdsk:ali",
        "q0": r"\efit_aeqdsk:qqmagx",
        "qstar": r"\efit_aeqdsk:qsta",
        "q95": r"\efit_aeqdsk:qsib",  # Not sure about this one
    }

    efit_derivs = {"dbetap_dt": "beta_p", "dli_dt": "li", "dwmhd_dt": "wmhd"}

    @staticmethod
    @physics_method(
        columns=[
            *efit_cols.keys(),
            *efit_cols_pre_2000.keys(),
            *efit_derivs.keys(),
            "v_surf",
            "v_loop_efit",
            "beta_n",
        ],
        tokamak=Tokamak.CMOD,
    )
    def _get_EFIT_parameters(params: PhysicsMethodParams):

        efit_time = params.mds_conn.get_data(
            r"\efit_aeqdsk:time", tree_name="_efit_tree", astype="float64"
        )  # [s]
        efit_data = dict()

        # Get data from each of the columns in efit_cols one at a time
        for param in CmodEfitMethods.efit_cols:
            try:
                # If shot before 2000 and the param is in efit_cols_pre_2000
                if (
                    params.shot_id <= 1000000000
                    and param not in CmodEfitMethods.efit_cols_pre_2000.keys()
                ):
                    efit_data[param] = params.mds_conn.get_data(
                        CmodEfitMethods.efit_cols_pre_2000[param],
                        tree_name="_efit_tree",
                        astype="float64",
                    )
                else:
                    efit_data[param] = params.mds_conn.get_data(
                        CmodEfitMethods.efit_cols[param],
                        tree_name="_efit_tree",
                        astype="float64",
                    )
            except:
                params.logger.warning(
                    f"[Shot {params.shot_id}]: Unable to get {param} from EFIT tree"
                )
                params.logger.debug(
                    f"[Shot {params.shot_id}]: {traceback.format_exc()}"
                )
                efit_data[param] = np.full(len(efit_time), np.nan)
                pass

        for deriv_param, param in CmodEfitMethods.efit_derivs.items():
            efit_data[deriv_param] = np.gradient(
                efit_data[param],
                efit_time,
                edge_order=1,
            )

        # Get data for V_surf := deriv(\ANALYSIS::EFIT_SSIBRY)*2*pi
        try:
            ssibry = params.mds_conn.get_data(
                r"\efit_geqdsk:ssibry", tree_name="_efit_tree", astype="float64"
            )
            efit_data["v_surf"] = np.gradient(ssibry, efit_time) * 2 * np.pi
        except:
            print("unable to get V_surf")
            efit_data["v_surf"] = np.full(len(efit_time), np.nan)
            pass

        # For shots before 2000, adjust units of aminor, compute beta_n and v_loop
        if params.shot_id <= 1000000000:

            # Adjust aminor units
            efit_data["aminor"] = efit_data["aminor"] / 100  # [cm] to [m]

            # Get data for v_loop --> deriv(\ANALYSIS::EFIT_SSIMAG)*$2pi (not totally
            #  sure on this one)
            try:  # TODO: confirm this
                ssimag = params.mds_conn.get_data(
                    r"\efit_geqdsk:ssimag", tree_name="_efit_tree", astype="float64"
                )
                efit_data["v_loop_efit"] = np.gradient(ssimag, efit_time) * 2 * np.pi
            except:
                print("unable to get v_loop_efit")
                efit_data["v_loop_efit"] = np.full(len(efit_time), np.nan)
                pass

            # Compute beta_n
            beta_t = params.mds_conn.get_data(
                r"\efit_aeqdsk:betat", tree_name="_efit_tree", astype="float64"
            )
            efit_data["beta_n"] = np.reciprocal(
                np.reciprocal(beta_t) + np.reciprocal(efit_data["beta_p"])
            )

        if not np.array_equal(params.times, efit_time):
            for param in efit_data:
                efit_data[param] = interp1(efit_time, efit_data[param], params.times)

        return efit_data

    @staticmethod
    def efit_check(params: PhysicsMethodParams):
        """
        # TODO: Get description from Jinxiang
        """
        values = []
        for expr in [
            r"_lf=\analysis::efit_aeqdsk:lflag",
            r"_l0=((sum(_lf,1) - _lf[*,20] - _lf[*,1])==0)",
            r"_n=\analysis::efit_fitout:nitera,(_l0 and (_n>4))",
        ]:
            values.append(params.mds_conn.get(expr, tree_name="analysis"))
        _n = values[2].data()
        valid_indices = np.nonzero(_n)
        (times,) = params.mds_conn.get_dims(
            r"\analysis::efit_aeqdsk:lflag", tree_name="analysis"
        )
        return valid_indices, times[valid_indices]
