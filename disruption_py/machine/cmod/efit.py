#!/usr/bin/env python3

"""
Module for retrieving and processing EFIT parameters for CMOD.
"""

import numpy as np

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
from disruption_py.inout.mds import mdsExceptions
from disruption_py.machine.generic.efit import GenericEFITMethods
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
    geqdsk_cols : dict
        A dictionary mapping parameter names to their corresponding EFIT data paths
    archive_dir : str
        Directory where 1 kHz EFIT data is archived.
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
        "rout": r"\efit_aeqdsk:rout/100",
    }

    efit_derivs = {"dbetap_dt": "beta_p", "dli_dt": "li", "dwmhd_dt": "wmhd"}

    geqdsk_cols = {
        # 0D signals
        "rmagx": r"\efit_g_eqdsk:rmaxis",
        "zmagx": r"\efit_g_eqdsk:zmaxis",
        "simagx": r"\efit_g_eqdsk:ssimag",
        "sibdry": r"\efit_g_eqdsk:ssibry",
        "bcentr": r"\efit_g_eqdsk:bcentr",
        "current": r"\efit_g_eqdsk:cpasma",
        # 1D profiles
        "fpol": r"\efit_g_eqdsk:fpol",
        "pres": r"\efit_g_eqdsk:pres",
        "ffprime": r"\efit_g_eqdsk:ffprim",
        "pprime": r"\efit_g_eqdsk:pprime",
        "qpsi": r"\efit_g_eqdsk:qpsi",
        # 2D flux grid
        "psirz": r"\efit_g_eqdsk:psirz",
        # Boundary
        "rbdry": r"\efit_g_eqdsk:rbbbs",
        "zbdry": r"\efit_g_eqdsk:zbbbs",
    }

    archive_dir = "/usr/local/mfe/disruptions/data/disruption-efit"

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
            The parameters containing the data connection and shot information.

        Returns
        -------
        dict
            A dictionary containing the retrieved EFIT parameters.
        """
        efit_time = params.get_data(r"\efit_aeqdsk:time", tree_name="_efit_tree")  # [s]
        efit_data = {}

        # Get data from each of the columns in efit_cols one at a time
        for param, path in CmodEfitMethods.efit_cols.items():
            try:
                efit_data[param] = params.get_data(
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
            The parameters containing the data connection and shot information.

        Returns
        -------
        tuple
            A tuple containing valid indices and corresponding times.
        """
        values = [
            params.data_conn.get(expr, tree_name="analysis")
            for expr in [
                r"_lf=\efit_aeqdsk:lflag",
                r"_l0=((sum(_lf,1) - _lf[*,20] - _lf[*,1])==0)",
                r"_n=\efit_fitout:nitera,(_l0 and (_n>4))",
            ]
        ]
        _n = values[2].data()
        valid_indices = np.nonzero(_n)
        (times,) = params.get_dims(r"\efit_aeqdsk:lflag", tree_name="analysis")
        return valid_indices, times[valid_indices]

    @staticmethod
    @physics_method(
        columns=[*geqdsk_cols.keys()],
        tokamak=Tokamak.CMOD,
    )
    def get_geqdsk_parameters(params: PhysicsMethodParams):
        """
        Create an Xarray dataset containing all signals necessary to recreate a GEQDSK file.
        Uses FreeQDSK canonical names. 2D psi grid stored as 'psirz' to distinguish from 'psi' coordinate.
        Data is in COCOS 1.
        """
        efit_time = params.get_data(r"\efit_a_eqdsk:atime", tree_name="_efit_tree")
        r_grid = params.get_data(r"\efit_g_eqdsk:rgrid", tree_name="_efit_tree")
        z_grid = params.get_data(r"\efit_g_eqdsk:zgrid", tree_name="_efit_tree")

        geqdsk_data = {}
        for param, path in CmodEfitMethods.geqdsk_cols.items():
            try:
                geqdsk_data[param] = params.get_data(path, tree_name="_efit_tree")
            except mdsExceptions.MdsException as e:
                params.logger.warning(repr(e))
                params.logger.opt(exception=True).debug(e)
                geqdsk_data[param] = None

        rlim = zlim = None
        try:
            rlim = params.get_data(r"\efit_g_eqdsk:xlim", tree_name="_efit_tree")
            zlim = params.get_data(r"\efit_g_eqdsk:ylim", tree_name="_efit_tree")
        except mdsExceptions.MdsException as e:
            params.logger.warning(repr(e))
            params.logger.opt(exception=True).debug(e)

        # MDS storage layouts (from C-Mod EFIT):
        #   0D scalars: (T,)
        #   1D profiles (fpol, pres, ffprime, pprime, qpsi): (psi_idx, T)
        #   2D psirz: (T, z, r)
        #   boundary (rbdry, zbdry): (boundary_idx, T)
        # make_geqdsk_dataset expects time on axis 0 for all signals.
        time_last_params = {"fpol", "pres", "ffprime", "pprime", "qpsi", "rbdry", "zbdry"}
        needs_interp = not np.array_equal(params.times, efit_time)
        for param in geqdsk_data:
            if geqdsk_data[param] is None:
                continue
            data = geqdsk_data[param]
            if data.ndim == 1:
                if needs_interp:
                    geqdsk_data[param] = interp1(efit_time, data, params.times)
            elif param in time_last_params:
                # (spatial_idx, T) -> interp last axis -> (spatial_idx, T') -> T
                if needs_interp:
                    geqdsk_data[param] = interp1(efit_time, data, params.times).T
                else:
                    geqdsk_data[param] = data.T
            else:
                # psirz: (T, z, r) - time already first, interp along axis 0
                if needs_interp:
                    geqdsk_data[param] = interp1(efit_time, data, params.times, axis=0)


        # This is EFIT, determine COCOS based on the sign of median Ip and B0 https://efit-ai.gitlab.io/efit/files.html
        sign_Ip = np.sign(np.median(geqdsk_data["current"]))
        sign_B0 = np.sign(np.median(geqdsk_data["bcentr"]))
        if sign_Ip > 0 and sign_B0 > 0:
            cocos_input = 1
        elif sign_Ip < 0 and sign_B0 > 0:
            cocos_input = 3
        elif sign_Ip > 0 and sign_B0 < 0:
            cocos_input = 5
        elif sign_Ip < 0 and sign_B0 < 0:
            cocos_input = 7
        else:
            params.logger.warning(
                "Unexpected sign combination for current and magnetic field. Assuming COCOS 1."
            )
            cocos_input = 1

        ds = GenericEFITMethods.make_geqdsk_dataset(
            shot_id=params.shot_id,
            times=params.times,
            r_grid=r_grid,
            z_grid=z_grid,
            rmagx=geqdsk_data["rmagx"],
            zmagx=geqdsk_data["zmagx"],
            simagx=geqdsk_data["simagx"],
            sibdry=geqdsk_data["sibdry"],
            bcentr=geqdsk_data["bcentr"],
            current=geqdsk_data["current"],
            fpol=geqdsk_data["fpol"],
            pres=geqdsk_data["pres"],
            ffprime=geqdsk_data["ffprime"],
            pprime=geqdsk_data["pprime"],
            qpsi=geqdsk_data["qpsi"],
            psirz=geqdsk_data["psirz"],
            rbdry=geqdsk_data["rbdry"],
            zbdry=geqdsk_data["zbdry"],
            cocos_input=cocos_input,
            rlim=rlim,
            zlim=zlim,
        )

        return ds