#!/usr/bin/env python3

"""
Module for retrieving and processing EFIT parameters for CMOD.
"""

import numpy as np
import xarray as xr

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
from disruption_py.inout.mds import mdsExceptions
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
    @physics_method(
        tokamak=Tokamak.CMOD,
    )
    def get_geqdsk_parameters(params: PhysicsMethodParams):
        """
        Retrieve parameters required to create a GEQDSKFile https://freeqdsk.readthedocs.io/en/stable/geqdsk.html
        Also see https://fusion.gat.com/theory/Efitgeqdsk
        This requires the following:
        
        1. shot number
        2. grid parameters: (nx, ny, rdim, zdim, rcenter, rleft, zmid) or (rgrid, zgrid)
        3. scalar parameters: rmagx, zmagx, simagx, sibdry, bccenter, cpasma, etc.
        4. 1D profile parameters: fpol, pres, ffprime, pprime, qpsi
        5. 2D parameters: psi or f or psirz
        6. boundary and limiter data (if available): nbdry, nlim, rbdry, zbdry, rlim, zlim

        Sometimes the C-Mod EFIT uses 'x' and 'y' instead of 'r' and 'z' for parameter names
        This function will convert them to the standard GEQDSK naming convention.

        Also, it doesn't make sense to get interpolated boundary/limiter data because the number of points can change over time.
        TODO(ZanderKeith): How should we handle this case?

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information.
        Returns
        -------
        xarray.Dataset
            A dataset containing the GEQDSK parameters.
        """
        # Get EFIT time array
        efit_time = params.mds_conn.get_data(r"\efit_a_eqdsk:atime", tree_name="_efit_tree")

        # Grid parameters (For C-Mod a lot are missing, but can be inferred from rgrid/zgrid)
        grid_params = {
            "nw": params.mds_conn.get_data(r"\efit_g_eqdsk:mw", tree_name="_efit_tree"),        # Number of horizontal grid points
            "nh": params.mds_conn.get_data(r"\efit_g_eqdsk:mh", tree_name="_efit_tree"),        # Number of vertical grid points
            "rdim": params.mds_conn.get_data(r"\efit_g_eqdsk:xdim", tree_name="_efit_tree"),    # Horizontal dimension of grid [m]
            "zdim": params.mds_conn.get_data(r"\efit_g_eqdsk:zdim", tree_name="_efit_tree"),    # Vertical dimension of grid [m]
        }
        rgrid = params.mds_conn.get_data(r"\efit_g_eqdsk:rgrid", tree_name="_efit_tree")  # R grid [m]
        zgrid = params.mds_conn.get_data(r"\efit_g_eqdsk:zgrid", tree_name="_efit_tree")  # Z grid [m]

        # Scalar parameters
        scalar_params = {
            # Plasma parameters
            "cpasma": params.mds_conn.get_data(r"\efit_g_eqdsk:cpasma", tree_name="_efit_tree"),  # Plasma current [A]
            "bcentr": params.mds_conn.get_data(r"\efit_g_eqdsk:bcentr", tree_name="_efit_tree"),  # Vacuum toroidal field [T]

            # Magnetic axis
            # Note that ssimag is -1 * simagx from the aeqdsk, apparently
            "rmaxis": params.mds_conn.get_data(r"\efit_g_eqdsk:rmaxis", tree_name="_efit_tree"), # R of magnetic axis [m]
            "zmaxis": params.mds_conn.get_data(r"\efit_g_eqdsk:zmaxis", tree_name="_efit_tree"), # Z of magnetic axis [m]
            "ssimag": params.mds_conn.get_data(r"\efit_g_eqdsk:ssimag", tree_name="_efit_tree"), # Flux at magnetic axis [Wb]

            # Plasma boundary
            "ssibry": params.mds_conn.get_data(r"\efit_g_eqdsk:ssibry", tree_name="_efit_tree"),   # Flux at plasma boundary [Wb]
        }

        # 1D profile parameters
        profile_params = {
            "fpol": params.mds_conn.get_data(r"\efit_g_eqdsk:fpol", tree_name="_efit_tree"),     # F = R*Bt [T*m]
            "pres": params.mds_conn.get_data(r"\efit_g_eqdsk:pres", tree_name="_efit_tree"),     # Pressure [Pa]
            "ffprim": params.mds_conn.get_data(r"\efit_g_eqdsk:ffprim", tree_name="_efit_tree"), # F*dF/dpsi
            "pprime": params.mds_conn.get_data(r"\efit_g_eqdsk:pprime", tree_name="_efit_tree"), # dp/dpsi
            "qpsi": params.mds_conn.get_data(r"\efit_g_eqdsk:qpsi", tree_name="_efit_tree"),     # Safety factor q
        }

        # 2D flux map
        psi_2d = params.mds_conn.get_data(r"\efit_g_eqdsk:psirz", tree_name="_efit_tree")  # 2D flux map [Wb]

        # Boundary and limiter data
        try:
            rbbbs = params.mds_conn.get_data(r"\efit_g_eqdsk:rbbbs", tree_name="_efit_tree")  # R of boundary points [m]
            zbbbs = params.mds_conn.get_data(r"\efit_g_eqdsk:zbbbs", tree_name="_efit_tree")  # Z of boundary points [m]
            nbbbs = params.mds_conn.get_data(r"\efit_g_eqdsk:nbbbs", tree_name="_efit_tree")  # Number of boundary points
        except:
            # Set empty boundary if not available
            rbbbs = np.array([])
            zbbbs = np.array([])
            nbbbs = 0

        # For some reason there is a difference in shape between xlim,ylim, and lim
        try:
            xlim = params.mds_conn.get_data(r"\efit_g_eqdsk:xlim", tree_name="_efit_tree")   # R of limiter points [m]
            ylim = params.mds_conn.get_data(r"\efit_g_eqdsk:ylim", tree_name="_efit_tree")   # Z of limiter points [m]
        except:
            # Set empty limiter if not available
            xlim = np.array([])
            ylim = np.array([])
        nlim = xlim.shape[0] if xlim.ndim > 0 else 0   # Number of limiter points

        # Create data variables for the dataset
        data_vars = {}

        # Add scalar parameters
        for param, values_raw in scalar_params.items():
            values_interp = interp1(efit_time, values_raw, params.times)
            data_vars[param] = xr.DataArray(
                values_interp,
                dims=["idx"],
                coords={"time": ("idx", params.times)},
                attrs={"description": f"GEQDSK parameter {param}"}
            )

        # Add 1D profile parameters.
        for param, values_raw in profile_params.items():
            if values_raw.ndim == 2:  # Time-dependent profiles
                # For some reason these are stored as (psi_idx, time) so transpose them
                values_interp = np.array([interp1(efit_time, values_raw[i, :], params.times) for i in range(values_raw.shape[0])]).T
                data_vars[param] = xr.DataArray(
                    values_interp,
                    dims=["idx", "psi_idx"],
                    coords={
                        "time": ("idx", params.times),
                        "psi_idx": np.arange(values_interp.shape[1])
                    },
                    attrs={"description": f"GEQDSK 1D profile {param}"}
                )
            else:  # Single profile
                data_vars[param] = xr.DataArray(
                    values_raw,
                    dims=["psi_idx"],
                    coords={"psi_idx": np.arange(len(values_raw))},
                    attrs={"description": f"GEQDSK 1D profile {param}"}
                )

        # Add 2D flux map
        if psi_2d.ndim == 3:  # Time-dependent 2D data
            psi_2d_interp = np.array([[interp1(efit_time, psi_2d[:, i, j], params.times) for j in range(psi_2d.shape[2])] for i in range(psi_2d.shape[1])]).T
            data_vars["psirz"] = xr.DataArray(
                psi_2d_interp,
                dims=["idx", "r_idx", "z_idx"],
                coords={
                    "time": ("idx", params.times),
                    "r_idx": ("r_idx", np.arange(psi_2d.shape[1])),
                    "z_idx": ("z_idx", np.arange(psi_2d.shape[2])),
                    "rgrid": ("r_idx", rgrid[0]),
                    "zgrid": ("z_idx", zgrid[0])
                },
                attrs={"description": "2D poloidal flux map", "units": "Wb"}
            )
        else:  # Single time slice
            data_vars["psirz"] = xr.DataArray(
                psi_2d,
                dims=["r_idx", "z_idx"],
                coords={
                    "r_idx": np.arange(psi_2d.shape[1]),
                    "z_idx": np.arange(psi_2d.shape[2]),
                    "rgrid": ("r_idx", rgrid[0]),
                    "zgrid": ("z_idx", zgrid[0])
                },
                attrs={"description": "2D poloidal flux map", "units": "Wb"}
            )

        # Add boundary data
        if len(rbbbs) > 0:
            if rbbbs.ndim == 2:  # Time-dependent boundary
                # Transpose to (time, boundary_idx)
                rbbbs_interp = np.array([interp1(efit_time, rbbbs[i, :], params.times) for i in range(rbbbs.shape[0])]).T
                data_vars["rbbbs"] = xr.DataArray(
                    rbbbs_interp,
                    dims=["idx", "boundary_idx"],
                    coords={
                        "time": ("idx", params.times),
                        "boundary_idx": np.arange(rbbbs_interp.shape[1])
                    },
                    attrs={"description": "R coordinates of plasma boundary", "units": "m"}
                )
                zbbbs_interp = np.array([interp1(efit_time, zbbbs[i, :], params.times) for i in range(zbbbs.shape[0])]).T
                data_vars["zbbbs"] = xr.DataArray(
                    zbbbs_interp,
                    dims=["idx", "boundary_idx"],
                    coords={
                        "time": ("idx", params.times),
                        "boundary_idx": np.arange(zbbbs_interp.shape[1])
                    },
                    attrs={"description": "Z coordinates of plasma boundary", "units": "m"}
                )
                nbbbs_interp = interp1(efit_time, nbbbs, params.times)
                data_vars["nbbbs"] = xr.DataArray(
                    nbbbs_interp,
                    dims=["idx"],
                    coords={"time": ("idx", params.times)},
                    attrs={"description": "Number of boundary points"}
                )
            else:  # Single boundary
                data_vars["rbbbs"] = xr.DataArray(
                    rbbbs,
                    dims=["boundary_idx"],
                    coords={"boundary_idx": np.arange(len(rbbbs))},
                    attrs={"description": "R coordinates of plasma boundary", "units": "m"}
                )
                data_vars["zbbbs"] = xr.DataArray(
                    zbbbs,
                    dims=["boundary_idx"],
                    coords={"boundary_idx": np.arange(len(zbbbs))},
                    attrs={"description": "Z coordinates of plasma boundary", "units": "m"}
                )
                data_vars["nbbbs"] = xr.DataArray(
                    nbbbs,
                    dims=[],
                    attrs={"description": "Number of boundary points"}
                )

        # Add limiter data and change the name to match GEQDSK convention
        if len(xlim) > 0:
            data_vars["rlim"] = xr.DataArray(
                xlim,
                dims=["limiter_idx"],
                coords={"limiter_idx": np.arange(len(xlim))},
                attrs={"description": "R coordinates of limiter", "units": "m"}
            )
            data_vars["zlim"] = xr.DataArray(
                ylim,
                dims=["limiter_idx"],
                coords={"limiter_idx": np.arange(len(ylim))},
                attrs={"description": "Z coordinates of limiter", "units": "m"}
            )
        data_vars["nlim"] = xr.DataArray(
            nlim,
            dims=[],
            attrs={"description": "Number of limiter points"}
        )

        data_vars.update(grid_params)

        # Create the dataset
        geqdsk_dataset = xr.Dataset(
            data_vars=data_vars,
            coords={
                "time": ("idx", params.times),
                "shot": ("idx", np.repeat(params.shot_id, len(params.times), axis=0)),
            },
            attrs={
                "description": "GEQDSK equilibrium parameters for C-Mod",
                "source": "EFIT reconstruction",
            }
        )

        return geqdsk_dataset


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
