#!/usr/bin/env python3

"""
Module for retrieving and processing EFIT parameters for DIII-D.
"""

import numpy as np
import xarray as xr

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
from disruption_py.machine.tokamak import Tokamak


class D3DEfitMethods:
    """
    A class for retrieving and processing EFIT parameters from DIII-D.
    """

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

    efit_derivs = {"dbetap_dt": "beta_p", "dli_dt": "li", "dwmhd_dt": "wmhd"}
    rt_efit_cols = {
        "beta_p_rt": r"\efit_a_eqdsk:betap",
        "li_rt": r"\efit_a_eqdsk:li",
        "q95_rt": r"\efit_a_eqdsk:q95",
        "wmhd_rt": r"\efit_a_eqdsk:wmhd",
        "chisq_rt": r"\efit_a_eqdsk:chisq",
    }

    @staticmethod
    @physics_method(
        columns=[*efit_cols.keys(), *efit_derivs.keys()],
        tokamak=Tokamak.D3D,
    )
    def get_efit_parameters(params: PhysicsMethodParams):
        """
        Retrieve EFIT parameters and their time derivatives.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information.

        Returns
        -------
        dict
            A dictionary containing the EFIT parameters and their derivatives.
        """
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

        for param in efit_data:
            efit_data[param][invalid_indices] = np.nan
        for deriv_param, param in D3DEfitMethods.efit_derivs.items():
            efit_data[deriv_param] = np.gradient(efit_data[param], efit_time)
        if not np.array_equal(params.times, efit_time):
            for param in efit_data:
                efit_data[param] = interp1(efit_time, efit_data[param], params.times)
        return efit_data

    @staticmethod
    @physics_method(
        columns=[*rt_efit_cols.keys()],
        tokamak=Tokamak.D3D,
    )
    def get_rt_efit_parameters(params: PhysicsMethodParams):
        """
        Retrieve real-time EFIT parameters.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information.

        Returns
        -------
        dict
            A dictionary containing the real-time EFIT parameters.
        """
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
        invalid_indices = np.where(efit_data["chisq_rt"] > 50)

        for param in efit_data:
            efit_data[param][invalid_indices] = np.nan
        if not np.array_equal(params.times, efit_time):
            for param in efit_data:
                efit_data[param] = interp1(efit_time, efit_data[param], params.times)
        return efit_data

    @staticmethod
    @physics_method(
        tokamak=Tokamak.D3D,
    )
    def get_geqdsk_parameters(params: PhysicsMethodParams):
        """
        Retrieve parameters required to create an GeqdskDataDict https://freeqdsk.readthedocs.io/en/stable/geqdsk.html
        Also see https://fusion.gat.com/theory/Efitgeqdsk

        This requires the following:
        
        1. shot
        2. grid parameters
        - Either (nx, ny, rdim, zdim, rcenter, rleft, zmid) or (rgrid, zgrid)
        3. scalar parameters
        - rmagx, zmagx, simagx, sibdry, bccenter, cpasma, 
        4. 1D profile parameters
        - fpol, pres, ffprime, pprime, qpsi
        5. 2D parameters
        - psi or f or psirz
        6. boundary and limiter data (if available)
        - nbdry, nlim, rbdry, zbdry, rlim, zlim

        Some signals are missing on DIII-D:
        - rcentr

        Also, since it doesn't really make sense to get interpolated 'nbbbs', we're not including it in the dataset.
        Later you can look at the nonzero values in 'rbbbs' or 'zbbbs' to determine how many boundary points there are at each time.

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
        efit_time = params.mds_conn.get_data(
            r"\efit_a_eqdsk:atime", tree_name="_efit_tree"
        ) / 1.0e3  # [ms] -> [s]

        # Grid parameters
        # For DIII-D, rcentr and rleft are missing but can be inferred from rgrid
        grid_params = {
            "nw": params.mds_conn.get_data(r"\efit_g_eqdsk:mw", tree_name="_efit_tree"),        # Number of horizontal grid points
            "nh": params.mds_conn.get_data(r"\efit_g_eqdsk:mh", tree_name="_efit_tree"),        # Number of vertical grid points
            "rdim": params.mds_conn.get_data(r"\efit_g_eqdsk:xdim", tree_name="_efit_tree"),    # Horizontal dimension of grid [m]
            "zdim": params.mds_conn.get_data(r"\efit_g_eqdsk:zdim", tree_name="_efit_tree"),    # Vertical dimension of grid [m]
            "zmid": params.mds_conn.get_data(r"\efit_g_eqdsk:zmid", tree_name="_efit_tree"),     # Z at center of grid [m]
        }
        rgrid = params.mds_conn.get_data(r"\efit_g_eqdsk:r", tree_name="_efit_tree"),  # R grid [m]
        zgrid = params.mds_conn.get_data(r"\efit_g_eqdsk:z", tree_name="_efit_tree"),  # Z grid [m]
        
        # Scalar parameters
        scalar_params = {
            # Plasma parameters
            "cpasma": params.mds_conn.get_data(r"\efit_g_eqdsk:cpasma", tree_name="_efit_tree"),  # Plasma current [A]
            "bcentr": params.mds_conn.get_data(r"\efit_g_eqdsk:bcentr", tree_name="_efit_tree"),  # Vacuum toroidal field [T]

            # Magnetic axis
            "rmaxis": params.mds_conn.get_data(r"\efit_g_eqdsk:rmaxis", tree_name="_efit_tree"), # R of magnetic axis [m]
            "zmaxis": params.mds_conn.get_data(r"\efit_g_eqdsk:zmaxis", tree_name="_efit_tree"), # Z of magnetic axis [m]
            "simag": params.mds_conn.get_data(r"\efit_g_eqdsk:ssimag", tree_name="_efit_tree"),   # Flux at magnetic axis [Wb]
            
            # Plasma boundary
            "sibry": params.mds_conn.get_data(r"\efit_g_eqdsk:ssibry", tree_name="_efit_tree"),   # Flux at plasma boundary [Wb]
            "xdim": params.mds_conn.get_data(r"\efit_g_eqdsk:xdim", tree_name="_efit_tree"),     # X-point R location [m]
            "zdim": params.mds_conn.get_data(r"\efit_g_eqdsk:zdim", tree_name="_efit_tree"),     # X-point Z location [m]
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
            
        try:
            rlim = params.mds_conn.get_data(r"\efit_g_eqdsk:xlim", tree_name="_efit_tree")   # R of limiter points [m]
            zlim = params.mds_conn.get_data(r"\efit_g_eqdsk:ylim", tree_name="_efit_tree")   # Z of limiter points [m]
            nlim = params.mds_conn.get_data(r"\efit_g_eqdsk:nlim", tree_name="_efit_tree")   # Number of limiter points
        except:
            # Set empty limiter if not available
            rlim = np.array([])
            zlim = np.array([])
            nlim = 0
        
        # Create time-dependent data arrays
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
                # These are stored as (time, psi_idx)
                values_interp = np.array([interp1(efit_time, values_raw[:, i], params.times) for i in range(values_raw.shape[1])]).T
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
                # These are stored as (time, boundary_idx)
                rbbbs_interp = np.array([interp1(efit_time, rbbbs[:, i], params.times) for i in range(rbbbs.shape[1])]).T
                data_vars["rbbbs"] = xr.DataArray(
                    rbbbs_interp,
                    dims=["idx", "boundary_idx"],
                    coords={
                        "time": ("idx", params.times),
                        "boundary_idx": np.arange(rbbbs_interp.shape[1])
                    },
                    attrs={"description": "R coordinates of plasma boundary", "units": "m"}
                )
                zbbbs_interp = np.array([interp1(efit_time, zbbbs[:, i], params.times) for i in range(zbbbs.shape[1])]).T
                data_vars["zbbbs"] = xr.DataArray(
                    zbbbs_interp,
                    dims=["idx", "boundary_idx"],
                    coords={
                        "time": ("idx", params.times),
                        "boundary_idx": np.arange(zbbbs_interp.shape[1])
                    },
                    attrs={"description": "Z coordinates of plasma boundary", "units": "m"}
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
        
        # Add limiter data and change the name to match GEQDSK convention
        if len(rlim) > 0:
            data_vars["rlim"] = xr.DataArray(
                rlim,
                dims=["limiter_idx"],
                coords={"limiter_idx": np.arange(len(rlim))},
                attrs={"description": "R coordinates of limiter", "units": "m"}
            )
            data_vars["zlim"] = xr.DataArray(
                zlim,
                dims=["limiter_idx"],
                coords={"limiter_idx": np.arange(len(zlim))},
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
