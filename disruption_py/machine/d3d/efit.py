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
            r"\efit_g_eqdsk:atime", tree_name="_efit_tree"
        ) / 1.0e3  # [ms] -> [s]
        
        # Scalar parameters
        scalar_params = {
            # Grid dimensions
            "nw": params.mds_conn.get_data(r"\efit_g_eqdsk:mw", tree_name="_efit_tree"),  # Number of horizontal grid points
            "nh": params.mds_conn.get_data(r"\efit_g_eqdsk:mh", tree_name="_efit_tree"),  # Number of vertical grid points
            
            # Plasma parameters
            "cpasma": params.mds_conn.get_data(r"\efit_g_eqdsk:cpasma", tree_name="_efit_tree"),  # Plasma current [A]
            "bcentr": params.mds_conn.get_data(r"\efit_g_eqdsk:bcentr", tree_name="_efit_tree"),  # Vacuum toroidal field [T]
            "rcentr": params.mds_conn.get_data(r"\efit_g_eqdsk:rcentr", tree_name="_efit_tree"),  # Reference major radius [m]
            
            # Grid boundaries
            "rdim": params.mds_conn.get_data(r"\efit_g_eqdsk:rdim", tree_name="_efit_tree"),    # Horizontal extent of grid [m]
            "zdim": params.mds_conn.get_data(r"\efit_g_eqdsk:zdim", tree_name="_efit_tree"),    # Vertical extent of grid [m]
            "rcentr": params.mds_conn.get_data(r"\efit_g_eqdsk:rcentr", tree_name="_efit_tree"), # R at center of grid [m]
            "rleft": params.mds_conn.get_data(r"\efit_g_eqdsk:rleft", tree_name="_efit_tree"),   # R at left edge of grid [m]
            "zmid": params.mds_conn.get_data(r"\efit_g_eqdsk:zmid", tree_name="_efit_tree"),     # Z at center of grid [m]
            
            # Magnetic axis
            "rmaxis": params.mds_conn.get_data(r"\efit_g_eqdsk:rmaxis", tree_name="_efit_tree"), # R of magnetic axis [m]
            "zmaxis": params.mds_conn.get_data(r"\efit_g_eqdsk:zmaxis", tree_name="_efit_tree"), # Z of magnetic axis [m]
            "simag": params.mds_conn.get_data(r"\efit_g_eqdsk:simag", tree_name="_efit_tree"),   # Flux at magnetic axis [Wb]
            
            # Plasma boundary
            "sibry": params.mds_conn.get_data(r"\efit_g_eqdsk:sibry", tree_name="_efit_tree"),   # Flux at plasma boundary [Wb]
            "xdim": params.mds_conn.get_data(r"\efit_g_eqdsk:xdim", tree_name="_efit_tree"),     # X-point R location [m]
            "zdim": params.mds_conn.get_data(r"\efit_g_eqdsk:ydim", tree_name="_efit_tree"),     # X-point Z location [m]
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
            rlim = params.mds_conn.get_data(r"\efit_g_eqdsk:rlim", tree_name="_efit_tree")   # R of limiter points [m]
            zlim = params.mds_conn.get_data(r"\efit_g_eqdsk:zlim", tree_name="_efit_tree")   # Z of limiter points [m]
            nlim = params.mds_conn.get_data(r"\efit_g_eqdsk:nlim", tree_name="_efit_tree")   # Number of limiter points
        except:
            # Set empty limiter if not available
            rlim = np.array([])
            zlim = np.array([])
            nlim = 0
        
        # Create coordinate arrays
        nw = scalar_params["nw"][0] if hasattr(scalar_params["nw"], '__len__') else scalar_params["nw"]
        nh = scalar_params["nh"][0] if hasattr(scalar_params["nh"], '__len__') else scalar_params["nh"]
        
        # Create time-dependent data arrays
        data_vars = {}
        
        # Add scalar parameters
        for param, values in scalar_params.items():
            data_vars[param] = xr.DataArray(
                values,
                dims=["time"],
                coords={"time": efit_time},
                attrs={"description": f"GEQDSK parameter {param}"}
            )
        
        # Add 1D profile parameters
        for param, values in profile_params.items():
            if values.ndim == 2:  # Time-dependent profiles
                data_vars[param] = xr.DataArray(
                    values,
                    dims=["time", "psi_index"],
                    coords={
                        "time": efit_time,
                        "psi_index": np.arange(values.shape[1])
                    },
                    attrs={"description": f"GEQDSK 1D profile {param}"}
                )
            else:  # Single profile
                data_vars[param] = xr.DataArray(
                    values,
                    dims=["psi_index"],
                    coords={"psi_index": np.arange(len(values))},
                    attrs={"description": f"GEQDSK 1D profile {param}"}
                )
        
        # Add 2D flux map
        if psi_2d.ndim == 3:  # Time-dependent 2D data
            data_vars["psirz"] = xr.DataArray(
                psi_2d,
                dims=["time", "r_index", "z_index"],
                coords={
                    "time": efit_time,
                    "r_index": np.arange(psi_2d.shape[1]),
                    "z_index": np.arange(psi_2d.shape[2])
                },
                attrs={"description": "2D poloidal flux map", "units": "Wb"}
            )
        else:  # Single time slice
            data_vars["psirz"] = xr.DataArray(
                psi_2d,
                dims=["r_index", "z_index"],
                coords={
                    "r_index": np.arange(psi_2d.shape[0]),
                    "z_index": np.arange(psi_2d.shape[1])
                },
                attrs={"description": "2D poloidal flux map", "units": "Wb"}
            )
        
        # Add boundary data
        if len(rbbbs) > 0:
            if rbbbs.ndim == 2:  # Time-dependent boundary
                data_vars["rbbbs"] = xr.DataArray(
                    rbbbs,
                    dims=["time", "boundary_index"],
                    coords={
                        "time": efit_time,
                        "boundary_index": np.arange(rbbbs.shape[1])
                    },
                    attrs={"description": "R coordinates of plasma boundary", "units": "m"}
                )
                data_vars["zbbbs"] = xr.DataArray(
                    zbbbs,
                    dims=["time", "boundary_index"],
                    coords={
                        "time": efit_time,
                        "boundary_index": np.arange(zbbbs.shape[1])
                    },
                    attrs={"description": "Z coordinates of plasma boundary", "units": "m"}
                )
            else:  # Single boundary
                data_vars["rbbbs"] = xr.DataArray(
                    rbbbs,
                    dims=["boundary_index"],
                    coords={"boundary_index": np.arange(len(rbbbs))},
                    attrs={"description": "R coordinates of plasma boundary", "units": "m"}
                )
                data_vars["zbbbs"] = xr.DataArray(
                    zbbbs,
                    dims=["boundary_index"],
                    coords={"boundary_index": np.arange(len(zbbbs))},
                    attrs={"description": "Z coordinates of plasma boundary", "units": "m"}
                )
        
        # Add limiter data
        if len(rlim) > 0:
            data_vars["rlim"] = xr.DataArray(
                rlim,
                dims=["limiter_index"],
                coords={"limiter_index": np.arange(len(rlim))},
                attrs={"description": "R coordinates of limiter", "units": "m"}
            )
            data_vars["zlim"] = xr.DataArray(
                zlim,
                dims=["limiter_index"],
                coords={"limiter_index": np.arange(len(zlim))},
                attrs={"description": "Z coordinates of limiter", "units": "m"}
            )
        
        # Add number of boundary and limiter points
        data_vars["nbbbs"] = xr.DataArray(
            nbbbs,
            dims=["time"] if hasattr(nbbbs, '__len__') and len(nbbbs) > 1 else [],
            coords={"time": efit_time} if hasattr(nbbbs, '__len__') and len(nbbbs) > 1 else {},
            attrs={"description": "Number of plasma boundary points"}
        )
        data_vars["nlim"] = xr.DataArray(
            nlim,
            dims=[],
            attrs={"description": "Number of limiter points"}
        )
        
        # Create the dataset
        geqdsk_dataset = xr.Dataset(
            data_vars,
            attrs={
                "shot": params.shot_id,
                "description": "GEQDSK equilibrium parameters for DIII-D",
                "source": "EFIT reconstruction",
            }
        )
        
        # Interpolate to requested timebase if different
        if not np.array_equal(params.times, efit_time):
            # Only interpolate time-dependent variables
            time_dependent_vars = [var for var in geqdsk_dataset.data_vars 
                                 if "time" in geqdsk_dataset[var].dims]
            for var in time_dependent_vars:
                if geqdsk_dataset[var].ndim == 1:  # 1D time series
                    geqdsk_dataset[var] = xr.DataArray(
                        interp1(efit_time, geqdsk_dataset[var].values, params.times),
                        dims=["time"],
                        coords={"time": params.times},
                        attrs=geqdsk_dataset[var].attrs
                    )
                elif geqdsk_dataset[var].ndim == 2 and "time" in geqdsk_dataset[var].dims:  # 2D with time
                    interpolated_data = np.array([
                        interp1(efit_time, geqdsk_dataset[var].values[:, i], params.times)
                        for i in range(geqdsk_dataset[var].shape[1])
                    ]).T
                    other_dim = [dim for dim in geqdsk_dataset[var].dims if dim != "time"][0]
                    geqdsk_dataset[var] = xr.DataArray(
                        interpolated_data,
                        dims=["time", other_dim],
                        coords={
                            "time": params.times,
                            other_dim: geqdsk_dataset[var].coords[other_dim]
                        },
                        attrs=geqdsk_dataset[var].attrs
                    )
        
        return geqdsk_dataset



