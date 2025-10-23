#!/usr/bin/env python3

"""
Module for retrieving and processing EFIT parameters for CMOD.
"""

import numpy as np
from MDSplus import mdsExceptions

import xarray as xr

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
from disruption_py.machine.tokamak import Tokamak


class CmodEfitMethods:
    """
    Class for retrieving and processing EFIT parameters for CMOD.

    Attributes
    ----------
    efit_parameter_cols : dict
        A dictionary mapping parameter names to their corresponding aeqdsk EFIT data paths.
    efit_equilibrium_cols : dict
        A dictionary mapping equilibrium names to their corresponding geqdsk EFIT data paths.
    efit_derivs : dict
        A dictionary mapping derivative parameter names to their corresponding base parameters.
    """

    efit_parameter_cols = {
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

    # This should contain everything which is necessary to recreate a geqdsk
    # https://freeqdsk.readthedocs.io/en/stable/geqdsk.html#freeqdsk.geqdsk.write
    efit_equilibrium_cols = {
        "rgrid": {"path": r"\efit_geqdsk:rgrid", "val": "grid"},
        "zgrid": {"path": r"\efit_geqdsk:zgrid", "val": "grid"},
        "psirz": {"path": r"\efit_geqdsk:psirz", "val": "2D"},
        "psin": {"path": r"\efit_geqdsk:psin", "val": "grid"},
        "rhovn": {"path": r"\efit_geqdsk:rhovn", "val": "1D"},
        "ssibry": {"path": r"\efit_geqdsk:ssibry", "val": "scalar"},
        "ssimag": {"path": r"\efit_geqdsk:ssimag", "val": "scalar"},
    }

    @staticmethod
    @physics_method(
        columns=[
            *efit_parameter_cols.keys(),
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

        # Get data from each of the columns in efit_parameter_cols one at a time
        for param, path in CmodEfitMethods.efit_parameter_cols.items():
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
    def _get_efit_equilibrium(params: PhysicsMethodParams, timebase):
        efit_time = params.mds_conn.get_data(
            r"\efit_aeqdsk:time", tree_name="_efit_tree"
        )
        efit_equilibrium = {}
        for param, info in CmodEfitMethods.efit_equilibrium_cols.items():
            path = info["path"]
            val = info["val"]
            try:
                param_data = params.mds_conn.get_data(
                    path=path,
                    tree_name="_efit_tree",
                )
                if val == "grid":
                    # Just the grid, doesn't change with time
                    efit_equilibrium[param] = param_data
                elif val in ["1D", "scalar"]:
                    # 1D data, interpolate in time
                    interp_data = interp1(efit_time, param_data, timebase)
                    efit_equilibrium[param] = interp_data
                elif val in ["1D", "2D"]:
                    # TODO(ZanderKeith): There's some tomfoolery here.
                    # For something like psirz, the data has shape (time, grid, grid), but I'm unsure if it goes (time, r, z) or (time, z, r).
                    # For now I'm assuming it starts as (time, z, r) since that's how it would be transposed from how it's built in mdsplus (r, z, time)
                    # That means returning from this function it will be in (r, z, time)
                    interp_data = interp1(efit_time, param_data.T, timebase)
                    efit_equilibrium[param] = interp_data
            except mdsExceptions.MdsException as e:
                params.logger.warning(repr(e))
                params.logger.opt(exception=True).debug(e)
                efit_equilibrium[param] = np.full(len(timebase), np.nan)

        return efit_equilibrium
    
    @staticmethod
    @physics_method(
        columns=[
            *efit_equilibrium_cols.keys(),
        ],
        tokamak=Tokamak.CMOD,
    )
    def get_efit_equilibrium(params: PhysicsMethodParams):
        """
        Retrieve EFIT equilibrium data for CMOD.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDS connection and shot information.

        Returns
        -------
        dict
            A dictionary containing the retrieved EFIT equilibrium data.
        """

        return CmodEfitMethods._get_efit_equilibrium(params, params.times)
    
    @staticmethod
    @physics_method(
        tokamak=Tokamak.CMOD,
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

        Sometimes the C-Mod EFIT uses 'x' and 'y' instead of 'r' and 'z' for parameter names
        This function will convert them to the standard GEQDSK naming convention.

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
        
        # Scalar parameters
        scalar_params = {
            # Grid dimensions
            "nw": params.mds_conn.get_data(r"\efit_g_eqdsk:mw", tree_name="_efit_tree"),  # Number of horizontal grid points
            "nh": params.mds_conn.get_data(r"\efit_g_eqdsk:mh", tree_name="_efit_tree"),  # Number of vertical grid points
            
            # Plasma parameters
            "cpasma": params.mds_conn.get_data(r"\efit_g_eqdsk:cpasma", tree_name="_efit_tree"),  # Plasma current [A]
            "bcentr": params.mds_conn.get_data(r"\efit_g_eqdsk:bcentr", tree_name="_efit_tree"),  # Vacuum toroidal field [T]
            
            # Magnetic axis
            "rmaxis": params.mds_conn.get_data(r"\efit_g_eqdsk:rmaxis", tree_name="_efit_tree"), # R of magnetic axis [m]
            "zmaxis": params.mds_conn.get_data(r"\efit_g_eqdsk:zmaxis", tree_name="_efit_tree"), # Z of magnetic axis [m]
            
            # Plasma boundary
            # Note that ssimag is -1 * simagx from the aeqdsk, apparently
            "ssimag": params.mds_conn.get_data(r"\efit_g_eqdsk:ssimag", tree_name="_efit_tree"),   # Flux at magnetic axis [Wb]
            "ssibry": params.mds_conn.get_data(r"\efit_g_eqdsk:ssibry", tree_name="_efit_tree"),   # Flux at plasma boundary [Wb]
        }

        # Grid boundaries (For C-Mod a lot are missing, but can be inferred from rgrid/zgrid)
        grid_params = {
            "rdim": params.mds_conn.get_data(r"\efit_g_eqdsk:xdim", tree_name="_efit_tree"),     # Horizontal dimension of grid [m]
            "zdim": params.mds_conn.get_data(r"\efit_g_eqdsk:zdim", tree_name="_efit_tree"),     # Vertical dimension of grid [m]
            "rgrid": params.mds_conn.get_data(r"\efit_g_eqdsk:rgrid", tree_name="_efit_tree"),  # R grid [m]
            "zgrid": params.mds_conn.get_data(r"\efit_g_eqdsk:zgrid", tree_name="_efit_tree"),  # Z grid [m]
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
                    dims=["psi_index", "time"],
                    coords={
                        "time": efit_time,
                        "psi_index": np.arange(values.shape[0])
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
                    dims=["boundary_index", "time"],
                    coords={
                        "time": efit_time,
                        "boundary_index": np.arange(rbbbs.shape[0])
                    },
                    attrs={"description": "R coordinates of plasma boundary", "units": "m"}
                )
                data_vars["zbbbs"] = xr.DataArray(
                    zbbbs,
                    dims=["boundary_index", "time"],
                    coords={
                        "time": efit_time,
                        "boundary_index": np.arange(zbbbs.shape[0])
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
        
        # Add limiter data and change the name to match GEQDSK convention
        if len(xlim) > 0:
            data_vars["rlim"] = xr.DataArray(
                xlim,
                dims=["limiter_index"],
                coords={"limiter_index": np.arange(len(xlim))},
                attrs={"description": "R coordinates of limiter", "units": "m"}
            )
            data_vars["zlim"] = xr.DataArray(
                ylim,
                dims=["limiter_index"],
                coords={"limiter_index": np.arange(len(ylim))},
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

        data_vars.update(
            {"rdim": grid_params["rdim"],
             "zdim": grid_params["zdim"],}
        )
        data_vars["rgrid"] = xr.DataArray(
            grid_params["rgrid"],
            dims=["r_index"],
            coords={"r_index": np.arange(len(grid_params["rgrid"]))},
            attrs={"description": "R grid points", "units": "m"}
        )
        data_vars["zgrid"] = xr.DataArray(
            grid_params["zgrid"],
            dims=["z_index"],
            coords={"z_index": np.arange(len(grid_params["zgrid"]))},
            attrs={"description": "Z grid points", "units": "m"}
        )
        
        # Create the dataset
        geqdsk_dataset = xr.Dataset(
            data_vars,
            attrs={
                "description": "GEQDSK equilibrium parameters for C-Mod",
                "source": "EFIT reconstruction",
            }
        )
        
        # Interpolate time-dependent variables to params.times, and replace the 

        # Only interpolate time-dependent variables
        time_dependent_vars = [var for var in geqdsk_dataset.data_vars 
                                if "time" in geqdsk_dataset[var].dims]
        time_independent_vars = [var for var in geqdsk_dataset.data_vars 
                                  if "time" not in geqdsk_dataset[var].dims]
        for var in time_dependent_vars:
            if geqdsk_dataset[var].ndim == 1:  # scalar time series
                geqdsk_dataset[var] = xr.DataArray(
                    interp1(efit_time, geqdsk_dataset[var].values, params.times),
                    dims=["time"],
                    coords={"time": params.times},
                    attrs=geqdsk_dataset[var].attrs
                )
            elif geqdsk_dataset[var].ndim == 2 and "time" in geqdsk_dataset[var].dims:  # 1D with time
                interpolated_data = np.array([interp1(efit_time, geqdsk_dataset[var].values[i, :], params.times)for i in range(geqdsk_dataset[var].shape[0])])
                other_dim = [dim for dim in geqdsk_dataset[var].dims if dim != "time"][0]
                geqdsk_dataset[var] = xr.DataArray(
                    interpolated_data,
                    dims=[other_dim, "time"],
                    coords={
                        "time": params.times,
                        other_dim: geqdsk_dataset[var].coords[other_dim]
                    },
                    attrs=geqdsk_dataset[var].attrs
                )
            elif geqdsk_dataset[var].ndim == 3 and "time" in geqdsk_dataset[var].dims:  # 2D with time
                interpolated_data = np.array([[interp1(efit_time, geqdsk_dataset[var].values[:, i, j], params.times) for j in range(geqdsk_dataset[var].shape[2])] for i in range(geqdsk_dataset[var].shape[1])])
                other_dims = [dim for dim in geqdsk_dataset[var].dims if dim != "time"]
                geqdsk_dataset[var] = xr.DataArray(
                    interpolated_data,
                    dims=[*other_dims, "time"],
                    coords={
                        "time": params.times,
                        other_dims[0]: geqdsk_dataset[var].coords[other_dims[0]],
                        other_dims[1]: geqdsk_dataset[var].coords[other_dims[1]],
                    },
                    attrs=geqdsk_dataset[var].attrs
                )

            

        # Modify dataset so "time" coordinate is on dimension "idx"
        geqdsk_dataset = geqdsk_dataset.rename_dims({"time": "idx"})
        geqdsk_dataset = geqdsk_dataset.assign_coords({"time": ("idx", params.times), "shot": ("idx", params.shot_id)})
        
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
