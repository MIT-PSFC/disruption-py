#!/usr/bin/env python3

"""
Module for generic EFIT methods.
"""

import numpy as np
import xarray as xr
from loguru import logger

class GenericEFITMethods:
    """
    Class to hold generic EFIT methods.
    """

    @staticmethod
    def make_geqdsk_dataset(
        shot_id: int,
        times: np.ndarray,
        # Grid info
        r_grid: np.ndarray,
        z_grid: np.ndarray,
        # 0D signals
        rmagx: float,
        zmagx: float,
        simagx: float,
        sibdry: float,
        bcentr: float,
        current: float,
        # 1D signals
        fpol: np.ndarray,
        pres: np.ndarray,
        ffprime: np.ndarray,
        pprime: np.ndarray,
        qpsi: np.ndarray,
        # 2D signals
        psirz: np.ndarray,
        # Boundary
        rbdry: np.ndarray,
        zbdry: np.ndarray,
        # COCOS info
        cocos_input: int,
        # Limiter (optional)
        rlim: np.ndarray | None = None,
        zlim: np.ndarray | None = None,
    ) -> xr.Dataset:
        """
        Create an Xarray dataset containing all signals necessary to recreate a GEQDSK file.
        Dataset has signals that have FREEQDSK names and are in COCOS 1

        Parameters
        ----------
        r_grid : np.ndarray
            Radial grid points.
        z_grid : np.ndarray
            Vertical grid points.
        rmagx : float
            R coordinate of magnetic axis.
        zmagx : float
            Z coordinate of magnetic axis.
        simagx : float
            Poloidal flux at magnetic axis.
        sibdry : float
            Poloidal flux at plasma boundary.
        bcentr : float
            Toroidal magnetic field at the center.
        current : float
            Plasma current.
        qpsi : np.ndarray
            Safety factor profile as a function of poloidal flux.
        fpol : np.ndarray
            Poloidal current function profile as a function of poloidal flux.
        pres : np.ndarray
            Pressure profile as a function of poloidal flux.
        ffprime : np.ndarray
            Derivative of fpol with respect to poloidal flux.
        pprime : np.ndarray
            Derivative of pressure with respect to poloidal flux.
        psirz : np.ndarray
            2D array of poloidal flux values on the (R,Z) grid.
        rbdry : np.ndarray
            R coordinates of the plasma boundary points.
        zbdry : np.ndarray
            Z coordinates of the plasma boundary points.
        cocos_input : int
            Input COCOS for the dataset.
        rlim : np.ndarray, optional
            R coordinates of the limiter points.
        zlim : np.ndarray, optional
            Z coordinates of the limiter points.
        """

        # Compute grid center and dimensions
        rcentr = r_grid[len(r_grid) // 2]
        rleft = r_grid[0]
        rdim = r_grid[-1] - r_grid[0]
        zmid = z_grid[len(z_grid) // 2]
        zdim = z_grid[-1] - z_grid[0]

        ds_geqdsk = xr.Dataset(
            data_vars={
                # 0D signals (names match make_eq_field expectations)
                "rmagx": ("idx", rmagx),
                "zmagx": ("idx", zmagx),
                "simagx": ("idx", simagx),
                "sibdry": ("idx", sibdry),
                "bcentr": ("idx", bcentr),
                "current": ("idx", current),
                # 1D signals
                "qpsi": (("idx", "psi_idx"), qpsi),
                "fpol": (("idx", "psi_idx"), fpol),
                "pres": (("idx", "psi_idx"), pres),
                "ffprime": (("idx", "psi_idx"), ffprime),
                "pprime": (("idx", "psi_idx"), pprime),
                # 2D signals
                "psirz": (("idx", "r_grid", "z_grid"), psirz),
                # Boundary
                "rbdry": (("idx", "boundary_idx"), rbdry),
                "zbdry": (("idx", "boundary_idx"), zbdry),
                # Limiter (static - not time-dependent)
                **({
                    "rlim": ("limiter_idx", rlim),
                    "zlim": ("limiter_idx", zlim),
                } if rlim is not None and zlim is not None else {}),
            },
            coords={
                "time": ("idx", times),
                "shot": ("idx", np.repeat(shot_id, len(times), axis=0)),
                # Grid info
                "rcentr": rcentr,
                "rleft": rleft,
                "zmid": zmid,
                "rdim": rdim,
                "zdim": zdim,
                "r_grid": r_grid,
                "z_grid": z_grid,
                # profile, LCFS, and limiter indices
                "psi_idx": np.arange(qpsi.shape[1]),
                "boundary_idx": np.arange(rbdry.shape[1]),
                **({
                    "limiter_idx": np.arange(len(rlim)),
                } if rlim is not None else {}),
            },
            attrs={
                "cocos": cocos_input,
            }
        )

        return ds_geqdsk
    
