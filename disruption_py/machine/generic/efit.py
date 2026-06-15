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
                **(
                    {
                        "rlim": ("limiter_idx", rlim),
                        "zlim": ("limiter_idx", zlim),
                    }
                    if rlim is not None and zlim is not None
                    else {}
                ),
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
                **(
                    {
                        "limiter_idx": np.arange(len(rlim)),
                    }
                    if rlim is not None
                    else {}
                ),
            },
            attrs={
                "cocos": cocos_input,
            },
        )

        return ds_geqdsk

    @staticmethod
    def efit_quality(
        shot_id: int,
        efit_time: np.ndarray,
        current: np.ndarray,
        wmhd: np.ndarray,
        chisq: np.ndarray,
        jflag: np.ndarray,
        num_iter: np.ndarray,
        max_iter: int,
        lflag: np.ndarray | None = None,
        terror: np.ndarray | None = None,
        cerror: np.ndarray | None = None,
    ):
        """
        Determine if specific timeslices are good and if the shot as a whole is good.

        Timeslice logic:
        - A timeslice is considered "perfect" (0) if it has no issues.
        - A timeslice is considered "good" (1) if it has minor issues, such as:
            - chisq > 10
            - jflag = 0 (indicating a non-converged solution, but it may still be usable)
            - lflag != 0 (indicating specific issues, only available on C-Mod)
            - terror > 0.005 (C-Mod convergence error threshold)
            - cerror > 0.01 (DIII-D convergence error threshold)
        - A timeslice is considered "bad" (2) if it has major issues, such as:
            - chisq > 50
            - num_iter >= max_iter (indicating failure to converge within the maximum iterations)
            - wmhd < 0

        Shot logic:
        - A shot is considered "perfect" (0) if all timeslices are perfect.
        - A shot is considered "good" (1) if it has no more than 20% bad timeslices and at least 80% good or perfect timeslices.
        - A shot is considered "bad" (2) if it has more than 20% bad timeslices or less than 80% good or perfect timeslices.

        Parameters:
            shot_id: int
                ID of the shot.
            efit_time: np.ndarray
                Time array corresponding to EFIT outputs.
            chisq: np.ndarray
                Chi-squared values from EFIT.
            jflag: np.ndarray
                JFLAG values from EFIT. 1 indicates good, 0 indicates bad
            num_iter: np.ndarray
                Number of iterations taken to converge
            max_iter: int
                Maximum number of iterations allowed to reach convergence
            lflag: np.ndarray | None
                2D array of specific error flags over time. Only available on C-Mod
            terror: np.ndarray | None
                Convergence error flag. Only available on C-Mod
            cerror: np.ndarray | None
                Convergence error flag. Only available on DIII-D

        Returns:
            ds_quality: An xarray Dataset containing the quality of each timeslice and the overall shot quality.
            - efit_quality_timeslice: 0 indicates a perfect timeslice (no issues), 1 indicates a good timeslice (minor issues), and 2 indicates a bad timeslice (major issues).
            - efit_quality_shot: 0 indicates a perfect shot (no issues), 1 indicates a good shot (minor issues), and 2 indicates a bad shot (major issues).
            - all non-nan input parameters
        """
        num_timesteps = len(efit_time)

        is_bad = np.asarray(num_iter) >= max_iter
        is_bad |= np.asarray(chisq) > 50
        is_bad |= np.asarray(wmhd) < 0

        is_minor = np.zeros(num_timesteps, dtype=bool)
        is_minor |= np.asarray(chisq) > 10
        is_minor |= np.asarray(jflag) == 0
        # If jflag is 1, we know that lflag is all cleared, no need to check it.
        if terror is not None:
            te = np.asarray(terror)
            is_minor |= np.isfinite(te) & (te > 0.005)
        if cerror is not None:
            ce = np.asarray(cerror)
            is_minor |= np.isfinite(ce) & (ce > 0.01)

        timeslice_quality = np.zeros(num_timesteps, dtype=int)
        timeslice_quality[is_minor] = 1
        timeslice_quality[is_bad] = 2

        bad_fraction = np.sum(is_bad) / num_timesteps if num_timesteps > 0 else 0.0
        if np.all(timeslice_quality == 0):
            shot_quality = 0
        elif bad_fraction <= 0.20:
            shot_quality = 1
        else:
            shot_quality = 2

        logger.debug(
            "Shot {shot}: {bad:.1%} bad timeslices, shot quality={q}",
            shot=shot_id,
            bad=bad_fraction,
            q=shot_quality,
        )

        def _include(arr: np.ndarray | None) -> np.ndarray | None:
            if arr is None:
                return None
            a = np.asarray(arr, dtype=float)
            return a if np.any(np.isfinite(a)) else None

        input_vars_1d = {
            "current": _include(current),
            "wmhd": _include(wmhd),
            "chisq": _include(chisq),
            "jflag": _include(jflag),
            "num_iter": _include(num_iter),
            "terror": _include(terror),
            "cerror": _include(cerror),
        }

        lflag_included = _include(lflag)
        lflag_var = {}
        if lflag_included is not None:
            lf = np.asarray(lflag_included)
            if lf.ndim == 2:
                # (n_flags, n_times) -> store as (n_times, n_flags)
                lflag_var = {"lflag": (("idx", "lflag_idx"), lf.T)}
            else:
                lflag_var = {"lflag": ("idx", lf)}

        return xr.Dataset(
            data_vars={
                "efit_quality_timeslice": ("idx", timeslice_quality),
                "efit_quality_shot": (
                    "idx",
                    np.full(num_timesteps, shot_quality, dtype=int),
                ),
                **{k: ("idx", v) for k, v in input_vars_1d.items() if v is not None},
                **lflag_var,
            },
            coords={
                "time": ("idx", efit_time),
                "shot": ("idx", np.repeat(shot_id, num_timesteps)),
                **(
                    {"max_iter": float(max_iter)}
                    if np.isfinite(float(max_iter))
                    else {}
                ),
            },
        )
