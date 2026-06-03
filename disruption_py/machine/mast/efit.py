#!/usr/bin/env python3

"""
Module for retrieving and processing EFIT parameters for MAST.
"""

import numpy as np
import xarray as xr

from disruption_py.config import config as get_config
from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
from disruption_py.machine.generic.efit import GenericEFITMethods
from disruption_py.machine.mast.util import MastUtilMethods
from disruption_py.machine.tokamak import Tokamak


class MastEfitMethods:
    """
    Class for retrieving and processing EFIT parameters for MAST.
    """

    efit_properties = {
        "a_minor": "minor_radius",
        "beta_n": "beta_tor_normal",
        "beta_p": "beta_pol",
        "beta_t": "beta_tor",
        "bphi_rmag": "bphi_rmag",
        "bvac_rmag": "bvac_rmag",
        "kappa": "elongation",
        "li": "li",
        "q95": "q95",
        "rmagx": "magnetic_axis_r",
        "rmagz": "magnetic_axis_z",
        "tribot": "triangularity_lower",
        "tritop": "triangularity_upper",
        "v_loop_dynamic": "vloop_dynamic",
        "v_loop_static": "vloop_static",
        "wmhd": "wmhd",
    }

    # Output columns for get_geqdsk_parameters.
    # rmagx excluded: already returned by get_efit_parameters; including it causes
    # a merge conflict since level1 and level2 sources give slightly different values.
    geqdsk_output_cols = [
        "zmagx", "simagx", "sibdry", "bcentr", "current",
        "fpol", "pres", "ffprime", "pprime", "qpsi",
        "psirz", "rbdry", "zbdry",
    ]

    @staticmethod
    @physics_method(columns=list(efit_properties.keys()), tokamak=Tokamak.MAST)
    def get_efit_parameters(params: PhysicsMethodParams):
        """
        Retrieve EFIT parameters for MAST.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing theconnection and shot information.

        Returns
        -------
        dict
            A dictionary containing the retrieved EFIT parameters.
        """
        eq_time = params.get_data("equilibrium/time")
        times = params.times

        outputs = {}
        for key, prop in MastEfitMethods.efit_properties.items():
            signal = params.get_data(f"equilibrium/{prop}")
            item = MastUtilMethods.interpolate_1d(eq_time, signal, times)
            outputs[key] = item

        return outputs

    @staticmethod
    @physics_method(
        columns=geqdsk_output_cols,
        tokamak=Tokamak.MAST,
    )
    def get_geqdsk_parameters(params: PhysicsMethodParams):
        """
        Create an Xarray dataset containing all signals necessary to recreate a GEQDSK file.
        Uses FreeQDSK canonical names. 2D psi grid stored as 'psirz' to distinguish from 'psi' coordinate.
        Data is in COCOS 1.

        Data is fetched from the level1 zarr efm group which contains the full EFM output including
        1D psi profiles (fpol, pres, ffprime, pprime, qpsi) not present in the level2 zarr.

        level1 efm layout (time FIRST for all multi-dim signals):
          0D scalars: (time,)
          1D profiles (fpsi_c, ppsi_c, ffprime, pprime, qpsi_c): (time, psi_norm)
          2D psirz: (time, profile_z, profile_r) -> transposed to (time, profile_r, profile_z)
          boundary (lcfs_r, lcfs_z): (time, lcfs_coords)
        """
        efm_cfg = get_config(Tokamak.MAST).inout.efm_zarr
        efm_url = f"{efm_cfg.endpoint_url}/{efm_cfg.file_path}/{params.shot_id}.{efm_cfg.file_ext}"
        efm = xr.open_zarr(efm_url, group=efm_cfg.group)

        times = params.times
        eq_time = efm.coords["time"].values

        # Grid from psirz dimension coordinates
        # psirz shape: (time, profile_z, profile_r) -> transpose to (time, profile_r, profile_z)
        # so that make_geqdsk_dataset dim labeling (idx, r_grid, z_grid) is shape-consistent
        r_grid = efm.coords["profile_r"].values   # (nr=129,)
        z_grid = efm.coords["profile_z"].values   # (nz=65,)
        psirz_raw = efm["psirz"].values.transpose(0, 2, 1)  # (T_eq, nr, nz)

        # 0D scalars
        rmagx = efm["magnetic_axis_r"].values
        zmagx = efm["magnetic_axis_z"].values
        simagx = efm["psi_axis"].values
        sibdry = efm["psi_boundary"].values
        bcentr = efm["bvac_rmag"].values
        current = efm["plasma_current_c"].values

        # 1D psi profiles: (time, psi_norm) - time already first
        fpol = efm["fpsi_c"].values
        pres = efm["ppsi_c"].values
        ffprime = efm["ffprime"].values
        pprime = efm["pprime"].values
        qpsi = efm["qpsi_c"].values

        # Boundary: (time, lcfs_coords) - time already first
        rbdry_raw = efm["lcfs_r"].values
        zbdry_raw = efm["lcfs_z"].values

        # Limiter (static)
        rlim = efm["limiterr"].values
        zlim = efm["limiterz"].values

        needs_interp = not np.array_equal(times, eq_time)
        if needs_interp:
            rmagx = interp1(eq_time, rmagx, times)
            zmagx = interp1(eq_time, zmagx, times)
            simagx = interp1(eq_time, simagx, times)
            sibdry = interp1(eq_time, sibdry, times)
            bcentr = interp1(eq_time, bcentr, times)
            current = interp1(eq_time, current, times)
            fpol = interp1(eq_time, fpol, times, axis=0)
            pres = interp1(eq_time, pres, times, axis=0)
            ffprime = interp1(eq_time, ffprime, times, axis=0)
            pprime = interp1(eq_time, pprime, times, axis=0)
            qpsi = interp1(eq_time, qpsi, times, axis=0)
            psirz = interp1(eq_time, psirz_raw, times, axis=0)
            rbdry = interp1(eq_time, rbdry_raw, times, axis=0)
            zbdry = interp1(eq_time, zbdry_raw, times, axis=0)
        else:
            psirz = psirz_raw
            rbdry = rbdry_raw
            zbdry = zbdry_raw

        sign_Ip = np.sign(np.nanmedian(current))
        sign_B0 = np.sign(np.nanmedian(bcentr))
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
            times=times,
            r_grid=r_grid,
            z_grid=z_grid,
            rmagx=rmagx,
            zmagx=zmagx,
            simagx=simagx,
            sibdry=sibdry,
            bcentr=bcentr,
            current=current,
            fpol=fpol,
            pres=pres,
            ffprime=ffprime,
            pprime=pprime,
            qpsi=qpsi,
            psirz=psirz,
            rbdry=rbdry,
            zbdry=zbdry,
            cocos_input=cocos_input,
            rlim=rlim,
            zlim=zlim,
        )

        # rmagx excluded to avoid merge conflict with get_efit_parameters
        return ds.drop_vars("rmagx")
