#!/usr/bin/env python3

import logging
import sys
import traceback
import warnings

import numpy as np
import scipy as sp

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
from disruption_py.machine.cmod.efit import CmodEfitMethods
from disruption_py.machine.tokamak import Tokamak


class CmodDraftPhysicsMethods:

    @staticmethod
    def get_edge_parameters(times, p_Te, p_ne, edge_rho_min=0.85, edge_rho_max=0.95):
        """
        Compute the edge Temperature and edge Density signal from the TS.

        Parameters
        ----------
        times : array_like
            The times at which to calculate the edge parameters.
        p_Te : BivariatePlasmaProfile
            The Te measurements [keV] in terms of the time and rho of the measurment.
        p_ne : BivariatePlasmaProfile
            The ne measurements [keV] in terms of the time and rho of the measurment.
        edge_rho_min : float [0,1]
            The rho that defines the minimum of the "edge" region
        edge_rho_max : float [0,1]
            The rho that defines the maximum of the "edge" region

        Returns
        -------
        Te_edge : array_like
            The edge temperature (averaged over the edge region) on the requested
            timebase.
        ne_edge : array_like
            The edge density (averaged over the edge region) on the requested timebase.

        Original Authors
        ----------------
        Andrew Maris (maris@mit.edu)


        """
        # Base of rho to interpolate onto
        rhobase = np.arange(0, 1, 0.001)

        # Linear interpolate on time and rho
        Te_interpolator = sp.interpolate.LinearNDInterpolator(
            (p_Te.X[:, 0], p_Te.X[:, 1]), p_Te.y
        )
        ne_interpolator = sp.interpolate.LinearNDInterpolator(
            (p_ne.X[:, 0], p_ne.X[:, 1]), p_ne.y
        )

        # Create mesh to compute interpolation over
        timebase_mesh, rhobase_mesh = np.meshgrid(times, rhobase)
        # Compute interpolated values
        Te_interp = Te_interpolator(timebase_mesh, rhobase_mesh)
        ne_interp = ne_interpolator(timebase_mesh, rhobase_mesh)

        plotting = False
        if plotting:
            import matplotlib.pyplot as plt

            plt.ion()
            plt.pcolormesh(timebase_mesh, rhobase_mesh, Te_interp, shading="auto")
            plt.plot(p_Te.X[:, 0], p_Te.X[:, 1], "ok", label="input point")
            plt.legend()
            plt.colorbar()
            plt.show(block=True)

        # Compute Te_edge
        # Make mask for rho in edge region
        rhobase_mesh_mask = (rhobase_mesh >= edge_rho_min) & (
            rhobase_mesh <= edge_rho_max
        )

        # Assert that rho values are indeed in desired range
        assert np.all(rhobase_mesh[rhobase_mesh_mask] >= edge_rho_min)
        assert np.all(rhobase_mesh[rhobase_mesh_mask] <= edge_rho_max)

        # Use mask to get only edge values
        Te_interp_edge = np.where(rhobase_mesh_mask, Te_interp, np.nan)
        ne_interp_edge = np.where(rhobase_mesh_mask, ne_interp, np.nan)

        # Compute edge quantities
        # Catch warning about taking nanmean of an empty array. This is ok because
        # we want it to return nan for empty arrays
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            Te_edge = np.nanmean(Te_interp_edge, axis=0)
            ne_edge = np.nanmean(ne_interp_edge, axis=0)
        assert len(Te_edge) == len(times)
        assert len(ne_edge) == len(times)

        if plotting:
            plt.plot(times, Te_edge, label="Te_edge")
            plt.plot(times, ne_edge, label="ne_edge")
            plt.legend()
            plt.show(block=True)

        output = {"te_edge": Te_edge, "ne_edge": ne_edge}
        return output

    @staticmethod
    @physics_method(
        columns=["te_edge", "ne_edge"],
        tokamak=Tokamak.CMOD,
    )
    def _get_edge_parameters(params: PhysicsMethodParams):
        nan_output = {
            "te_edge": [np.nan],
            "ne_edge": [np.nan],
        }
        try:
            # sys.path.append("/home/sciortino/usr/python3modules/eqtools3")
            sys.path.append("/home/sciortino/usr/python3modules/profiletools3")
            sys.path.append("/home/sciortino/usr/python3modules/gptools3")
            # import eqtools
            import profiletools
        except Exception as e:
            logging.warning("Could not import profiletools or eqtools")
            logging.debug(traceback.format_exc())
            pass

        # Ignore shots on the blacklist
        if CmodPhysicsMethods.is_on_blacklist(params.shot_id):
            return nan_output
        # Range of rho to interpolate over
        rhobase = np.arange(0, 1, 0.001)
        # Get mina and max time from TS tree
        node_path = ".yag_new.results.profiles"
        try:
            (ts_time,) = params.mds_conn.get_dims(
                node_path + ":te_rz", tree_name="electrons"
            )
        except:
            return nan_output

        t_min = np.max([0.1, np.min(ts_time)])
        t_max = np.max(ts_time)

        # Get core and edge Thomson profiles over rho := sqrtpsinorm
        p_Te = profiletools.Te(
            params.shot_id,
            include=["CTS", "ETS"],
            abscissa="sqrtpsinorm",
            t_min=t_min,
            t_max=t_max,
            remove_zeros=True,
        )
        p_ne = profiletools.ne(
            params.shot_id,
            include=["CTS", "ETS"],
            abscissa="sqrtpsinorm",
            t_min=t_min,
            t_max=t_max,
            remove_zeros=True,
        )

        # try:
        #    equal_R = p_ne.X[:,1] == p_Te.X[:,1]
        #    assert np.sum(equal_R) == len(p_ne.X[:,1])
        # except:
        #    raise ValueError('Edge Thomson rhobase differs between ne and Te')
        #    return None, None

        # consider only flux surface on which points were measured, regardless of
        # LFS or HFS
        p_Te.X = np.abs(p_Te.X)
        p_ne.X = np.abs(p_ne.X)

        # set some minimum uncertainties. Recall that units in objects are 1e20m^{-3}
        # and keV
        p_ne.y[p_ne.y <= 0.0] = 0.01  # 10^18 m^-3
        p_Te.y[p_Te.y <= 0.01] = 0.01  # 10 eV
        p_ne.err_y[p_ne.err_y <= 0.01] = 0.01  # 10^18 m^-3
        p_Te.err_y[p_Te.err_y <= 0.02] = 0.02  # 20 eV

        # points in the pedestal that have x uncertainties larger than 0.1 don't help at all
        # do this filtering here because filtering of err_X only works before time-averaging
        p_ne.remove_points(np.logical_and(p_ne.X[:, 1] >= 0.85, p_ne.err_X[:, 1] > 0.1))
        p_Te.remove_points(np.logical_and(p_Te.X[:, 1] >= 0.85, p_Te.err_X[:, 1] > 0.1))

        # cleanup of low Te values
        # TS Te should be >15 eV inside near SOL
        p_Te.remove_points(np.logical_and(p_Te.X[:, 0] < 1.03, p_Te.y < 0.015))

        output = CmodPhysicsMethods.get_edge_parameters(params.times, p_Te, p_ne)
        return output

    @staticmethod
    def get_H98():
        pass

    # TODO: Finish
    @staticmethod
    @physics_method(
        columns=["h98", "wmhd", "btor", "dwmhd_dt", "p_input"],
        tokamak=Tokamak.CMOD,
    )
    def _get_H98(params: PhysicsMethodParams):
        """
        Prepare to compute H98 by getting tau_E

        Scaling from eq. 20, ITER Physics Basis Chapter 2
        https://iopscience.iop.org/article/10.1088/0029-5515/39/12/302/pdf
        (in s, MA, T, MW, 10^19 m^−3, AMU, m)
        Original Authors
        ----------------
        Andrew Maris (maris@mit.edu)

        """

        # Get parameters for calculating confinement time
        powers_df = CmodPhysicsMethods._get_power(params=params)
        efit_df = CmodEfitMethods.get_efit_parameters(params=params)
        density_df = CmodPhysicsMethods._get_densities(params=params)
        ip_df = CmodPhysicsMethods._get_ip_parameters(params=params)

        # Get BT

        btor, t_mag = params.mds_conn.get_data_with_dims(
            r"\btor", tree_name="magnetics"
        )  # tmag: [s]
        # Toroidal power supply takes time to turn on, from ~ -1.8 and should be
        # on by t=-1. So pick the time before that to calculate baseline
        baseline_indices = np.where(t_mag <= -1.8)
        btor = btor - np.mean(btor[baseline_indices])
        btor = np.abs(interp1(t_mag, btor, params.times))

        ip = np.abs(ip_df.ip) / 1.0e6  # [A] -> [MA]
        n_e = density_df.n_e / 1.0e19  # [m^-3] -> [10^19 m^-3]
        p_input = powers_df.p_input / 1.0e6  # [W] -> [MW]
        dwmhd_dt = efit_df.dwmhd_dt / 1.0e6  # [W] -> [MW]
        wmhd = efit_df.wmhd / 1.0e6  # [J] -> [MJ]
        R0 = efit_df.rmagx / 100  # [cm] -> [m]
        # Estimate confinement time
        tau = wmhd / (p_input - dwmhd_dt)

        # Compute 1998 tau_E scaling, taking A (atomic mass) = 2
        tau_98 = (
            0.0562
            * (n_e**0.41)
            * (2**0.19)
            * (ip**0.93)
            * (R0**1.39)
            * (efit_df.a_minor**0.58)
            * (efit_df.kappa**0.78)
            * (btor**0.15)
            * (p_input**-0.69)
        )
        H98 = tau / tau_98
        output = {
            "h98": H98,
            "wmhd": wmhd,
            "btor": btor,
            "dwmhd_dt": dwmhd_dt,
            "p_input": p_input,
        }
        return output

    # TODO: Calculate v_mid
    @staticmethod
    @physics_method(columns=["v_0"], tokamak=Tokamak.CMOD)
    def get_rotation_velocity(params: PhysicsMethodParams):
        nan_output = {"v_0": [np.nan]}
        data = resources.files(disruption_py.data)
        file = data.joinpath("lock_mode_calib_shots.txt")
        with resources.as_file(file) as fio:
            calibrated = pd.read_csv(fio)
        # Check to see if shot was done on a day where there was a locked
        # mode HIREX calibration by cross checking with list of calibrated
        # runs. If not calibrated, return NaN outputs.
        if params.shot_id not in calibrated:
            return nan_output
        try:
            intensity, time = params.mds_conn.get_data_with_dims(
                ".hirex_sr.analysis.a:int", tree_name="spectroscopy", astype="float64"
            )
            vel, hirextime = params.mds_conn.get_data_with_dims(
                ".hirex_sr.analysis.a:vel", tree_name="spectroscopy", astype="float64"
            )
        except mdsExceptions.TreeFOPENR:
            params.logger.warning(
                "[Shot %s]: Failed to open necessary trees for rotational velocity calculations.",
                params.shot_id,
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
            return nan_output
        output = CmodPhysicsMethods._get_rotation_velocity(
            params.times, intensity, time, vel, hirextime
        )
        return output
