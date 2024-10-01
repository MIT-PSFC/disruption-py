#!/usr/bin/env python3

"""Module for processing Thomson electron density measurements."""

import numpy as np
import scipy as sp
from MDSplus import mdsExceptions

from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
from disruption_py.core.utils.misc import safe_cast
from disruption_py.machine.cmod.efit import CmodEfitMethods


class CmodThomsonDensityMeasure:
    """
    A helper class for Thomson electron density measurements.
    """

    @staticmethod
    def compare_ts_tci(params: PhysicsMethodParams, nlnum=4):
        """
        Helper function used for comparing electron density measurements from
        Thomson scattering with the two-color interferometer (TCI).

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information.
        nlnum : int, optional
            The number of the TCI channel to compare (default is 4).

        Returns
        -------
        tuple
            Tuple containing the Thomson and TCI electron density measurements.
        """
        nl_ts1 = [1e32]
        nl_ts2 = [1e32]
        nl_tci1 = [1e32]
        nl_tci2 = [1e32]
        ts_time1 = [1e32]
        ts_time2 = [1e32]
        (tci_time,) = params.mds_conn.get_dims(
            ".YAG_NEW.RESULTS.PROFILES:NE_RZ", tree_name="electrons"
        )
        tci, tci_t = params.mds_conn.get_data_with_dims(
            f".TCI.RESULTS:NL_{nlnum:02d}", tree_name="electrons"
        )
        nlts, nlts_t = CmodThomsonDensityMeasure._integrate_ts_tci(params, nlnum)
        t0 = np.amin(nlts_t)
        t1 = np.amax(nlts_t)
        nyag1, nyag2, indices1, indices2 = CmodThomsonDensityMeasure._parse_yags(params)
        time1, time2 = -1, -1
        if nyag1 > 0:
            ts_time1 = tci_time[indices1]
            (valid_indices,) = np.where((ts_time1 >= t0) & (ts_time1 <= t1))
            if valid_indices.size > 0:
                nl_tci1 = interp1(tci_t, tci, ts_time1[valid_indices])
                nl_ts1 = interp1(nlts_t, nlts, ts_time1[valid_indices])
                time1 = ts_time1[valid_indices]
        if nyag2 > 0:
            ts_time2 = tci_time[indices2]
            (valid_indices,) = np.where((ts_time2 >= t0) & (ts_time2 <= t1))
            if valid_indices.size > 0:
                nl_tci1 = interp1(tci_t, tci, ts_time2[valid_indices])
                nl_ts1 = interp1(nlts_t, nlts, ts_time2[valid_indices])
                time2 = ts_time2[valid_indices]
        return nl_ts1, nl_ts2, nl_tci1, nl_tci2, time1, time2

    @staticmethod
    def _parse_yags(params: PhysicsMethodParams):
        """
        Parse YAG laser data to determine indices and counts.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information.

        Returns
        -------
        tuple
            Tuple containing the number of YAGs and their indices.
        """
        nyag1 = params.mds_conn.get_data(r"\knobs:pulses_q", tree_name="electrons")
        nyag2 = params.mds_conn.get_data(r"\knobs:pulses_q_2", tree_name="electrons")
        indices1 = -1
        indices2 = -1
        dark = params.mds_conn.get_data(r"\n_dark_prior", tree_name="electrons")
        ntotal = params.mds_conn.get_data(r"\n_total", tree_name="electrons")
        nt = ntotal - dark
        if nyag1 == 0:
            if nyag2 != 0:
                indices2 = np.arange(nyag2)
        else:
            if nyag2 == 0:
                indices1 = np.arange(nyag1)
            else:
                if nyag1 == nyag2:
                    indices1 = 2 * np.arange(nyag1)
                    indices2 = indices1 + 1
                else:
                    if nyag1 == nyag2:
                        indices1 = 2 * np.arange(nyag1)
                        indices2 = indices1 + 1
                    else:
                        indices1 = 2 * np.arange(nyag1) + (nyag1 > nyag2)
                        indices2 = np.concatenate(
                            (
                                2 * np.arange(nyag2) + (nyag1 < nyag2),
                                2 * nyag2 + np.arange(nyag1 - nyag2 - 1),
                            )
                        )
        (v_ind1,) = np.where(indices1 < nt)
        if nyag1 > 0 and v_ind1.size > 0:
            indices1 = indices1[v_ind1]
        else:
            indices1 = -1
        (v_ind2,) = np.where(indices2 < nt)
        if nyag2 > 0 and v_ind2.size > 0:
            indices2 = indices2[v_ind2]
        else:
            indices2 = -1
        return nyag1, nyag2, indices1, indices2

    @staticmethod
    def _integrate_ts_tci(params: PhysicsMethodParams, nlnum):
        """
        Integrate Thomson electron density measurement to the line integrated electron
        density for comparison with two-color interferometer (TCI) measurement results.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information.
        nlnum : int
            The number of the TCI channel to integrate.

        Returns
        -------
        tuple
            Tuple containing the integrated electron density and corresponding times.
        """
        t, z, n_e, n_e_sig = CmodThomsonDensityMeasure._map_ts2tci(params, nlnum)
        if z[0, 0] == 1e32:
            return None, None  # TODO: Log and maybe return nan arrs
        nts = len(t)
        nlts_t = t
        nlts = np.full(t.shape, np.nan)
        for i in range(nts):
            (ind,) = np.where(
                (np.abs(z[i, :]) < 0.5)
                & (n_e[i, :] > 0)
                & (n_e[i, :] < 1e21)
                & (n_e[i, :] / n_e_sig[i, :] > 2)
            )
            if len(ind) < 3:
                nlts[i] = 0
            else:
                x = z[i, ind]
                y = n_e[i, ind]
                _, ind_uniq = np.unique(x, return_index=True)
                y = y[ind_uniq]
                nlts[i] = np.trapz(y, x)
        return nlts, nlts_t

    @staticmethod
    def _map_ts2tci(params: PhysicsMethodParams, nlnum):
        """
        Map Thomson density measurements to TCI measurements.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information.
        nlnum : int
            The number of the TCI channel to map.

        Returns
        -------
        tuple
            Tuple containing time, z-coordinates, electron density, and density errors.
        """
        core_mult = 1.0
        edge_mult = 1.0
        t = [1e32]
        z = [1e32]
        n_e = [1e32]
        n_e_sig = [1e32]
        flag = 1
        valid_indices, efit_times = CmodEfitMethods.efit_check(params)
        ip = params.mds_conn.get_data(r"\ip", "cmod")
        if np.mean(ip) > 0:
            flag = 0
        efit_times = params.mds_conn.get_data(
            r"\efit_aeqdsk:time", tree_name="_efit_tree", astype="float64"
        )
        t1 = np.amin(efit_times)
        t2 = np.amax(efit_times)
        psia, psia_t = params.mds_conn.get_data_with_dims(
            r"\efit_aeqdsk:SIBDRY", tree_name="_efit_tree"
        )
        psi_0 = params.mds_conn.get(r"\efit_aeqdsk:SIMAGX", tree_name="_efit_tree")
        nets_core, nets_core_t = params.mds_conn.get_data_with_dims(
            ".YAG_NEW.RESULTS.PROFILES:NE_RZ", tree_name="electrons"
        )
        nets_core_err = params.mds_conn.get_data(
            ".YAG_NEW.RESULTS.PROFILES:NE_ERR", tree_name="electrons"
        )
        zts_core = params.mds_conn.get_data(
            ".YAG_NEW.RESULTS.PROFILES:Z_SORTED", tree_name="electrons"
        )
        mts_core = len(zts_core)
        zts_edge = params.mds_conn.get_data(r"\fiber_z")
        mts_edge = len(zts_edge)
        try:
            nets_edge = params.mds_conn.get_data(r"\ts_ne")
            nets_edge_err = params.mds_conn.get_data(r"\ts_ne_err")
        except mdsExceptions.MdsException:
            nets_edge = np.zeros((len(nets_core[:, 1]), mts_edge))
            nets_edge_err = nets_edge + 1e20
        mts = mts_core + mts_edge
        rts = params.mds_conn.get(".YAG.RESULTS.PARAM:R") + np.zeros((1, mts))
        rtci = params.mds_conn.get_data(".tci.results:rad")
        nts = len(nets_core_t)
        zts = np.zeros((1, mts))
        zts[:, :mts_core] = zts_core
        zts[:, mts_core:] = zts_edge
        nets = np.zeros((nts, mts))
        nets_err = np.zeros((nts, mts))
        nets[:, :mts_core] = (nets_core * core_mult).T
        nets_err[:, :mts_core] = (nets_core_err * core_mult).T
        nets[:, mts_core:] = (nets_edge * edge_mult).T
        nets_err[:, mts_core:] = (nets_edge_err * edge_mult).T
        valid_indices = np.where((nets_core_t >= t1) & (nets_core_t <= t2))
        if len(valid_indices) == 0:
            return t, z, n_e, n_e_sig
        nets_core_t = nets_core_t[valid_indices]
        nets = nets[valid_indices]
        nets_err = nets_err[valid_indices]
        psits = CmodThomsonDensityMeasure._efit_rz2psi(params, rts, zts, nets_core_t)
        mtci = 101
        ztci = -0.4 + 0.8 * np.arange(0, mtci) / (mtci - 1)
        rtci = rtci[nlnum] + np.zeros((1, mtci))
        psitci = CmodThomsonDensityMeasure._efit_rz2psi(params, rtci, ztci, nets_core_t)
        psia = interp1(psia_t, psia, nets_core_t)
        psi_0 = interp1(psia_t, psi_0, nets_core_t)
        nts = len(nets_core_t)
        for i in range(nts):
            psits[:, i] = (psits[:, i] - psi_0[i]) / (psia[i] - psi_0[i])
            psitci[:, i] = (psitci[:, i] - psi_0[i]) / (psia[i] - psi_0[i])
        zmapped = np.zeros((nts, 2 * mts)) + 1e32
        nemapped = zmapped.copy()
        nemapped_err = zmapped.copy()
        for i in range(nts):
            index = np.argmin(psitci[i, :]) if flag else np.argmax(psitci[i, :])
            psi_val = psitci[i, index]
            for j in range(mts):
                if (flag and psits[j, i] >= psi_val) or (
                    not flag and psits[j, i] <= psi_val
                ):
                    a1 = interp1(psitci[:index, i], ztci[:index], psits[j, i])
                    a2 = interp1(psitci[index:, i], ztci[index:], psits[j, i])
                    zmapped[i, [j, j + mts]] = [a1, a2]
                    nemapped[i, [j, j + mts]] = nets[i, j]
                    nemapped_err[i, [j, j + mts]] = nets_err[i, j]
            sorted_indices = np.argsort(zmapped[i, :])
            zmapped[i, :] = zmapped[i, sorted_indices]
            nemapped[i, :] = nemapped[i, sorted_indices]
            nemapped_err[i, :] = nemapped_err[i, sorted_indices]
        z = zmapped
        n_e = nemapped
        n_e_sig = nemapped_err
        t = nets_core_t
        return t, z, n_e, n_e_sig

    @staticmethod
    def _efit_rz2psi(params: PhysicsMethodParams, r, z, t, tree="analysis"):
        """
        Interpolate the magnetic flux function (psi) from R and Z coordinates.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information.
        r : array_like
            Radial coordinates.
        z : array_like
            Vertical coordinates.
        t : array_like
            Time points for interpolation.
        tree : str, optional
            The MDSplus tree name (default is "analysis").

        Returns
        -------
        np.ndarray
            Interpolated psi values corresponding to the input R and Z coordinates.
        """
        r = r.flatten()
        z = z.flatten()
        psi = np.full((len(r), len(t)), np.nan)
        z = safe_cast(z, "float32")
        psirz, rgrid, zgrid, times = params.mds_conn.get_data_with_dims(
            r"\efit_geqdsk:psirz", tree_name=tree, dim_nums=[0, 1, 2]
        )
        rgrid, zgrid = np.meshgrid(rgrid, zgrid)

        points = np.array(
            [rgrid.flatten(), zgrid.flatten()]
        ).T  # This transposes the array to shape (n, 2)
        for i, time in enumerate(t):
            # Find the index of the closest time
            time_idx = np.argmin(np.abs(times - time))
            # Extract the corresponding psirz slice and transpose it
            psirz = np.transpose(psirz[time_idx, :, :])
            # Perform cubic interpolation on the psirz slice
            values = psirz.flatten()
            psi[:, i] = sp.interpolate.griddata(points, values, (r, z), method="cubic")

        return psi
