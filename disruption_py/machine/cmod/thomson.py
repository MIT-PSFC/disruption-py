#!/usr/bin/env python3

"""Module for processing Thomson electron density measurements."""

import numpy as np
import scipy as sp

import jax
jax.config.update("jax_enable_x64", True) # The data is in f32, but doing computations in f64 improves stability
import jax.numpy as jnp
import gpjax as gpx

from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
from disruption_py.inout.mds import mdsExceptions
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
            r"\efit_aeqdsk:time", tree_name="_efit_tree"
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


class CmodThomsonGPProfiles:
    """
    Helper class for using Gaussian Processes to fit Thomson scattering profiles.
    """

    @jax.jit
    def gp_profile_single_timestep(
        signal_data: jnp.ndarray,
        signal_error: jnp.ndarray,
        signal_rho: jnp.ndarray,
        result_rho: jnp.ndarray,
        kernel: gpx.kernels.AbstractKernel,
        edge_rho: float,
        edge_value: float,
        edge_error: float,
        epsilon: float,
        outlier_penalty: float,
        outlier_cutoff: float,
        jitter: float,
        ):
        """
        Fit a single timestep using Gaussian Processes.
        Assumes there are no missing / nan values in the input data.

        Since this is a 1D profile of a physical plasma quantity, we can make the following assumptions to constrain the fit:
        1. The signal is always positive.
        2. The signal has zero gradient at the core (rho=0).
        3. The signal goes to some relatively known value at the edge (rho=edge_rho).
        - For ne, Te, J, etc., this value should be 0
        - For q it could be taken as something like 10

        Implementing these assumptions is done as follows:
        1. Use a log-transform on the signal values.
        2. Using JAX, we compute covariences in both the kernel function and its derivative, and add a low-uncertainty data point at (0, signal(0)).
        3. A synthetic data point is added at (edge_rho, edge_value) with uncertainty edge_error.

        The final piece which makes the profiles 'look good' is fitting in rho^2 space.
        This essentially spreads out the edge region, allowing the GP to vary more quickly there while remaining smooth in the core.

        Parameters
        ----------
        signal_data : jnp.ndarray
            The signal data to fit.
        signal_error : jnp.ndarray
            The standard deviations associated with the signal data.
        signal_rho : jnp.ndarray
            The normalized radius (rho) values corresponding to the signal data.
        result_rho : jnp.ndarray
            The rho values where the fitted profile should be evaluated.
        kernel : gpx.kernels.AbstractKernel
            The kernel to use for the Gaussian Process. Must be initialized outside this jitted function.
        edge_rho : float
            The rho value at the edge where the signal is expected to approach edge_value
        edge_value : float
            The expected signal value at edge_rho.
        edge_error : float
            The uncertainty associated with the edge_value.
        epsilon : float
            A small value added to avoid log(0) issues.
        outlier_penalty : float
            The error value to assign to identified outliers to effectively ignore them.
        outlier_cutoff : float
            The rho value below which outliers will be identified.
        jitter : float
            A small jitter value added to the diagonal of the covariance matrix for numerical stability.

        Returns
        -------
        jnp.ndarray, jnp.ndarray
            The fitted signal profile and its standard deviation evaluated at result_rho.
        """

        # Ensure inputs are jnp arrays with a flat shape. They will be reshaped as needed for given matrix operations.
        signal_data = jnp.array(signal_data).flatten()
        signal_error = jnp.array(signal_error).flatten()
        signal_rho = jnp.array(signal_rho).flatten()
        result_rho = jnp.array(result_rho).flatten()

        # Handle NaN and invalid values by replacing them with dummy values and giving them huge error bars
        # This keeps array shapes constant for JIT compilation
        valid_mask = ~jnp.isnan(signal_data) & ~jnp.isnan(signal_error) & ~jnp.isnan(signal_rho) & (signal_data > 0) & (signal_error > 0)

        # Replace invalid data with 0 at the edge point
        signal_data = jnp.where(valid_mask, signal_data, edge_value)
        signal_error = jnp.where(valid_mask, signal_error, outlier_penalty)  # Give invalid points huge error
        signal_rho = jnp.where(valid_mask, signal_rho, edge_rho)  # Dummy rho value for invalid points

        # Normalize inputs by average value to improve numerical stability
        # Use nanmean to handle any remaining NaNs
        norm_factor = jnp.mean(signal_data)
        observations_norm = signal_data / norm_factor
        errors_norm = signal_error / norm_factor

        # Transform to log-space to enforce positivity
        observations_log = jnp.log(observations_norm + epsilon)
        # Error in log space std_y = |dy/dx| * std_x = (1/x) * std_x
        # Intuition is that in log space we care about relative error, not absolute error
        errors_log = errors_norm / (observations_norm + epsilon)

        # Transform rho to rho^2 space to allow more variation near edge while maintaining smoothness in core
        rho2 = signal_rho ** 2

        # Perform a quick GP fit without boundary conditions to identify outliers
        K_rough = kernel.gram(rho2.reshape(-1,1)).to_dense()  # K(X, X), prior covariance at observed points
        noise_rough = jnp.diag(errors_log ** 2)  # Observation noise, GPJax uses variance
        K_total_rough = K_rough + noise_rough + jnp.eye(len(rho2)) * jitter  # Add small jitter for numerical stability
        # Make predictions at observed points
        alpha = jnp.linalg.solve(K_total_rough, observations_log)
        predictions_rough = K_rough @ alpha
        # Compute prediction variance at training points
        predictions_var = jnp.diag(K_rough - K_rough @ jnp.linalg.solve(K_total_rough, K_rough))
        predictions_std = jnp.sqrt(jnp.maximum(predictions_var, 0))
        # Identify outliers as core points where residual > 3 * predicted stddev
        residuals = jnp.abs(observations_log - predictions_rough)
        outlier_mask = (residuals > (3 * predictions_std)) & (signal_rho < outlier_cutoff)
        
        # Make new arrays without outliers
        errors_clean = jnp.where(outlier_mask, outlier_penalty, errors_log)
        observations_clean = observations_log  # Value won't matter due to large error
        rho2_clean = rho2

        # Add synthetic edge point
        y = jnp.concatenate(
            [observations_clean, jnp.array([jnp.log(edge_value / norm_factor + epsilon)])]
        )
        y_std = jnp.concatenate(
            [errors_clean, jnp.array([edge_error / norm_factor / (edge_value / norm_factor + epsilon)])]
        )
        x = jnp.concatenate(
            [rho2_clean, jnp.array([edge_rho**2])]
        )

        # Now it's time to get funky
        x_deriv_bc = jnp.array([0.0])  # Derivative boundary condition at rho=0
        y_deriv_bc = jnp.array([0.0])  # We want d(signal)/d(rho) = 0 at rho=0
        y_deriv_bc_std = jnp.array([1e-6])  # Small uncertainty on derivative BC, since the gradient really should be 0 there
        
        def _kernel_dx2(kernel, x1, x2):
            """Compute dk(x1, x2)/dx2 with autodiff"""
            def k_func(x2_arg):
                return kernel(x1, x2_arg).squeeze()
            grad_fn = jax.grad(k_func)
            return grad_fn(x2)
        
        def _kernel_dx1_dx2(kernel, x1, x2):
            """Compute d2k(x1, x2)/dx1dx2 with autodiff"""
            def dk_dx2(x1_arg):
                def k_func(x2_arg):
                    return kernel(x1_arg, x2_arg).squeeze()
                grad_fn = jax.grad(k_func)
                return grad_fn(x2).squeeze()
            grad_fn = jax.grad(dk_dx2)
            return grad_fn(x1)
        
        # Build augmented kernel matrix
        # K = [[K_ff,  K_fd],
        #      [K_df,  K_dd]]
        # where f = function values, d = derivatives

        # K_ff: value-value covariances
        K_ff = kernel.gram(x.reshape(-1, 1)).to_dense()
        # K_fd: value-derivative cross-covariances
        K_fd = jax.vmap(lambda x1: jax.vmap(lambda x2: _kernel_dx2(kernel, x1, x2))(x_deriv_bc.reshape(-1, 1)))(x.reshape(-1, 1)).reshape(-1, 1)
        # K_df: derivative-value cross-covariances (just the transpose of K_fd)
        K_df = K_fd.T
        # K_dd: derivative-derivative covariances
        K_dd = jax.vmap(lambda x1: jax.vmap(lambda x2: _kernel_dx1_dx2(kernel, x1, x2))(x_deriv_bc.reshape(-1, 1)))(x_deriv_bc.reshape(-1, 1)).reshape(-1, 1)

        K_aug = jnp.block([
            [K_ff, K_fd],
            [K_df, K_dd]
        ])
        
        # Add observation noise
        noise_aug = jnp.diag(jnp.concatenate([y_std**2, y_deriv_bc_std**2]))
        K_total_aug = K_aug + noise_aug + jnp.eye(K_aug.shape[0]) * jitter  # Add small jitter for numerical stability

        # Prepare augmented observations
        y_aug = jnp.concatenate([y, y_deriv_bc])

        # Rho values to predict at
        x_pred = result_rho**2

        # Build cross-covariance matrix between prediction points and training points
        K_star_f = kernel.cross_covariance(x_pred.reshape(-1, 1), x.reshape(-1, 1)) 
        K_star_d = jax.vmap(lambda x1: jax.vmap(lambda x2: _kernel_dx2(kernel, x1, x2).squeeze())(x_deriv_bc.reshape(-1, 1)))(x_pred.reshape(-1, 1)).reshape(x_pred.shape[0], 1)
        K_star = jnp.hstack([K_star_f, K_star_d])

        # Make predictions at new points
        alpha_aug = jnp.linalg.solve(K_total_aug, y_aug)
        y_pred = K_star @ alpha_aug
        # Compute prediction variance at training points
        y_pred_var = jnp.diag(
            kernel.gram(x_pred.reshape(-1, 1)).to_dense()
            - K_star @ jnp.linalg.solve(K_total_aug, K_star.T)
        )
        y_error = jnp.sqrt(jnp.maximum(y_pred_var, 0))

        # Transform back to normal space in rho and signal
        signal_pred = jnp.exp(y_pred) * norm_factor - epsilon
        signal_error = y_error * jnp.exp(y_pred) * norm_factor

        # If insufficient valid data (< 3 points), return NaN
        n_valid = jnp.sum(valid_mask)
        signal_pred = jnp.where(n_valid >= 3, signal_pred, jnp.nan)
        signal_error = jnp.where(n_valid >= 3, signal_error, jnp.nan)

        return signal_pred, signal_error

    def gp_profile(
        signal_data, 
        signal_error, 
        signal_rho, 
        result_rho, 
        kernel=gpx.kernels.Matern32(),
        edge_rho=1.2,
        edge_value=1e-3,
        edge_error=0.5,
        epsilon=1e-3,
        outlier_penalty=1e6,
        outlier_cutoff=0.8,
        jitter=1e-6,
    ):
        """
        Fit Gaussian Process profiles over multiple timesteps.

        Parameters
        ----------
        signal_data : jnp.ndarray
            The signal data to fit. Expected shape is (n_times, n_measurements).
        signal_error : jnp.ndarray
            The standard deviations associated with the signal data. Expected shape is (n_times, n_measurements).
        signal_rho : jnp.ndarray
            The normalized radius (rho) values corresponding to the signal data. Expected shape is (n_times, n_measurements).
        result_rho : jnp.ndarray
            The rho values where the fitted profile should be evaluated. Expected shape is (n_result_points,).
        kernel : gpx.kernels.Kernel, optional
            The kernel to use for the Gaussian Process (default is Matern32).
        edge_rho : float, optional
            The rho value at the edge where the signal is expected to approach edge_value (default is 1.2).
        edge_value : float, optional
            The expected signal value at edge_rho (default is 1e-3).
        edge_error : float, optional
            The uncertainty associated with the edge_value (default is 0.5).
        epsilon : float, optional
            A small value added to avoid log(0) issues (default is 1e-3).
        outlier_penalty : float, optional
            The error value to assign to identified outliers to effectively ignore them (default is 1e6).
        outlier_cutoff : float, optional
            The rho value below which outliers will be identified (default is 0.8).
        jitter : float, optional
            A small jitter value added to the diagonal of the covariance matrix for numerical stability (default is 1e-6).

        Returns
        -------
        jnp.ndarray, jnp.ndarray
            The fitted signal profiles and their standard deviations evaluated at result_rho.
            Shapes are both (n_times, n_result_points).
        """

        signal_pred_all, signal_error_all = jax.vmap(
            CmodThomsonGPProfiles.gp_profile_single_timestep,
            in_axes=(0, 0, 0, None, None, None, None, None, None, None, None, None),  # 12 args: data, error, rho (all vmapped), result_rho, edge_rho, edge_value, edge_error, epsilon, kernel
            out_axes=(0, 0)
        )(
            signal_data,
            signal_error,
            signal_rho,
            result_rho,
            kernel,
            edge_rho,
            edge_value,
            edge_error,
            epsilon,
            outlier_penalty,
            outlier_cutoff,
            jitter,
        )

        return signal_pred_all, signal_error_all

        
        