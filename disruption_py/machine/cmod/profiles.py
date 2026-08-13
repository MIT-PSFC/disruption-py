#!/usr/bin/env python3

"""
Module for the C-Mod edge/core profile machinery.

Implements:
- `_EfitFluxGrid`: (R, Z, t) -> rho and (R, Z, t) -> |B_T| mappings from the EFIT flux grid
- `_TsProfiles`: Combined core+edge Thomson and GPC2 ECE profile readers
- `CmodEdgeCoreProfiles`:
    * Raw edge/core-average and Osborne modified-tanh profile fitting
    * Dimensionless variables calculation:
        ** collisionality nu*,
        ** normalized gyroradius rho*,
        ** toroidal beta beta_T
    Check Maris 2025, eqs. (2)-(4), plus the L-mode density limit proximity metric of
    eq. (7).

"Edge" is the band rho in [0.85, 0.95] with rho = sqrt(psi_norm)
(paper section 2.2); "core" is the analogous band rho in [0.0, 0.2].

The decorated `@physics_method` wrappers live in
`disruption_py.machine.cmod.physics`; like `disruption_py.machine.cmod.thomson`,
this module is deliberately not a method holder.

Units: densities in m^-3, temperatures in eV (converted to J internally),
pressures in Pa, magnetic fields in T, currents in A, lengths in m.

References
-------
- publication: AD Maris, et al. (2025), "Correlation of the L-mode density
limit with edge collisionality", Nuclear Fusion 65 016051, DOI: [10.1088/
1741-4326/ad90f0](https://doi.org/10.1088/1741-4326/ad90f0)
"""

from dataclasses import dataclass, field

import numpy as np
import scipy.constants as const
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import curve_fit

from disruption_py.core.physics_method.errors import (
    CalculationError,
    MismatchCalculationError,
)
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
from disruption_py.inout.mds import mdsExceptions

# rho = sqrt(psi_norm) bands. The edge band is Maris 2025 section 2.2; the core
# band mirrors the existing DisruptionPy core convention (core_bound_factor =
# 0.2 in `CmodPhysicsMethods._get_te_profile_params_ece`).
RHO_BANDS: dict[str, tuple[float, float]] = {
    "edge": (0.85, 0.95),
    "core": (0.0, 0.2),
}

# Preferred EFIT trees. efit21 is the latest C-Mod reconstruction (high time
# resolution, ~1700 slices vs ~90 for analysis); analysis is the robust
# fallback. Pinned by literal name rather than the "_efit_tree" nickname so
# that the flux grid, the geometry, and B_T come from one self-consistent
# reconstruction regardless of the `efit_nickname_setting`.
EFIT_TREE_PREFERENCE: tuple[str, ...] = ("efit21", "analysis")

# Deuterium ion mass; Maris 2025 assumes a deuterium majority plasma.
M_DEUTERON = const.physical_constants["deuteron mass"][0]  # [kg]

# Plasma currents below this magnitude make the cylindrical safety factor (and
# hence nu*) meaningless; such time slices are returned as NaN.
MIN_IP = 1e4  # [A]

# C-Mod MDSplus node paths.
_TS_CORE_NODE = ".yag_new.results.profiles"
_TS_CORE_R = ".yag.results.param:r"  # vertical-system major radius [m]
_ECE_NODE = ".gpc_2.results"
_EFIT_GEQDSK = r"\efit_geqdsk"
_EFIT_AEQDSK = r"\efit_aeqdsk"

# LSVM-D boundary coefficients, Maris 2025 eq. (7): nu*_edge = 3.5 *
# beta_T,edge^-0.40, with beta_T,edge as a fraction (NOT a percentage).
LSVM_D_COEFF = 3.5
LSVM_D_EXPONENT = -0.40


@dataclass
class _EfitFluxGrid:
    """
    EFIT equilibrium quantities needed for the rho and |B_T| mappings.

    Arrays use SI units. `psirz` is stored as `[nt, nz, nr]` (poloidal flux
    per radian, the EFIT convention; only differences/normalizations are used,
    so the COCOS sign and 2*pi factor cancel). `fpol` is `F = R*B_phi` on a
    uniform normalized-poloidal-flux grid `[0, 1]` of length `npsi`
    (`[nt, npsi]`). Array layouts follow the EFIT g-file convention; the whole
    grid is loaded eagerly by `CmodEdgeCoreProfiles._load_flux_grid`.
    """

    times: np.ndarray  # [nt], [s]
    rgrid: np.ndarray  # [nr], [m]
    zgrid: np.ndarray  # [nz], [m]
    psirz: np.ndarray  # [nt, nz, nr], [Wb/rad]
    psi_axis: np.ndarray  # [nt], [Wb/rad]
    psi_bdry: np.ndarray  # [nt], [Wb/rad]
    fpol: np.ndarray | None  # [nt, npsi], F = R*B_phi [T*m]
    bcentr: np.ndarray | None  # [nt], vacuum field at rcentr [T]
    rcentr: np.ndarray | None  # [nt], reference major radius [m]
    _splines: dict = field(default_factory=dict, init=False, repr=False)

    def _time_index(self, t: float) -> int:
        """Return the index of the EFIT time slice nearest to `t`."""
        return int(np.argmin(np.abs(self.times - t)))

    def _spline(self, it: int) -> RectBivariateSpline:
        """
        Return the cached bicubic spline of psi(Z, R) for time index `it`.

        Splines are cached per time index because the profile reduction
        evaluates the same EFIT slice for every diagnostic chord.
        """
        if it not in self._splines:
            self._splines[it] = RectBivariateSpline(
                self.zgrid, self.rgrid, self.psirz[it], kx=3, ky=3
            )
        return self._splines[it]

    def psinorm(self, r, z, t: float) -> np.ndarray:
        r"""
        Normalized poloidal flux $\psi_N = (\psi - \psi_a)/(\psi_b - \psi_a)$.

        Bivariate-spline interpolation of psi on the regular (R, Z) grid at the
        EFIT slice nearest `t`. $\psi_N$ is sign-robust (independent of the
        COCOS sign of psi) and equals 0 on the magnetic axis and 1 on the LCFS.

        Parameters
        ----------
        r : array_like
            Major radius (or radii) [m]; broadcast against `z`.
        z : array_like
            Vertical position(s) [m]; broadcast against `r`.
        t : float
            Time [s]; the nearest EFIT slice is used.

        Returns
        -------
        np.ndarray
            Normalized poloidal flux at each (R, Z) point.
        """
        r = np.atleast_1d(np.asarray(r, dtype=float)).ravel()
        z = np.atleast_1d(np.asarray(z, dtype=float)).ravel()
        if r.size == 1 and z.size > 1:
            r = np.full(z.shape, r.item())
        elif z.size == 1 and r.size > 1:
            z = np.full(r.shape, z.item())

        psi = self._spline(self._time_index(t)).ev(z, r)

        psi_a = float(np.interp(t, self.times, self.psi_axis))
        psi_b = float(np.interp(t, self.times, self.psi_bdry))
        denom = psi_b - psi_a
        if not np.isfinite(denom) or denom == 0:
            return np.full(r.shape, np.nan)
        return (psi - psi_a) / denom

    def rz2rho(self, r, z, t: float) -> np.ndarray:
        r"""
        Map (R, Z, t) to $\rho = \sqrt{\psi_N}$, clipped at 0.

        Clipping the (unphysical) negative $\psi_N$ excursions to zero before
        the square root is the standard C-Mod convention. Returns NaN where
        $\psi_N$ is undefined.
        """
        psin = self.psinorm(r, z, t)
        rho = np.sqrt(np.clip(psin, 0.0, None))
        rho[~np.isfinite(psin)] = np.nan
        return rho

    def rz2bt(self, r, z, t: float) -> np.ndarray:
        r"""
        Toroidal field magnitude $|B_T| = |F(\psi_N)| / R$ at (R, Z).

        Inside the LCFS, $F = R B_\phi$ is interpolated from `fpol` on its
        normalized flux grid; outside the LCFS (and everywhere if `fpol` is
        missing) the vacuum value $F = B_0 R_0$ from `bcentr`/`rcentr` is used
        ($B_T \sim 1/R$). Returns positive magnitudes.

        References
        -------
        - the Grad-Shafranov identity $B_\phi = F(\psi)/R$, with `fpol` read
        from the EFIT g-file; the vacuum-outside-LCFS blending and the
        `bcentr`/`rcentr` fallback are specific to this implementation
        """
        r = np.atleast_1d(np.asarray(r, dtype=float)).ravel()
        z = np.atleast_1d(np.asarray(z, dtype=float)).ravel()
        it = self._time_index(t)

        f_vac = None
        if self.bcentr is not None and self.rcentr is not None:
            f_vac = abs(
                float(np.interp(t, self.times, self.bcentr))
                * float(np.interp(t, self.times, self.rcentr))
            )

        if self.fpol is not None:
            fpol_slice = np.abs(self.fpol[it])
            psin_grid = np.linspace(0.0, 1.0, len(fpol_slice))
            psin = self.psinorm(r, z, t)
            f = np.interp(np.clip(psin, 0.0, 1.0), psin_grid, fpol_slice)
            outside = ~np.isfinite(psin) | (psin > 1.0)
            f_edge = f_vac if f_vac is not None else fpol_slice[-1]
            f = np.where(outside, f_edge, f)
        elif f_vac is not None:
            f = np.full(r.shape, f_vac)
        else:
            return np.full(r.shape, np.nan)

        return f / r


@dataclass
class _TsProfiles:
    """
    Combined core+edge Thomson profiles. te in eV, ne in m^-3, R/z in m.
    """

    time: np.ndarray  # [nt], [s]
    r_major: float  # vertical-system major radius [m]
    z: np.ndarray  # [nchord], [m]
    te: np.ndarray  # [nchord, nt], [eV]
    ne: np.ndarray  # [nchord, nt], [m^-3]
    te_err: np.ndarray  # [nchord, nt], [eV]
    ne_err: np.ndarray  # [nchord, nt], [m^-3]


class CmodEdgeCoreProfiles:
    """
    A helper class for the C-Mod edge/core profile variables (Maris 2025).

    Only `compute_all`, `region_columns` and `all_columns` are public; every
    other member is an internal helper. The `@physics_method` wrappers live in
    `CmodPhysicsMethods`.
    """

    @staticmethod
    def region_columns(region: str) -> list[str]:
        """
        Return the 13 column names produced for one rho band.

        Parameters
        ----------
        region : str
            The region label, one of the `RHO_BANDS` keys ("edge" or "core").

        Returns
        -------
        list of str
            The 13 column names for the region.
        """
        return [
            f"n_{region}",
            f"t_{region}",
            f"p_{region}",
            f"t_{region}_ece",
            f"n_{region}_fit",
            f"t_{region}_fit",
            f"p_{region}_fit",
            f"nu_star_{region}",
            f"rho_star_{region}",
            f"beta_t_{region}",
            f"nu_star_{region}_fit",
            f"rho_star_{region}_fit",
            f"beta_t_{region}_fit",
        ]

    @staticmethod
    def all_columns() -> list[str]:
        """
        Return every column produced by `compute_all` (27 in total).

        Returns
        -------
        list of str
            The 13 edge columns, the 13 core columns, and `lsvm_d`.
        """
        cols = []
        for region in RHO_BANDS:
            cols.extend(CmodEdgeCoreProfiles.region_columns(region))
        return cols + ["lsvm_d"]

    @staticmethod
    def compute_all(params: PhysicsMethodParams) -> dict:
        """
        Compute every edge/core profile column on `params.times`.

        Reads the EFIT flux grid, the geometry, and the Thomson and ECE
        profiles once, then reduces them over both rho bands. On unrecoverable
        read failures (no EFIT flux grid, no Thomson data) this logs a warning
        and returns all-NaN columns rather than raising, so that the seven
        physics methods sharing this (cached) computation degrade together
        without re-triggering the expensive reads.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the data connection, shot id and more.

        Returns
        -------
        dict
            A dictionary with all 27 edge/core columns (see `all_columns`),
            each an array of `len(params.times)`.
        """
        times = params.times
        empty = {
            col: np.full(len(times), np.nan)
            for col in CmodEdgeCoreProfiles.all_columns()
        }

        # EFIT flux grid (rho / B_T mapping); without it nothing can be done.
        try:
            grid, efit_tree = CmodEdgeCoreProfiles._load_flux_grid_preferred(params)
        except (mdsExceptions.MdsException, ValueError, CalculationError) as e:
            params.logger.warning(f"profiles: EFIT flux grid load failed: {repr(e)}")
            return empty
        params.logger.debug(f"profiles: using EFIT tree '{efit_tree}'")

        # Geometry from the SAME tree as the flux grid; I_p from the standard
        # (tree-independent) reader; B_T at the magnetic axis from the
        # equilibrium itself.
        try:
            a_efit, r0_efit, kappa_efit, t_efit = CmodEdgeCoreProfiles._load_geometry(
                params, efit_tree
            )
        except mdsExceptions.MdsException as e:
            params.logger.warning(
                f"profiles: geometry read from '{efit_tree}' failed: {repr(e)}"
            )
            return empty
        a_minor = interp1(t_efit, a_efit, times)
        r_0 = interp1(t_efit, r0_efit, times)
        kappa = interp1(t_efit, kappa_efit, times)

        # pylint: disable-next=import-outside-toplevel,cyclic-import
        from disruption_py.machine.cmod.physics import CmodPhysicsMethods

        try:
            ip = np.asarray(
                CmodPhysicsMethods.get_ip_parameters(params=params).get("ip"),
                dtype=float,
            )  # [A]
        except mdsExceptions.MdsException as e:
            params.logger.warning(f"profiles: I_p read failed: {repr(e)}")
            ip = np.full(len(times), np.nan)
        bt = interp1(
            t_efit, CmodEdgeCoreProfiles._bt_on_axis(grid, r0_efit, t_efit), times
        )  # [T]

        # Thomson core+edge ne/Te (full profile; the rho bands select the
        # regions below), mapped to rho once.
        try:
            ts = CmodEdgeCoreProfiles._read_ts_profiles(params)
        except (mdsExceptions.MdsException, MismatchCalculationError) as e:
            params.logger.warning(f"profiles: Thomson read failed: {repr(e)}")
            return empty
        ts_rho = CmodEdgeCoreProfiles._map_ts_rho(grid, ts)
        ts_valid = ts.ne > 0  # drop empty/garbage chords

        # GPC2 ECE Te (best-effort cross-check; absent before ~2000).
        try:
            ece_te, ece_time, ece_r = CmodEdgeCoreProfiles._read_ece_te(params)
            ece_rho = CmodEdgeCoreProfiles._map_ece_rho(grid, ece_r, ece_time)
        except (mdsExceptions.MdsException, ValueError, IndexError) as e:
            params.logger.warning(f"profiles: ECE read failed: {repr(e)}")
            ece_te, ece_time, ece_rho = None, None, None

        result = {}
        for region, (lo, hi) in RHO_BANDS.items():
            result.update(
                CmodEdgeCoreProfiles._compute_region(
                    params=params,
                    region=region,
                    lo=lo,
                    hi=hi,
                    ts=ts,
                    ts_rho=ts_rho,
                    ts_valid=ts_valid,
                    ece_te=ece_te,
                    ece_time=ece_time,
                    ece_rho=ece_rho,
                    a_minor=a_minor,
                    r_0=r_0,
                    kappa=kappa,
                    ip=ip,
                    bt=bt,
                )
            )
        result["lsvm_d"] = CmodEdgeCoreProfiles._lsvm_d(
            result["nu_star_edge"], result["beta_t_edge"]
        )
        return result

    # ------------------------------------------------------------------
    # EFIT flux grid and geometry loaders
    # ------------------------------------------------------------------
    @staticmethod
    def _load_flux_grid(params: PhysicsMethodParams, tree: str) -> _EfitFluxGrid:
        """
        Load the EFIT flux grid from one EFIT tree.

        Reads psirz (with its three dim arrays), the a-file scalars
        simagx/sibdry (interpolated onto the psirz timebase), and best-effort
        fpol/bcentr/rcentr for the |B_T| mapping.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the data connection, shot id and more.
        tree : str
            Literal EFIT tree name (e.g. "efit21", "analysis").

        Returns
        -------
        _EfitFluxGrid
            The eagerly-loaded flux grid.

        References
        -------
        - original source: node-set pattern of [xTomo `read_efit_psi`](https://
        github.com/chandrarn/xTomo/blob/main/src/xtomo/xtomo_mds.py#L270), by
        R Chandra (no license declared); I/O layer rewritten for the
        DisruptionPy connection, fpol/bcentr/rcentr added
        """
        psirz_raw, d_0, d_1, d_2 = params.get_data_with_dims(
            f"{_EFIT_GEQDSK}:psirz", tree_name=tree, dim_nums=[0, 1, 2]
        )  # [Wb/rad], (dims)
        psirz, rgrid, zgrid, times = CmodEdgeCoreProfiles._orient_psirz(
            np.asarray(psirz_raw),
            np.asarray(d_0, dtype=float),
            np.asarray(d_1, dtype=float),
            np.asarray(d_2, dtype=float),
        )

        psi_axis, t_axis = params.get_data_with_dims(
            f"{_EFIT_AEQDSK}:simagx", tree_name=tree, dim_nums=[0]
        )  # [Wb/rad], [s]
        psi_bdry, t_bdry = params.get_data_with_dims(
            f"{_EFIT_AEQDSK}:sibdry", tree_name=tree, dim_nums=[0]
        )  # [Wb/rad], [s]
        psi_axis = np.interp(
            times, np.asarray(t_axis, dtype=float), np.asarray(psi_axis, dtype=float)
        )
        psi_bdry = np.interp(
            times, np.asarray(t_bdry, dtype=float), np.asarray(psi_bdry, dtype=float)
        )

        fpol = CmodEdgeCoreProfiles._try_load_time_profile(
            params, f"{_EFIT_GEQDSK}:fpol", tree, times
        )  # [T*m]
        bcentr = CmodEdgeCoreProfiles._try_load_scalar(
            params, f"{_EFIT_AEQDSK}:bcentr", tree, times
        )  # [T]
        rcentr = CmodEdgeCoreProfiles._try_load_scalar(
            params, f"{_EFIT_AEQDSK}:rcentr", tree, times
        )  # [m]

        return _EfitFluxGrid(
            times=times,
            rgrid=rgrid,
            zgrid=zgrid,
            psirz=psirz,
            psi_axis=psi_axis,
            psi_bdry=psi_bdry,
            fpol=fpol,
            bcentr=bcentr,
            rcentr=rcentr,
        )

    @staticmethod
    def _load_flux_grid_preferred(
        params: PhysicsMethodParams, trees: tuple = EFIT_TREE_PREFERENCE
    ):
        """
        Load the EFIT flux grid from the first available tree in `trees`.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the data connection, shot id and more.
        trees : tuple of str, optional
            Tree names to try in order, by default `EFIT_TREE_PREFERENCE`.

        Returns
        -------
        tuple
            `(grid, tree_used)`; raises the last error if no tree works.
        """
        last_error = CalculationError("no EFIT tree available")
        for tree in trees:
            try:
                return CmodEdgeCoreProfiles._load_flux_grid(params, tree), tree
            except (mdsExceptions.MdsException, ValueError) as error:
                params.logger.debug(
                    f"profiles: EFIT tree '{tree}' unavailable: {repr(error)}"
                )
                last_error = error
        raise last_error

    @staticmethod
    def _load_geometry(params: PhysicsMethodParams, tree: str):
        """
        Read minor radius, axis major radius and elongation from one EFIT tree.

        Read here (rather than via `CmodEfitMethods.get_efit_parameters`, which
        is pinned to the `_efit_tree` nickname) so the geometry matches the
        same reconstruction as the rho/B_T flux grid.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the data connection, shot id and more.
        tree : str
            Literal EFIT tree name, matching the flux grid.

        Returns
        -------
        tuple
            `(a_minor [m], r_0 [m], kappa, efit_time [s])` on the tree's own
            timebase.

        References
        -------
        - referenced source: the `aout/100` / `rmagx/100` node reads of
        `CmodPhysicsMethods.get_peaking_factors` and `get_te_profile_params_ece`
        in `disruption_py/machine/cmod/physics.py`
        """
        a_minor, t = params.get_data_with_dims(
            f"{_EFIT_AEQDSK}:aout/100", tree_name=tree, dim_nums=[0]
        )  # [m], [s]
        r_0 = params.get_data(f"{_EFIT_AEQDSK}:rmagx/100", tree_name=tree)  # [m]
        kappa = params.get_data(f"{_EFIT_AEQDSK}:eout", tree_name=tree)  # []
        return (
            np.asarray(a_minor, dtype=float),
            np.asarray(r_0, dtype=float),
            np.asarray(kappa, dtype=float),
            np.asarray(t, dtype=float),
        )

    @staticmethod
    def _orient_psirz(psirz, d_0, d_1, d_2):
        """
        Orient psirz to `[nt, nz, nr]` from arbitrary MDSplus axis ordering.

        EFIT/MDSplus axis ordering of psirz is not guaranteed, so each axis is
        detected from the three dim-of arrays instead of assuming a fixed
        order: the R and Z flux grids are square (equal length, typically 65 or
        129) so the time dim is the odd one out; among the spatial grids, R is
        strictly positive while Z straddles 0. In the fully ambiguous square
        case, axis 0 is assumed to be time (the C-Mod convention).

        Parameters
        ----------
        psirz : np.ndarray
            3D poloidal flux array in any axis order.
        d_0, d_1, d_2 : np.ndarray
            The dim-of arrays of psirz axes 0, 1, 2.

        Returns
        -------
        tuple
            `(psirz [nt, nz, nr], rgrid [nr], zgrid [nz], times [nt])`.
        """
        if psirz.ndim != 3:
            raise ValueError(f"psirz expected 3D, got shape {psirz.shape}")
        dims = [np.asarray(d_0, float), np.asarray(d_1, float), np.asarray(d_2, float)]
        lengths = [len(d) for d in dims]
        shape = psirz.shape

        # Choose the time dim: two grids share a length (square R/Z), so time
        # is the unique-length dim.
        counts = {n: lengths.count(n) for n in lengths}
        if any(c == 1 for c in counts.values()) and not all(
            c == 1 for c in counts.values()
        ):
            t_dim = next(i for i, n in enumerate(lengths) if counts[n] == 1)
        else:
            # Fall back: time is the dim that does not look like a spatial
            # grid (C-Mod R and |Z| are well below 2 m).
            def looks_spatial(v):
                return np.nanmax(np.abs(v)) < 2.0

            cand = [i for i in range(3) if not looks_spatial(dims[i])]
            t_dim = cand[0] if cand else int(np.argmax(lengths))

        spatial = [i for i in range(3) if i != t_dim]
        # R is the all-positive grid; Z straddles zero.
        if np.nanmin(dims[spatial[0]]) < 0 <= np.nanmax(dims[spatial[0]]):
            z_dim, r_dim = spatial[0], spatial[1]
        else:
            r_dim, z_dim = spatial[0], spatial[1]

        times, rgrid, zgrid = dims[t_dim], dims[r_dim], dims[z_dim]

        # Map each data axis to its dim by matching lengths; when R and Z share
        # a length this is ambiguous and axis 0 is assumed to be time.
        nt, nz, nr = len(times), len(zgrid), len(rgrid)
        data_t = next(
            (ax for ax in range(3) if shape[ax] == nt and nt not in (nz, nr)), None
        )
        if data_t is None:
            data_t = 0
        psirz = np.moveaxis(psirz, data_t, 0)
        # Remaining two axes -> (Z, R); transpose if swapped.
        if psirz.shape[1] == nr and psirz.shape[2] == nz and nr != nz:
            psirz = np.transpose(psirz, (0, 2, 1))
        return np.ascontiguousarray(psirz, dtype=float), rgrid, zgrid, times

    @staticmethod
    def _try_load_time_profile(params: PhysicsMethodParams, path, tree, times):
        """
        Load a [time, profile] node oriented to `[nt, npsi]`; None on failure.

        Each profile column is interpolated onto `times` (the psirz timebase)
        if the node has its own timebase.
        """
        try:
            data, t = params.get_data_with_dims(path, tree_name=tree, dim_nums=[1])
        except mdsExceptions.MdsException:
            return None
        data = np.asarray(data, dtype=float)
        t = np.asarray(t, dtype=float)
        if data.ndim != 2:
            return None
        if data.shape[0] != len(t) and data.shape[1] == len(t):
            data = data.T
        if data.shape[0] != len(times):
            out = np.empty((len(times), data.shape[1]))
            for j in range(data.shape[1]):
                out[:, j] = np.interp(times, t, data[:, j])
            return out
        return data

    @staticmethod
    def _try_load_scalar(params: PhysicsMethodParams, path, tree, times):
        """Load a scalar-vs-time node interpolated onto `times`; None on failure."""
        try:
            val, t = params.get_data_with_dims(path, tree_name=tree, dim_nums=[0])
        except mdsExceptions.MdsException:
            return None
        return np.interp(
            times, np.asarray(t, dtype=float), np.asarray(val, dtype=float)
        )

    @staticmethod
    def _bt_on_axis(grid: _EfitFluxGrid, r0_efit, t_efit):
        """
        Return |B_T| at the magnetic axis on the EFIT timebase.

        At the axis psi_N = 0 by definition, so F = fpol[:, 0] exactly and
        |B_T| = |F|/R_0, with the vacuum `bcentr*rcentr` fallback when fpol is
        unavailable. This is equivalent to evaluating `rz2bt` at (R_0, 0)
        without the spline evaluation.

        Parameters
        ----------
        grid : _EfitFluxGrid
            The loaded flux grid.
        r0_efit : np.ndarray
            Magnetic-axis major radius on the EFIT timebase [m].
        t_efit : np.ndarray
            The EFIT timebase [s].

        Returns
        -------
        np.ndarray
            |B_T| at the axis for each EFIT time [T].
        """
        bt = np.full(len(t_efit), np.nan)
        for i, t in enumerate(t_efit):
            r0_val = r0_efit[i]
            if not np.isfinite(r0_val) or r0_val <= 0:
                continue
            it = int(np.argmin(np.abs(grid.times - t)))
            if grid.fpol is not None:
                bt[i] = abs(grid.fpol[it, 0]) / r0_val
            elif grid.bcentr is not None and grid.rcentr is not None:
                bt[i] = abs(grid.bcentr[it] * grid.rcentr[it]) / r0_val
        return bt

    # ------------------------------------------------------------------
    # Diagnostic readers
    # ------------------------------------------------------------------
    @staticmethod
    def _read_ts_profiles(params: PhysicsMethodParams) -> _TsProfiles:
        r"""
        Read the C-Mod core + edge Thomson ne/Te profiles.

        Core: `\.yag_new.results.profiles:{te_rz [keV], ne_rz [m^-3],
        z_sorted [m], te_err, ne_err}`. Edge (shot > 1e9 only): `\ts_te` [eV],
        `\ts_ne` [m^-3], `\fiber_z` [m]. Both share the vertical-system major
        radius `\.yag.results.param:r`. Raises `MismatchCalculationError` if
        the edge chord positions and profiles disagree in length.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the data connection, shot id and more.

        Returns
        -------
        _TsProfiles
            The combined profiles, with temperatures kept in eV.

        References
        -------
        - referenced source: the Thomson read block of
        `CmodPhysicsMethods.get_peaking_factors` in `disruption_py/machine/
        cmod/physics.py` (keV -> eV conversion, core/edge concatenation, edge
        chord-count consistency check); restructured as a standalone reader
        with best-effort error nodes
        """
        te_core, ts_time = params.get_data_with_dims(
            f"{_TS_CORE_NODE}:te_rz", tree_name="electrons"
        )  # [keV], [s]
        te_core = np.asarray(te_core, dtype=float) * 1e3  # [keV] -> [eV]
        ne_core = np.asarray(
            params.get_data(f"{_TS_CORE_NODE}:ne_rz", tree_name="electrons"),
            dtype=float,
        )  # [m^-3]
        z_core = np.asarray(
            params.get_data(f"{_TS_CORE_NODE}:z_sorted", tree_name="electrons"),
            dtype=float,
        )  # [m]
        te_core_err = (
            np.asarray(
                params.get_data(f"{_TS_CORE_NODE}:te_err", tree_name="electrons"),
                dtype=float,
            )
            * 1e3
        )  # [keV] -> [eV]
        ne_core_err = np.asarray(
            params.get_data(f"{_TS_CORE_NODE}:ne_err", tree_name="electrons"),
            dtype=float,
        )  # [m^-3]
        r_major = float(
            np.ravel(params.get_data(_TS_CORE_R, tree_name="electrons"))[0]
        )  # [m]
        ts_time = np.asarray(ts_time, dtype=float)

        try:
            te_edge = np.asarray(
                params.get_data(r"\ts_te", tree_name="electrons"), dtype=float
            )  # [eV]
            ne_edge = np.asarray(
                params.get_data(r"\ts_ne", tree_name="electrons"), dtype=float
            )  # [m^-3]
            z_edge = np.asarray(
                params.get_data(r"\fiber_z", tree_name="electrons"), dtype=float
            )  # [m]
            if len(z_edge) != te_edge.shape[0]:
                raise MismatchCalculationError(
                    f"len(z_edge) = {len(z_edge)} vs. "
                    f"te_edge.shape[0] = {te_edge.shape[0]}"
                )
            te_edge_err = CmodEdgeCoreProfiles._try_get(
                params, r"\ts_te_err", te_edge.shape
            )
            ne_edge_err = CmodEdgeCoreProfiles._try_get(
                params, r"\ts_ne_err", ne_edge.shape
            )
            z = np.concatenate([z_core, z_edge])
            te = np.concatenate([te_core, te_edge], axis=0)
            ne = np.concatenate([ne_core, ne_edge], axis=0)
            te_err = np.concatenate([te_core_err, te_edge_err], axis=0)
            ne_err = np.concatenate([ne_core_err, ne_edge_err], axis=0)
        except mdsExceptions.MdsException:
            params.logger.debug("profiles: no edge Thomson system; using core only")
            z, te, ne = z_core, te_core, ne_core
            te_err, ne_err = te_core_err, ne_core_err

        return _TsProfiles(
            time=ts_time,
            r_major=r_major,
            z=z,
            te=te,
            ne=ne,
            te_err=te_err,
            ne_err=ne_err,
        )

    @staticmethod
    def _try_get(params: PhysicsMethodParams, path, shape, fill=np.nan):
        """Best-effort node read from the electrons tree; `fill` on failure."""
        try:
            return np.asarray(params.get_data(path, tree_name="electrons"), dtype=float)
        except mdsExceptions.MdsException:
            return np.full(shape, fill)

    @staticmethod
    def _read_ece_te(params: PhysicsMethodParams):
        """
        Read the C-Mod GPC2 ECE Te(R, t) profiles.

        WARNING: unlike `CmodPhysicsMethods._get_te_profile_params_ece`, this
        reader applies NONE of the ECE quality filters (B_T threshold, LH
        heating, density cutoff, harmonic overlap, non-thermal rising tail,
        outlier rejection). The resulting `t_*_ece` columns are a rough
        cross-check of the Thomson band averages only and must not be used
        quantitatively.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the data connection, shot id and more.

        Returns
        -------
        tuple
            `(te [eV] [nchan, nt], te_time [s], r [m] [nchan, nt])`, with each
            channel's radius interpolated onto the Te timebase.

        References
        -------
        - referenced source: the GPC2 node reads and the radii-interpolation
        loop of `CmodPhysicsMethods.get_te_profile_params_ece` in
        `disruption_py/machine/cmod/physics.py`
        """
        te, te_time = params.get_data_with_dims(
            f"{_ECE_NODE}:gpc2_te", tree_name="electrons"
        )  # [keV], [s]
        te = np.asarray(te, dtype=float) * 1e3  # [keV] -> [eV]
        te_time = np.asarray(te_time, dtype=float)
        radii, radii_time = params.get_data_with_dims(
            f"{_ECE_NODE}:radii", tree_name="electrons"
        )  # [m], [s]
        radii = np.asarray(radii, dtype=float)
        radii_time = np.asarray(radii_time, dtype=float)

        # Interpolate each channel's radius onto the Te timebase (the radii
        # are sampled more slowly).
        r_on_te = np.full((radii.shape[0], len(te_time)), np.nan)
        for i in range(radii.shape[0]):
            if radii.shape[1] > 1:
                r_on_te[i, :] = np.interp(te_time, radii_time, radii[i, :])
            else:
                r_on_te[i, :] = radii[i, 0]
        return te, te_time, r_on_te

    @staticmethod
    def _map_ts_rho(grid: _EfitFluxGrid, ts: _TsProfiles) -> np.ndarray:
        """
        Map every Thomson chord and time slice to rho.

        Returns a `[nchord, nt]` array; NaN for times outside the EFIT range.
        """
        out = np.full((len(ts.z), len(ts.time)), np.nan)
        t_min, t_max = grid.times.min(), grid.times.max()
        r = np.full(ts.z.shape, ts.r_major)
        for i, t in enumerate(ts.time):
            if t < t_min or t > t_max:
                continue
            out[:, i] = grid.rz2rho(r, ts.z, float(t))
        return out

    @staticmethod
    def _map_ece_rho(grid: _EfitFluxGrid, r_on_te, te_time) -> np.ndarray:
        """
        Map every ECE channel radius (at Z = 0) and time slice to rho.

        Returns a `[nchan, nt]` array; NaN for times outside the EFIT range.
        """
        out = np.full(r_on_te.shape, np.nan)
        t_min, t_max = grid.times.min(), grid.times.max()
        for i, t in enumerate(te_time):
            if t < t_min or t > t_max:
                continue
            out[:, i] = grid.rz2rho(r_on_te[:, i], 0.0, float(t))
        return out

    # ------------------------------------------------------------------
    # Profile reducers (pure numeric)
    # ------------------------------------------------------------------
    @staticmethod
    def _band_average(values, rho, lo, hi) -> float:
        """Mean of the finite `values` with rho in [lo, hi]; NaN if empty."""
        mask = np.isfinite(values) & np.isfinite(rho) & (rho >= lo) & (rho <= hi)
        if not np.any(mask):
            return np.nan
        return float(np.mean(values[mask]))

    @staticmethod
    def _band_average_profiles(values, valid, rho, lo, hi) -> np.ndarray:
        """
        Band-average a `[nchord, nt]` profile over rho in [lo, hi] per slice.

        Method A of the port: the raw band average, matching the paper's
        C-Mod treatment (Maris 2025 section 2.2).

        Parameters
        ----------
        values : np.ndarray
            Profile values `[nchord, nt]`.
        valid : np.ndarray
            Boolean chord-validity mask `[nchord, nt]`.
        rho : np.ndarray
            rho of each chord and slice `[nchord, nt]`.
        lo, hi : float
            The rho band.

        Returns
        -------
        np.ndarray
            The per-slice band average `[nt]`.
        """
        n_slices = values.shape[1]
        out = np.full(n_slices, np.nan)
        for i in range(n_slices):
            vals = np.where(valid[:, i], values[:, i], np.nan)
            out[i] = CmodEdgeCoreProfiles._band_average(vals, rho[:, i], lo, hi)
        return out

    @staticmethod
    def _osborne_mtanh(x, c_0, c_1, c_2, c_3, c_4, c_5, c_6):
        """
        Osborne modified-tanh profile with cubic inboard term and flat SOL.

        `c_0` pedestal center, `c_1` pedestal full width, `c_2` pedestal top,
        `c_3` pedestal bottom, `c_4`/`c_5`/`c_6` inboard linear/quadratic/cubic
        terms.

        References
        -------
        - the modified-tanh pedestal form of T Osborne, as implemented for
        Alcator C-Mod in the JW Hughes IDL analysis scripts
        """
        z = 2.0 * (c_0 - x) / c_1
        p_1 = 1.0 + c_4 * z + c_5 * z**2 + c_6 * z**3
        e_1 = np.exp(z)
        e_2 = np.exp(-z)
        return 0.5 * (c_2 + c_3 + (c_2 - c_3) * (p_1 * e_1 - e_2) / (e_1 + e_2))

    @staticmethod
    def _fit_band_average(psin_pts, y_pts, lo2, hi2) -> float:
        """
        Fit y(psi_N) with the Osborne mtanh; mean over psi_N in [lo2, hi2].

        Method B of the port. Returns NaN if there are too few points, the
        retained points do not reach the band (no extrapolation), or the fit
        does not converge.

        Parameters
        ----------
        psin_pts : np.ndarray
            Normalized poloidal flux of each data point (NaN to exclude).
        y_pts : np.ndarray
            Data values at each point.
        lo2, hi2 : float
            The band in psi_N (i.e. rho^2) over which the fit is averaged.

        Returns
        -------
        float
            The band-averaged fit value, or NaN.

        Notes
        -------
        The pedestal-centre bound `c_0` in [0.85, 1.1] is the usual C-Mod
        range in psi_N. The remaining bounds are scaled by the pedestal top so
        that one code path serves both ne and Te.
        """
        mask = (
            np.isfinite(psin_pts)
            & np.isfinite(y_pts)
            & (y_pts > 0)
            & (psin_pts >= 0)
            & (psin_pts <= 1.3)
        )
        x = psin_pts[mask]
        y = y_pts[mask]
        if x.size < 7:
            return np.nan
        # No extrapolation: refuse to evaluate the fit over a band the data
        # does not even reach (the raw band average correctly reports NaN
        # there, and the fit column must not silently disagree).
        if x.min() > hi2 or x.max() < lo2:
            return np.nan
        order = np.argsort(x)
        x, y = x[order], y[order]
        ped_top = float(np.nanmax(y))
        guess = [0.98, 0.05, ped_top, 0.0, 0.0, 0.0, 0.0]
        lower = [0.85, 0.005, 0.0, -0.05 * ped_top, -np.inf, -np.inf, -np.inf]
        upper = [1.1, 0.3, 5.0 * ped_top + 1.0, 0.5 * ped_top + 1.0]
        upper += [np.inf, np.inf, np.inf]
        try:
            popt, _ = curve_fit(
                CmodEdgeCoreProfiles._osborne_mtanh,
                x,
                y,
                p0=guess,
                bounds=(lower, upper),
                maxfev=4000,
            )
        except (RuntimeError, ValueError):
            return np.nan
        xs = np.linspace(lo2, hi2, 21)
        vals = CmodEdgeCoreProfiles._osborne_mtanh(xs, *popt)
        vals = vals[np.isfinite(vals) & (vals > 0)]
        return float(np.mean(vals)) if vals.size else np.nan

    @staticmethod
    def _fit_profiles(values, valid, ts_rho, lo, hi) -> np.ndarray:
        """
        Osborne-fit band average of a `[nchord, nt]` profile, per time slice.

        Parameters
        ----------
        values : np.ndarray
            Profile values `[nchord, nt]`.
        valid : np.ndarray
            Boolean chord-validity mask `[nchord, nt]`.
        ts_rho : np.ndarray
            rho of each chord and slice `[nchord, nt]`.
        lo, hi : float
            The rho band; the fit is averaged over psi_N in [lo^2, hi^2].

        Returns
        -------
        np.ndarray
            The per-slice fitted band average `[nt]`.
        """
        n_slices = values.shape[1]
        out = np.full(n_slices, np.nan)
        psin = ts_rho**2
        lo2, hi2 = lo**2, hi**2
        for i in range(n_slices):
            x = np.where(valid[:, i], psin[:, i], np.nan)
            out[i] = CmodEdgeCoreProfiles._fit_band_average(x, values[:, i], lo2, hi2)
        return out

    # ------------------------------------------------------------------
    # Physics
    # ------------------------------------------------------------------
    @staticmethod
    def _coulomb_log(n_m3, t_ev):
        r"""
        Coulomb logarithm, NRL electron formula (valid for $T_e > 10$ eV).

        $\ln\Lambda = 24 - \ln(\sqrt{n_e[\mathrm{cm}^{-3}]}/T_e[\mathrm{eV}])$.
        Maris 2025 does not state the $\ln\Lambda$ expression used, so the NRL
        Plasma Formulary form is adopted; this is the one place where a
        bit-level comparison against the published values may differ (the
        boundary of eq. (9) carries $(\ln\Lambda)^{0.7}$ with < 10% spread
        across the database, so the effect is a few percent at most).

        References
        -------
        - NRL Plasma Formulary (electron-electron collisions)
        """
        n_cm3 = np.asarray(n_m3, dtype=float) * 1e-6
        t_ev = np.asarray(t_ev, dtype=float)
        with np.errstate(invalid="ignore", divide="ignore"):
            return 24.0 - np.log(np.sqrt(n_cm3) / t_ev)

    @staticmethod
    def _dimensionless(n_m3, t_ev, a_minor, r_0, kappa, ip, bt) -> dict:
        r"""
        Dimensionless band variables (Maris 2025 eqs. 2-4) and total pressure.

        With T supplied in eV and converted to J internally:

        $$
        q_\mathrm{cyl} = \frac{2\pi}{\mu_0}\frac{B_T a^2}{I_p R_0}
        \frac{1+\kappa^2}{2}, \qquad
        \nu^* = \frac{e^4\ln\Lambda}{2\pi\varepsilon_0^2}\frac{n}{T^2}
        \frac{q_\mathrm{cyl}R_0}{\epsilon^{3/2}} \;(\mathrm{eq.}\,2)
        $$
        $$
        \rho^* = \frac{\sqrt{m_D T}}{e B_T a} \;(\mathrm{eq.}\,3), \qquad
        \beta_T = \frac{4\mu_0 n T}{B_T^2} \;(\mathrm{eq.}\,4), \qquad
        p = 2 n T
        $$

        $\beta_T$ is a fraction (not a percent); $p$ is the total pressure
        assuming $T_i = T_e$ and $n_i = n_e$; the ion mass is deuterium.
        Slices with $|I_p| <$ `MIN_IP` are returned as NaN.

        Parameters
        ----------
        n_m3 : array_like
            Band electron density [m^-3].
        t_ev : array_like
            Band electron temperature [eV].
        a_minor, r_0 : array_like
            Minor radius and magnetic-axis major radius [m].
        kappa : array_like
            Elongation.
        ip : array_like
            Plasma current [A]; the magnitude is used.
        bt : array_like
            Toroidal field at the axis [T]; the magnitude is used.

        Returns
        -------
        dict
            A dictionary with `nu_star`, `rho_star`, `beta_t`, `pressure`
            [Pa], and `q_cyl`.

        References
        -------
        - publication: AD Maris, et al. (2025), Nuclear Fusion 65 016051,
        DOI: [10.1088/1741-4326/ad90f0](https://doi.org/10.1088/1741-4326/
        ad90f0), eqs. (2)-(4)
        """
        n = np.asarray(n_m3, dtype=float)
        t_ev = np.asarray(t_ev, dtype=float)
        t_j = t_ev * const.e  # [eV] -> [J]
        a_minor = np.asarray(a_minor, dtype=float)
        r_0 = np.asarray(r_0, dtype=float)
        kappa = np.asarray(kappa, dtype=float)
        ip = np.abs(np.asarray(ip, dtype=float))
        bt = np.abs(np.asarray(bt, dtype=float))
        ip = np.where(ip > MIN_IP, ip, np.nan)

        eps = a_minor / r_0
        ln_lambda = CmodEdgeCoreProfiles._coulomb_log(n, t_ev)

        with np.errstate(invalid="ignore", divide="ignore"):
            q_cyl = (
                (2.0 * np.pi / const.mu_0)
                * (bt * a_minor**2 / (ip * r_0))
                * ((1.0 + kappa**2) / 2.0)
            )
            nu_star = (
                (const.e**4 * ln_lambda)
                / (2.0 * np.pi * const.epsilon_0**2)
                * (n / t_j**2)
                * (q_cyl * r_0 / eps**1.5)
            )
            rho_star = np.sqrt(M_DEUTERON * t_j) / (const.e * bt * a_minor)
            beta_t = 4.0 * const.mu_0 * n * t_j / bt**2
            pressure = 2.0 * n * t_j

        return {
            "nu_star": nu_star,
            "rho_star": rho_star,
            "beta_t": beta_t,
            "pressure": pressure,
            "q_cyl": q_cyl,
        }

    @staticmethod
    def _lsvm_d(nu_star_edge, beta_t_edge) -> np.ndarray:
        r"""
        L-mode density limit proximity metric (Maris 2025 eq. 7).

        $\mathrm{lsvm}_d = \nu^*_\mathrm{edge} /
        (3.5\,\beta_{T,\mathrm{edge}}^{-0.40})$, with $\beta_{T,\mathrm{edge}}$
        as a fraction. Values above 1 place the discharge on the
        density-limit-unstable side of the boundary.

        References
        -------
        - publication: AD Maris, et al. (2025), Nuclear Fusion 65 016051,
        DOI: [10.1088/1741-4326/ad90f0](https://doi.org/10.1088/1741-4326/
        ad90f0), eq. (7)
        """
        nu_star_edge = np.asarray(nu_star_edge, dtype=float)
        beta_t_edge = np.asarray(beta_t_edge, dtype=float)
        with np.errstate(invalid="ignore", divide="ignore"):
            return nu_star_edge / (LSVM_D_COEFF * beta_t_edge**LSVM_D_EXPONENT)

    # ------------------------------------------------------------------
    # Per-region assembly
    # ------------------------------------------------------------------
    @staticmethod
    def _compute_region(
        params: PhysicsMethodParams,
        region: str,
        lo: float,
        hi: float,
        ts: _TsProfiles,
        ts_rho: np.ndarray,
        ts_valid: np.ndarray,
        ece_te,
        ece_time,
        ece_rho,
        a_minor,
        r_0,
        kappa,
        ip,
        bt,
    ) -> dict:
        """
        Compute the 13 columns for one rho band on `params.times`.

        Reuses the shared flux-grid mapping, Thomson/ECE reads and geometry;
        see `compute_all` for the orchestration.

        Returns
        -------
        dict
            The 13 `<name>_<region>` columns for this band.
        """
        times = params.times
        nan_col = np.full(len(times), np.nan)

        # Method A: raw band average.
        n_a = interp1(
            ts.time,
            CmodEdgeCoreProfiles._band_average_profiles(
                ts.ne, ts_valid, ts_rho, lo, hi
            ),
            times,
        )
        t_a = interp1(
            ts.time,
            CmodEdgeCoreProfiles._band_average_profiles(
                ts.te, ts_valid, ts_rho, lo, hi
            ),
            times,
        )

        # Method B: Osborne modified-tanh fit (best-effort; NaN on failure).
        try:
            n_f = interp1(
                ts.time,
                CmodEdgeCoreProfiles._fit_profiles(ts.ne, ts_valid, ts_rho, lo, hi),
                times,
            )
            t_f = interp1(
                ts.time,
                CmodEdgeCoreProfiles._fit_profiles(ts.te, ts_valid, ts_rho, lo, hi),
                times,
            )
        except (RuntimeError, ValueError) as e:
            params.logger.warning(f"profiles: {region} Osborne fit failed: {repr(e)}")
            n_f, t_f = nan_col.copy(), nan_col.copy()

        # ECE cross-check temperature in this band (unfiltered; see
        # `_read_ece_te`). Best-effort: a malformed ECE profile (e.g. a
        # channel-count mismatch between `:gpc2_te` and `:radii`) must not
        # take down the Thomson-derived columns.
        if ece_te is None:
            t_ece = nan_col.copy()
        else:
            try:
                t_ece = interp1(
                    ece_time,
                    CmodEdgeCoreProfiles._band_average_profiles(
                        ece_te, ece_te > 0, ece_rho, lo, hi
                    ),
                    times,
                )
            except (ValueError, IndexError) as e:
                params.logger.warning(
                    f"profiles: {region} ECE band average failed: {repr(e)}"
                )
                t_ece = nan_col.copy()

        # Dimensionless variables (Maris 2025 eqs. 2-4) + total pressure, for
        # both profile-reduction methods.
        dim = CmodEdgeCoreProfiles._dimensionless(n_a, t_a, a_minor, r_0, kappa, ip, bt)
        dim_f = CmodEdgeCoreProfiles._dimensionless(
            n_f, t_f, a_minor, r_0, kappa, ip, bt
        )
        return {
            f"n_{region}": n_a,
            f"t_{region}": t_a,
            f"p_{region}": dim["pressure"],
            f"t_{region}_ece": t_ece,
            f"n_{region}_fit": n_f,
            f"t_{region}_fit": t_f,
            f"p_{region}_fit": dim_f["pressure"],
            f"nu_star_{region}": dim["nu_star"],
            f"rho_star_{region}": dim["rho_star"],
            f"beta_t_{region}": dim["beta_t"],
            f"nu_star_{region}_fit": dim_f["nu_star"],
            f"rho_star_{region}_fit": dim_f["rho_star"],
            f"beta_t_{region}_fit": dim_f["beta_t"],
        }

    # ------------------------------------------------------------------
    # Self-validation
    # ------------------------------------------------------------------
    @staticmethod
    def _validate_rho_mapping(
        params: PhysicsMethodParams,
        grid: _EfitFluxGrid,
        tree: str = "analysis",
        t: float | None = None,
    ) -> dict:
        """
        Sanity-check the rho/|B_T| mapping against the EFIT equilibrium.

        Maps the LCFS boundary points (`rbbbs`/`zbbbs`) and confirms rho ~ 1,
        and the magnetic axis (`rmagx`/`zmagx`) and confirms rho ~ 0. Intended
        for interactive/diagnostic use (e.g. the branch plotting script).

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the data connection, shot id and more.
        grid : _EfitFluxGrid
            The loaded flux grid to validate.
        tree : str, optional
            The EFIT tree holding the boundary nodes, by default "analysis".
        t : float, optional
            Time [s]; by default the middle of the EFIT timebase.

        Returns
        -------
        dict
            A report with `rho_lcfs_mean`, `rho_lcfs_std`, `rho_axis` and
            `bt_axis` (or the corresponding read errors).
        """
        if t is None:
            t = float(grid.times[len(grid.times) // 2])
        it = int(np.argmin(np.abs(grid.times - t)))

        report: dict = {"time": t}
        try:
            rbbbs = np.asarray(
                params.get_data(f"{_EFIT_GEQDSK}:rbbbs", tree_name=tree), dtype=float
            )  # [m]
            zbbbs = np.asarray(
                params.get_data(f"{_EFIT_GEQDSK}:zbbbs", tree_name=tree), dtype=float
            )  # [m]
            if rbbbs.ndim == 2:
                rbbbs = rbbbs[min(it, rbbbs.shape[0] - 1)]
            if zbbbs.ndim == 2:
                zbbbs = zbbbs[min(it, zbbbs.shape[0] - 1)]
            valid = (rbbbs > 0) & np.isfinite(rbbbs) & np.isfinite(zbbbs)
            rho_lcfs = grid.rz2rho(rbbbs[valid], zbbbs[valid], t)
            report["rho_lcfs_mean"] = float(np.nanmean(rho_lcfs))
            report["rho_lcfs_std"] = float(np.nanstd(rho_lcfs))
        except (mdsExceptions.MdsException, ValueError, IndexError) as exc:
            report["lcfs_error"] = repr(exc)

        try:
            rmagx = float(
                np.ravel(params.get_data(f"{_EFIT_AEQDSK}:rmagx/100", tree_name=tree))[
                    it
                ]
            )  # [m]
            zmagx = float(
                np.ravel(params.get_data(f"{_EFIT_AEQDSK}:zmagx/100", tree_name=tree))[
                    it
                ]
            )  # [m]
            report["rho_axis"] = float(grid.rz2rho(rmagx, zmagx, t)[0])
            report["bt_axis"] = float(grid.rz2bt(rmagx, zmagx, t)[0])
        except (mdsExceptions.MdsException, ValueError, IndexError) as exc:
            report["axis_error"] = repr(exc)

        return report
