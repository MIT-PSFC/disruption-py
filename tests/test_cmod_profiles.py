#!/usr/bin/env python3

"""
Unit tests for the C-Mod edge/core profile machinery.

These exercise the pure math and the (R, Z, t) -> rho / |B_T| mapping on a
synthetic analytic equilibrium, so they run without MDSplus access.
"""

# The tests deliberately exercise the internal helpers of the profiles module.
# pylint: disable=protected-access

import numpy as np
import pytest
import scipy.constants as const

from disruption_py.config import config
from disruption_py.machine.cmod.physics import CmodPhysicsMethods
from disruption_py.machine.cmod.profiles import (
    LSVM_D_COEFF,
    LSVM_D_EXPONENT,
    M_DEUTERON,
    CmodEdgeCoreProfiles,
    _EfitFluxGrid,
)
from disruption_py.machine.tokamak import Tokamak

R_AXIS = 0.68
A_MINOR = 0.22


@pytest.fixture(name="grid")
def grid_fixture() -> _EfitFluxGrid:
    """Circular equilibrium psi = (R-R0)^2 + Z^2 so that rho = r / a_minor."""
    rgrid = np.linspace(0.40, 0.95, 65)
    zgrid = np.linspace(-0.40, 0.40, 65)
    rr, zz = np.meshgrid(rgrid, zgrid)  # [nz, nr]
    psi_slice = (rr - R_AXIS) ** 2 + zz**2
    times = np.array([0.5, 1.0, 1.5])
    return _EfitFluxGrid(
        times=times,
        rgrid=rgrid,
        zgrid=zgrid,
        psirz=np.stack([psi_slice] * 3, axis=0),
        psi_axis=np.zeros(3),
        psi_bdry=np.full(3, A_MINOR**2),
        fpol=None,
        bcentr=np.full(3, 5.4),
        rcentr=np.full(3, R_AXIS),
    )


def test_rho_axis_lcfs_and_band(grid):
    """rho is 0 on the axis, 1 on the LCFS, and linear in minor radius."""
    assert grid.rz2rho(R_AXIS, 0.0, 1.0)[0] == pytest.approx(0.0, abs=1e-5)
    assert grid.rz2rho(R_AXIS + A_MINOR, 0.0, 1.0)[0] == pytest.approx(1.0, abs=1e-5)
    assert grid.rz2rho(R_AXIS + 0.85 * A_MINOR, 0.0, 1.0)[0] == pytest.approx(
        0.85, abs=1e-5
    )
    assert grid.rz2rho(R_AXIS, 0.90 * A_MINOR, 1.0)[0] == pytest.approx(0.90, abs=1e-5)


def test_bt_inverse_r_falloff(grid):
    """|B_T| follows the vacuum 1/R falloff when fpol is unavailable."""
    bt = grid.rz2bt(R_AXIS + A_MINOR, 0.0, 1.0)[0]
    assert bt == pytest.approx(5.4 * R_AXIS / (R_AXIS + A_MINOR), abs=1e-5)


def test_spline_cache_is_reused(grid):
    """Repeated evaluations at one time slice build a single cached spline."""
    for _ in range(5):
        grid.psinorm(R_AXIS, 0.0, 1.0)
    assert len(grid._splines) == 1
    grid.psinorm(R_AXIS, 0.0, 0.5)
    assert len(grid._splines) == 2


def test_orient_psirz_scrambled_axes(grid):
    """Axis inference restores [nt, nz, nr] from a scrambled axis ordering."""
    scrambled = np.transpose(grid.psirz, (1, 0, 2))  # axes (Z, t, R)
    psirz, rgrid, zgrid, times = CmodEdgeCoreProfiles._orient_psirz(
        scrambled, grid.zgrid, grid.times, grid.rgrid
    )
    assert psirz.shape == (len(times), len(zgrid), len(rgrid))
    assert np.array_equal(times, grid.times)
    assert np.array_equal(rgrid, grid.rgrid)
    assert np.array_equal(zgrid, grid.zgrid)
    assert np.array_equal(psirz, grid.psirz)


def test_dimensionless_matches_analytic_formulas():
    """The dimensionless variables reproduce Maris 2025 eqs. (2)-(4) exactly."""
    t_j = 100.0 * const.e
    out = CmodEdgeCoreProfiles._dimensionless(
        n_m3=np.array([1e20]),
        t_ev=np.array([100.0]),
        a_minor=np.array([A_MINOR]),
        r_0=np.array([R_AXIS]),
        kappa=np.array([1.5]),
        ip=np.array([1.0e6]),
        bt=np.array([5.4]),
    )
    q_cyl = (
        (2 * np.pi / const.mu_0)
        * (5.4 * A_MINOR**2 / (1.0e6 * R_AXIS))
        * ((1 + 1.5**2) / 2)
    )
    assert out["q_cyl"][0] == pytest.approx(q_cyl, rel=1e-12)
    assert out["pressure"][0] == pytest.approx(2 * 1e20 * t_j, rel=1e-12)
    assert out["beta_t"][0] == pytest.approx(
        4 * const.mu_0 * 1e20 * t_j / 5.4**2, rel=1e-12
    )
    assert out["rho_star"][0] == pytest.approx(
        np.sqrt(M_DEUTERON * t_j) / (const.e * 5.4 * A_MINOR), rel=1e-12
    )
    assert np.isfinite(out["nu_star"][0])
    assert out["nu_star"][0] > 0
    assert 0 < out["beta_t"][0] < 1


def test_dimensionless_low_ip_is_nan():
    """Slices with |I_p| below MIN_IP give NaN nu*/q_cyl, not absurd values."""
    out = CmodEdgeCoreProfiles._dimensionless(
        n_m3=[1e20],
        t_ev=[100.0],
        a_minor=[A_MINOR],
        r_0=[R_AXIS],
        kappa=[1.5],
        ip=[0.0],
        bt=[5.4],
    )
    assert np.isnan(out["nu_star"][0])
    assert np.isnan(out["q_cyl"][0])
    assert np.isfinite(out["beta_t"][0])


def test_lsvm_d_is_one_on_the_boundary():
    """nu* = 3.5*beta^-0.40 maps exactly onto lsvm_d = 1 (Maris 2025 eq. 7)."""
    beta = np.array([1e-3, 3e-3, 1e-2])
    nu_limit = LSVM_D_COEFF * beta**LSVM_D_EXPONENT
    np.testing.assert_allclose(
        CmodEdgeCoreProfiles._lsvm_d(nu_limit, beta), 1.0, rtol=1e-12
    )
    doubled = CmodEdgeCoreProfiles._lsvm_d(2 * nu_limit[:1], beta[:1])
    assert doubled[0] == pytest.approx(2.0, rel=1e-12)


def test_osborne_mtanh_pedestal_shape():
    """The mtanh has the pedestal top/mid/bottom structure."""
    c_0, c_1, top, bottom = 0.98, 0.05, 100.0, 1.0
    x = np.array([c_0 - 3 * c_1, c_0, c_0 + 3 * c_1])
    vals = CmodEdgeCoreProfiles._osborne_mtanh(x, c_0, c_1, top, bottom, 0, 0, 0)
    assert vals[0] == pytest.approx(top, abs=0.5)
    assert vals[1] == pytest.approx((top + bottom) / 2, rel=1e-12)
    assert vals[2] == pytest.approx(bottom, abs=0.5)


def test_fit_band_average_recovers_band_mean():
    """A noise-free mtanh profile is fitted and band-averaged correctly."""
    truth = (0.95, 0.08, 4.0, 0.1, 0.05, 0.0, 0.0)
    psin = np.linspace(0.0, 1.2, 40)
    y = CmodEdgeCoreProfiles._osborne_mtanh(psin, *truth)
    lo2, hi2 = 0.85**2, 0.95**2
    got = CmodEdgeCoreProfiles._fit_band_average(psin, y, lo2, hi2)
    band = CmodEdgeCoreProfiles._osborne_mtanh(np.linspace(lo2, hi2, 21), *truth)
    assert got == pytest.approx(float(np.mean(band)), rel=0.05)


def test_fit_band_average_refuses_extrapolation():
    """No fit value is reported for a band the data does not reach."""
    truth = (0.95, 0.08, 4.0, 0.1, 0.0, 0.0, 0.0)
    psin = np.linspace(0.0, 0.5, 20)  # core-only coverage
    y = CmodEdgeCoreProfiles._osborne_mtanh(psin, *truth)
    assert np.isnan(CmodEdgeCoreProfiles._fit_band_average(psin, y, 0.85**2, 0.95**2))


def test_band_average_respects_band_and_validity():
    """The band average uses only finite values with rho inside the band."""
    vals = np.array([1.0, 2.0, 3.0, np.nan])
    rho = np.array([0.86, 0.90, 0.5, 0.9])
    assert CmodEdgeCoreProfiles._band_average(vals, rho, 0.85, 0.95) == pytest.approx(
        1.5
    )
    assert np.isnan(CmodEdgeCoreProfiles._band_average(vals, rho, 0.3, 0.4))


def test_coulomb_log_reasonable():
    """lnLambda is in the usual tokamak range for edge-like parameters."""
    val = CmodEdgeCoreProfiles._coulomb_log(np.array([1e20]), np.array([80.0]))[0]
    assert 8 < val < 20


def test_region_and_column_registry_consistency():
    """physics.profile_cols and profiles.all_columns agree on the 27 columns."""
    physics_cols = [
        col for cols in CmodPhysicsMethods.profile_cols.values() for col in cols
    ]
    assert len(physics_cols) == 26
    assert len(set(physics_cols)) == 26
    assert set(physics_cols + ["lsvm_d"]) == set(CmodEdgeCoreProfiles.all_columns())
    assert len(CmodEdgeCoreProfiles.all_columns()) == 27


def test_every_new_column_has_config_metadata():
    """Each of the 27 new columns has description/units (+ ordered validity)."""
    attrs = config(Tokamak.CMOD).physics.attributes
    for col in CmodEdgeCoreProfiles.all_columns():
        assert col in attrs, f"missing [cmod.physics.attributes.{col}]"
        entry = attrs[col]
        assert entry["description"]
        assert entry["units"] in {"m^-3", "eV", "Pa", "dimensionless"}
        if "validity" in entry:
            lo, hi = entry["validity"]
            assert lo < hi
