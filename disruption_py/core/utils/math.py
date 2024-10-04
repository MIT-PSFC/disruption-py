#!/usr/bin/env python3

"""
This module contains utility functions for various numerical operations.
"""

import copy
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.signal import lfilter, medfilt

pd.options.mode.chained_assignment = None


def interp1(x, y, new_x, kind="linear", bounds_error=False, fill_value=np.nan, axis=-1):
    """
    Interpolate a 1-D array.

    This function uses scipy.interpolate.interp1d to interpolate a 1-D array
    using the given x and y values and interpolation method. It also allows
    the user to specify the value of the new_x array, and whether or not to
    raise an error if the new x-values are outside the range of the original
    x-values.

    Parameters
    ----------
    x : array
        The x-values of the original array.
    y : array
        The y-values of the original array.
    new_x : array
        The x-values of the interpolated array.
    kind : str, optional
        The interpolation method to use. Options are 'linear', 'nearest',
        'zero', 'slinear', 'quadratic', and 'cubic'. Defaults to 'linear'.
    bounds_error : bool, optional
        If True, an error is raised if new_x is outside of the range of x.
        If False, the new x-values are set to fill_value. Defaults to False.
    fill_value : float, optional
        The value to use for new_x values outside of the range of x. This
        is only used if bounds_error is False. Defaults to nan.

    Returns
    -------
    _ : array
        The interpolated y-values.
    """
    set_interp = interp1d(
        x, y, kind=kind, bounds_error=bounds_error, fill_value=fill_value, axis=axis
    )
    return set_interp(new_x)


def exp_filter(x, w, strategy="fragmented"):
    """
    Implements an exponential filter.

    This function implements an exponential filter on the given array x. In the
    case of nan values in the input array, we default to using the last timestep
    that was not a nan value. In the fragmented strategy, any time we encounter
    invald values, we restart the filter at the next valid value.

    Parameters
    ----------
    x : array
        The array to filter.
    w : float
        The filter weight.
    strategy: str, optional
        Imputation strategy to be used, if any. Options are 'fragmented' or 'none'.
        Default is 'fragmented.'

    Returns
    -------
    _ : array
        The filtered array.
    """
    filtered_x = np.zeros(x.shape)
    filtered_x[0] = x[0]
    for i in range(1, len(x)):
        filtered_x[i] = w * x[i] + (1 - w) * filtered_x[i - 1]
        if strategy == "fragmented":
            if np.isnan(filtered_x[i - 1]):
                filtered_x[i] = x[i]
    return filtered_x


def smooth(arr: np.ndarray, window_size: int) -> np.ndarray:
    """
    Implements Matlab's smooth function https://www.mathworks.com/help/curvefit/smooth.html.

    Parameters
    ----------
    arr: np.ndarray
        Array to smooth

    window_size: int
        Size of the window to smooth over

    Returns
    -------
    np.ndarray
        Smoothed array
    """
    mid = np.convolve(arr, np.ones(window_size, dtype=int), "valid") / window_size
    b_weights = np.arange(1, window_size - 1, 2)
    start = np.cumsum(arr[: window_size - 1][::2] / b_weights)
    end = (np.cumsum(arr[:-window_size:-1])[::2] / b_weights)[::-1]
    return np.concatenate((start, mid, end))


def gauss_smooth(y, smooth_width, ends_type):
    """
    Smooth a dataset using a Gaussian window.

    Parameters
    ----------
    y : array_like
        The y coordinates of the dataset.
    smooth_width : int
        The width of the smoothing window.
    ends_type : int
        Determines how the "ends" of the signal are handled.
        0 -> ends are "zeroed"
        1 -> the ends are smoothed with progressively smaller smooths the closer to the end.

    Returns
    -------
    array_like
        The smoothed dataset.
    """
    w = np.round(smooth_width)
    w = int(w)  # Ensure w is an integer
    ly = len(y)
    s = np.zeros(ly)

    for i in range(ly):
        if i < w // 2:
            if ends_type == 0:
                s[i] = 0
            else:
                s[i] = np.mean(y[: i + w // 2])
        elif i >= ly - w // 2:
            if ends_type == 0:
                s[i] = 0
            else:
                s[i] = np.mean(y[i - w // 2 :])
        else:
            s[i] = np.mean(y[i - w // 2 : i + w // 2])
    return s


def gaussian_fit(*args):
    """
    Fits a Gaussian curve to a set of data points (x,y).
    Returns an array of coefficients, c, that describe the fit.

    Reference: https://stackoverflow.com/q/11507028

    Parameters
    ----------
    x : array
        The x-coordinates of the data points.
    y : array
        The y-coordinates of the data points.
    p0 : array
        The initial values of the parameters.

    Returns
    -------
    coeffs : array
        The coefficients of the fit.
    """

    coeffs, *_ = curve_fit(gauss, *args)
    return coeffs


def gauss(x, *params):
    """
    Gaussian function.

    Parameters
    ----------
    x : array
        The x-coordinates of the data points.
    params : array
        The 3 parameters of the Gaussian function.

    Returns
    -------
    out : array
        The Gaussian function evaluated at the given x-coordinates.
    """

    a, mu, sigma = params
    out = a * np.exp(-((x - mu) ** 2) / (2.0 * sigma**2))
    return out


def gaussian_fit_with_fixed_mean(mu, *args):
    """Same as gaussian_fit() but with mu as a fixed parameter"""
    coeffs, *_ = curve_fit(lambda x, a, sigma: gauss(x, a, mu, sigma), *args)
    return coeffs


def matlab_gsastd(
    x, y, derivative_mode, width, smooth_type=1, ends_type=0, slew_rate=0
):
    """
    Python reimplementation of the GSASTD routine originally by Alessandro Pau
    Fast non-causal differentiation of noisy data.

    Parameters
    ----------
    x : array_like
        The x coordinates of the dataset.
    y : array_like
        The y coordinates of the dataset.
    derivative_mode : int
        If 0, the standard deviation of the dataset is calculated. If 1, the
        standard deviation of the first derivative of the dataset is
        calculated.
    width : int
        The width of the smoothing window.
    smooth_type : int, optional
        Determines the type of smoothing to use.
        0 -> no smoothing.
        1 -> rectangular (sliding-average or boxcar)
        2 -> triangular (2 passes of sliding-average)
        3 -> pseudo-Gaussian (3 passes of sliding-average)
    ends_type : int, optional
        Determines how the "ends" of the signal are handled.
        0 -> ends are "zeroed"
        1 -> the ends are smoothed with progressively smaller smooths the closer to the end.
    slew_rate : float, optional
        The slew rate of the dataset. If specified, the dataset is
        processed to remove any points that exceed the slew rate.

    Returns
    -------
    array_like
        Smoothed dataset.

    References
    -------
    - https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/get_P_ohm_d3d.m

    Last major update by: William Wei on 7/24/2024
    """
    y_arr = y.copy()
    if slew_rate:
        for i in range(0, len(y_arr) - 1):
            diff = y_arr[i + 1] - y_arr[i]
            if abs(diff) > slew_rate:
                y_arr[i + 1] = y_arr[i] + np.sign(diff) * slew_rate
    if smooth_type == 0:
        width = 1
    if derivative_mode == 0:
        return matlab_fastsmooth(y_arr, width, smooth_type, ends_type)
    if derivative_mode == 1:
        return matlab_fastsmooth(matlab_deriv(x, y_arr), width, smooth_type, ends_type)
    raise ValueError("derivative_mode only takes 0 or 1 as input")


def matlab_deriv(x, y):
    """
    Calculate the derivative of a dataset. Used by GSASTD.

    Parameters
    ----------
    x: array
        The x-coordinates of the data points.
    y: array
        The y-coordinates of the data points.

    Returns
    -------
    d: array
        The derivative of the dataset.

    References:
    -------
    - https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/GSASTD.m
    """
    n = len(y)
    d = np.zeros(y.shape)
    d[0] = (y[1] - y[0]) / (x[1] - x[0])
    d[n - 1] = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2])
    for i in range(1, n - 1):
        d[i] = (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1])
    return d


def matlab_fastsmooth(y, w, smooth_type=1, ends_type=0):
    """
    Wrapper for looping over the dataset using the smooth function.
    Used by GSASTD.

    Parameters
    ----------
    y : array_like
        The y coordinates of the dataset.
    w : int
        The width of the smoothing window.
    smooth_type : int, optional
        Determines the type of smoothing to use.
        0 -> no smoothing.
            -- NOTE: In the original MATLAB function, smooth_type = 0 will still
                cause the function to smoooth the array once.
            -- (i.e. it's identical to smooth_type = 1).
            -- We keep the same behavior in this python implementation.
        1 -> rectangular (sliding-average or boxcar)
        2 -> triangular (2 passes of sliding-average)
        3 -> pseudo-Gaussian (3 passes of sliding-average)
    ends_type : int, optional
        Determines how the "ends" of the signal are handled.
        0 -> ends are "zeroed"
        1 -> the ends are smoothed with progressively smaller smooths the closer to the end.

    Returns
    -------
    array_like
        The smoothed dataset.

    References:
    -------
    - https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/GSASTD.m

    Last major update by William Wei on 7/24/2024
    """
    smoothed_y = matlab_sa(y, w, ends_type)
    for _ in range(smooth_type - 1):
        smoothed_y = matlab_sa(smoothed_y, w, ends_type)
    return smoothed_y


def matlab_sa(y, smooth_width, ends_type=0):
    """
    Compute the centered moving average of y
    Used by GSASTD

    Parameters
    ----------
    y : array_like
        The y coordinates of the dataset.
    smooth_width : int
        The width of the smoothing window.
    ends_type : int
        Determines how the "ends" of the signal are handled.
        0 -> ends are "zeroed"
        1 -> the ends are smoothed with progressively smaller smooths the closer to the end.

    Returns
    -------
    array_like
        The smoothed dataset.

    References:
    -------
    - https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/GSASTD.

    Last major update by William Wei on 7/24/2024
    """
    # NOTE: numpy behaviour is different than matlab and will round X.5 to nearest
    # even value instead of value farther away from 0

    w = matlab_round_int(smooth_width)
    sum_points = np.sum(y[:w])
    s = np.zeros(y.shape)
    half_w = matlab_round_int(w / 2.0)
    ly = len(y)
    for i in range(ly - w):
        s[i + half_w - 1] = sum_points
        sum_points = sum_points - y[i]
        sum_points = sum_points + y[i + w]
    s[i + half_w] = np.sum(y[ly - w : ly])
    y_smooth = s / w

    if ends_type == 1:
        start_point = matlab_round_int((smooth_width + 1) / 2)
        y_smooth[0] = (y[0] + y[1]) / 2
        for i in range(1, start_point):
            y_smooth[i] = np.mean(y[0 : 2 * i + 1])
            y_smooth[ly - i - 1] = np.mean(y[ly - 2 * i - 1 : ly])
        y_smooth[ly - 1] = (y[ly - 1] + y[ly - 2]) / 2

    return y_smooth


def matlab_round_int(x):
    """
    Custom rounding function. Round x.5 to the nearest integer with larger magnitude
    numpy behaviour is different than matlab and will round X.5 to nearest even value
    instead of value farther away from 0.
    """
    sign = np.sign(x)
    x = abs(x)
    if x % 1 == 0.5:
        return int(sign * np.ceil(x))
    return int(sign * round(x))


# TODO: Cover documentation with Cristina
def matlab_power(a):
    """
    Python reimplementation of powers_new.m
    Used for calculating the total radiated power from the bolometer array.

    -------
    References:
    - https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/powers_new.m

    Last major update by William Wei on 7/26/2024
    """

    # Multiplicative constants (kappa) to get the power radiating in the i^th viewing chord
    kappa = np.array(
        [
            1.976e8,
            2.060e8,
            2.146e8,
            2.319e8,
            2.277e8,
            2.773e8,
            2.845e8,
            2.877e8,
            2.780e8,
            2.692e8,
            2.519e8,
            2.509e8,
            1.004e8,
            0.935e8,
            0.893e8,
            0.790e8,
            0.711e8,
            0.683e8,
            0.861e8,
            0.834e8,
            0.812e8,
            0.860e8,
            0.941e8,
            0.991e8,
            1.605e8,
            1.588e8,
            1.592e8,
            1.561e8,
            1.496e8,
            0.700e8,
            0.698e8,
            0.699e8,
            0.719e8,
            0.774e8,
            0.865e8,
            0.692e8,
            0.683e8,
            0.698e8,
            0.739e8,
            1.742e8,
            1.771e8,
            1.898e8,
            1.935e8,
            1.994e8,
            2.035e8,
            2.033e8,
            2.030e8,
            0.848e8,
            4.632e8,
            4.370e8,
            4.230e8,
            3.792e8,
            3.565e8,
            3.166e8,
            2.836e8,
        ],
        dtype=np.float64,
    )
    # Multiplicative constants (etendu, or K from Leonard 1995 formula 2) to get the
    # brightness for the i^th viewing chord
    # [See bolo_geometric_values_d3d.xlsx or Tony Leonard's prad.pro routine]
    etendu = (
        np.array(
            [
                30206,
                29034,
                28066,
                27273,
                26635,
                40340,
                39855,
                39488,
                39235,
                39091,
                39055,
                39126,
                7972,
                8170,
                8498,
                7549,
                7129,
                6854,
                11162,
                11070,
                11081,
                11196,
                11419,
                11761,
                29321,
                28825,
                28449,
                28187,
                28033,
                7058,
                7140,
                7334,
                7657,
                8136,
                8819,
                7112,
                6654,
                6330,
                6123,
                29621,
                29485,
                29431,
                29458,
                29565,
                29756,
                30032,
                30397,
                6406,
            ],
            dtype=np.float64,
        )
        * 1.0e4
    )  # convert to [m^(-2)]

    @dataclass
    # pylint: disable-next=missing-class-docstring
    class Channel:
        label: str
        chanpwr: np.ndarray
        brightness: np.ndarray
        r: float
        z: float
        angle: float

    @dataclass
    # pylint: disable-next=missing-class-docstring
    class Power:
        pwrmix: np.ndarray
        divl: np.ndarray
        divu: np.ndarray
        chan: list

    b = Power(
        pwrmix=np.zeros((4096)),
        divl=np.zeros((4096)),
        divu=np.zeros((4096)),
        chan=[],
    )
    for i in range(48):
        b.chan.append(Channel("", np.zeros((4096)), np.zeros((4096)), 0.0, 0.0, 0.0))
        b.chan[i].chanpwr = kappa[i] * a.channels[i].pwr
        b.chan[i].brightness = etendu[i] * a.channels[i].pwr
        b.chan[i].r = a.channels[i].r
        b.chan[i].z = a.channels[i].z
        b.chan[i].angle = a.channels[i].angle
    b.pwrmix = 0.0
    b.divl = 0.0
    b.divu = 0.0
    # Calculate power radiated from lower divertor region
    for i in range(24, 31):
        b.divl = b.divl + b.chan[i].chanpwr
    # Calculate power radiated from upper divertor region
    for i in range(21, 24):
        b.divu = b.divu + b.chan[i].chanpwr
    # Calculate total radiated power (based on Tony Leonard's IDL code)
    for i in range(21):
        b.pwrmix = b.pwrmix + b.chan[i].chanpwr
    for i in range(5, 12):
        b.pwrmix = b.pwrmix - kappa[i] * b.divl / 7.0 / kappa[i + 43]
    b.pwrmix = b.pwrmix + b.divu + b.divl
    return b


def matlab_get_bolo(shot_id, bol_channels, bol_prm, bol_top, bol_time, drtau=50):
    """
    Python reimplementation of getbolo_new.m.
    Used for calculating the total radiated power from the bolometer array.

    -------
    References:
    - https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/getbolo_new.m

    Last major update by William Wei on 7/26/2024
    """
    drtau /= 1e3
    # NOTE: why set gam, tau & scrfact as 2D arrays?
    gam = np.zeros((1, 49))
    tau = np.zeros((1, 49))
    kappa = [
        1.2307e8,
        1.4719e8,
        1.5750e8,
        1.5316e8,
        1.6025e8,
        1.6418e8,
        1.7073e8,
        2.3007e8,
        2.2421e8,
        2.1928e8,
        2.1501e8,
        2.1228e8,
        2.0982e8,
        2.0762e8,
        2.0568e8,
        1.8170e8,
        1.4038e8,
        1.1152e8,
        0.98347e8,
        0.91370e8,
        0.86452e8,
        0.87789e8,
        0.93361e8,
        0.97920e8,
        5.202e8,
        5.072e8,
        4.865e8,
        4.527e8,
        4.180e8,
        4.154e8,
        4.137e8,
        4.248e8,
        4.459e8,
        4.727e8,
        5.195e8,
        6.633e7,
        6.544e7,
        6.699e7,
        7.100e7,
        1.749e8,
        1.860e8,
        1.939e8,
        2.016e8,
        2.095e8,
        2.169e8,
        2.256e8,
        2.311e8,
        10.39e7,
        4.4559e8,
        4.3989e8,
        4.1525e8,
        3.9185e8,
        3.6744e8,
        3.4002e8,
        3.1598e8,
    ]
    scrfact = np.ones((1, 48))
    if shot_id > 97000:
        kappa[24:35] = [
            1.5398e8,
            1.5013e8,
            1.4461e8,
            1.3399e8,
            1.2371e8,
            1.2297e8,
            1.2246e8,
            1.2574e8,
            1.3198e8,
            1.3993e8,
            1.5376e8,
        ]
        # transmission factor through ECH screen
        # improved ECH screen
        scrfact = [
            0.7990,
            0.8006,
            0.8021,
            0.8035,
            0.8048,
            0.8060,
            0.8070,
            0.8080,
            0.8089,
            0.8099,
            0.8105,
            0.8096,
            0.8080,
            0.8058,
            0.8035,
            0.8011,
            0.8038,
            0.8062,
            0.8083,
            0.8100,
            0.8096,
            0.8079,
            0.8061,
            0.8042,
            0.8051,
            0.8062,
            0.8074,
            0.8085,
            0.8096,
            0.8101,
            0.8079,
            0.8056,
            0.8031,
            0.8004,
            0.7973,
            0.7991,
            0.8019,
            0.8044,
            0.8066,
            0.8086,
            0.8096,
            0.8105,
            0.8099,
            0.8089,
            0.8080,
            0.8070,
            0.8060,
            0.8037,
        ]
    # Channel positions and angles of orientation
    # Reference: /fusion/pillar-archive/u/leonard/idl/bolo_geom.pro on D3D
    aperx = (
        np.array(
            [
                234.953,
                234.953,
                234.953,
                234.953,
                234.953,
                234.953,
                234.953,
                234.953,
                234.953,
                234.953,
                234.953,
                234.953,
                234.953,
                234.953,
                234.953,
                231.136,
                231.136,
                231.136,
                231.136,
                231.136,
                231.136,
                231.136,
                231.136,
                231.136,
                231.894,
                231.894,
                231.894,
                231.894,
                231.894,
                231.894,
                231.894,
                231.894,
                231.894,
                231.894,
                231.894,
                234.932,
                234.932,
                234.932,
                234.932,
                234.932,
                234.932,
                234.932,
                234.932,
                234.932,
                234.932,
                234.932,
                234.932,
                234.932,
            ]
        )
        * 1e-2
    )  # converted to meters
    apery = (
        np.array(
            [
                72.990,
                72.990,
                72.990,
                72.990,
                72.990,
                72.990,
                72.990,
                72.990,
                72.990,
                72.990,
                72.990,
                72.990,
                72.990,
                72.990,
                72.990,
                82.261,
                82.261,
                82.261,
                82.261,
                82.261,
                82.261,
                82.261,
                82.261,
                82.261,
                -77.254,
                -77.254,
                -77.254,
                -77.254,
                -77.254,
                -77.254,
                -77.254,
                -77.254,
                -77.254,
                -77.254,
                -77.254,
                -66.881,
                -66.881,
                -66.881,
                -66.881,
                -66.881,
                -66.881,
                -66.881,
                -66.881,
                -66.881,
                -66.881,
                -66.881,
                -66.881,
                -66.881,
            ]
        )
        * 1e-2
    )  # converted to meters
    angle = np.array(
        [
            269.449,
            265.689,
            261.923,
            258.163,
            254.406,
            250.901,
            247.898,
            244.898,
            241.888,
            238.875,
            235.879,
            232.862,
            227.997,
            221.241,
            214.504,
            208.235,
            201.115,
            194.026,
            187.671,
            182.184,
            176.689,
            171.193,
            165.700,
            160.222,
            213.690,
            210.190,
            206.690,
            203.180,
            199.690,
            194.440,
            187.440,
            180.440,
            173.420,
            166.420,
            159.430,
            155.960,
            149.230,
            142.420,
            135.780,
            129.600,
            126.600,
            123.600,
            120.600,
            117.600,
            114.600,
            111.600,
            108.600,
            101.910,
        ]
    )

    @dataclass
    # pylint: disable-next=missing-class-docstring
    class Channel:
        label: str
        r: float
        z: float
        angle: float
        ier: int
        pwr: np.ndarray
        raw: np.ndarray
        gam: float
        tau: float
        scrfact: float

    one_channel = Channel(
        label="",
        r=0.0,
        z=0.0,
        angle=0.0,
        ier=0,
        pwr=np.zeros((1, 4096)),
        raw=np.zeros((1, 16384)),
        gam=0.0,
        tau=0.0,
        scrfact=0.0,
    )
    channels = [copy.deepcopy(one_channel) for i in range(48)]

    @dataclass
    # pylint: disable-next=missing-class-docstring
    class Bolo:
        shot_id: int
        kappa: np.ndarray
        time: np.ndarray
        raw_time: np.ndarray
        ntimes: int
        tot_pwr: np.ndarray
        channels: list

    bolo_shot = Bolo(
        shot_id=shot_id,
        kappa=kappa,
        time=np.zeros((1, 4096)),
        raw_time=np.zeros((1, 16384)),
        ntimes=0,
        tot_pwr=np.zeros((1, 4096)),
        channels=channels,
    )
    # TODO: Find explanation for this
    if shot_id > 79400:
        gam = bol_prm[:49]
        tau = bol_prm[49:98]
    # The bolometry signals and their common timebase should each have 163840
    # samples, digitized at 10 kHz.  Tony Leonard's IDL procedure down-samples
    # these to 16384 samples using IDL's REBIN function, and then again to just
    # 4096 samples.  But REBIN does a forward averaging, which introduces some
    # amount of non-causality.  Instead, for the timebase I will generate a new
    # array with 16384 points, covering the same range of time.  For the
    # bolometry signals, I will use nearest-neighbor linear interpolation to
    # down-sample.

    # Sometimes the bolo data is garbage.  I need to execute a few checks to
    # decide if the data is valid or not.  If not, set 'ier' = 1 in the
    # appropriate field for each channel.
    if (
        len(bol_time) <= 16384
        or bol_time[-1] <= bol_time[0]
        or min(np.diff(bol_time)) < -1e-5
    ):
        print("Bolo data is garbage")
        for i in range(48):
            bolo_shot.channels[i].ier = 1
        return bolo_shot

    # time = np.linspace(np.min(bol_time[0]), np.min(bol_time[-1]), 16384)
    time = np.linspace(bol_time[0], bol_time[-1], 16384)
    dt = time[1] - time[0]
    window_size = matlab_round_int(drtau / dt)
    smoothing_kernel = (1.0 / window_size) * np.ones(window_size)

    bolo_shot.ntimes = int(len(time) / 4)
    bolo_shot.time = np.linspace(np.min(time), np.max(time), bolo_shot.ntimes)
    bolo_shot.raw_time = time

    for i in range(48):
        bolo_shot.channels[i].label = bol_channels[i]
        data = interp1(bol_time, bol_top[i], bolo_shot.raw_time)

        bolo_shot.channels[i].ier = 0
        bolo_shot.channels[i].raw = data
        bolo_shot.channels[i].gam = gam[i + 1]
        bolo_shot.channels[i].tau = tau[i + 1]
        bolo_shot.channels[i].scrfact = scrfact[i]
        bolo_shot.channels[i].r = aperx[i]
        bolo_shot.channels[i].z = apery[i]
        bolo_shot.channels[i].angle = angle[i]

        # Subtract baseline offset
        temp = data - np.mean(data[:20])
        # Filter signal using causal moving average filter (i.e. boxcar)
        # NOTE: lfilter gives closer results to MATLAB than np.convolve
        # temp_filtered = np.convolve(temp, smoothing_kernel, "same")
        temp_filtered = lfilter(smoothing_kernel, 1, temp)
        dr_dt = np.gradient(temp_filtered, time)
        # Calculate power on each detector, P_d(t) [as given in Leonard et al,
        # Rev. Sci. Instr. (1995)]
        # BUG: MATLAB's medfilt1 function does not require window_size to be odd.
        # This could introduce discrepancy with this python implementation.
        bolo_shot.channels[i].pwr = medfilt(
            (gam[i + 1] * temp_filtered + tau[i + 1] * dr_dt) / scrfact[i],
            window_size + (not window_size % 2),
        )

    return bolo_shot


def matlab_gradient_1d_vectorized(f, h, **_kwargs):
    """
    Compute the gradient for a 1D array using vectorized operations.

    :param f: Input 1D array
    :param h: Spacing array with the same length as f
    :return: Gradient of f
    """
    f = np.array(f)
    h = np.array(h)

    h_diff = np.diff(h)
    f_diff = np.diff(f)
    # Combine into a single gradient array
    g = np.empty_like(f)
    g[0] = f_diff[0] / h_diff[0]
    g[-1] = f_diff[-1] / h_diff[-1]
    g[1:-1] = (f[2:] - f[0:-2]) / (h[2:] - h[0:-2])

    return g
