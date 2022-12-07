from scipy.interpolate import interp1d, interp2d
from scipy.optimize import curve_fit
import numpy as np


def interp1(x, y, new_x, kind='linear', bounds_error=True, fill_value=0):
    set_interp = interp1d(
        x, y, kind=kind, bounds_error=bounds_error, fill_value=fill_value)
    return set_interp(new_x)


def interp2(X, Y, V, Xq, Yq, kind='linear'):
    set_interp = interp2d(X, Y, V, kind=kind)
    return set_interp(Xq, Yq)


def smooth(arr, window_size):
    """ 
    Implements Matlab's smooth function https://www.mathworks.com/help/curvefit/smooth.html.
    """
    mid = np.convolve(arr, np.ones(window_size, dtype=int),
                      'valid')/window_size
    b_weights = np.arange(1, window_size - 1, 2)
    start = np.cumsum(arr[:window_size-1][::2]/b_weights)
    end = (np.cumsum(arr[:-window_size:-1])[::2]/b_weights)[::-1]
    return np.concatenate((start, mid, end))

# Credit to: https://stackoverflow.com/questions/11507028/fit-a-gaussian-function


def gaussian_fit(x, y):
    coeffs, var_matrix = curve_fit(guass, x, y)
    return coeffs


def gauss(x, *params):
    """ Guassian function"""
    z, mu, sigma = p
    return z*np.exp(-(x-mu)**2/(2.0*sigma**2))

# Alessandro Pau (JET & AUG) has given Cristina a robust routine that
# performs time differentiation with smoothing, while preserving causality.
# It can be useful for differentiating numerous signals such as Ip, Vloop,
# etc.  It is called 'GSASTD'. We will use this routine in place of Matlab's
# 'gradient' and smoothing/filtering routines for certain signals.


def gsastd(x, y, derivative_mode, width, smooth_type=1, ends_type=0, slew_rate=None):
    if slew_rate is not None:
        for i in range(1, len(y)-1):
            diff = y[i+1] - y[i]
            if abs(diff) > slew_rate:
                y[i+1] = y[i] + np.sign(diff)*slew_rate
    if smooth_type == 0:
        width = 0
    if derivative_mode == 0:
        return fastsmooth(y, width, smooth_type, ends_type)
    elif derivative_mode == 1:
        return fastsmooth(deriv(x, y), width, smooth_type, ends_type)
    raise ValueError("derivative_mode only takes 0 or 1 as input")


def deriv(x, y):
    n = len(y)
    d = np.zeros(y.shape)
    d[0] = (y[1]-y[0])/(x[1]-x[0])
    d[n] = (y[n]-y[n-1])/(x[n]-x[n-1])
    for i in range(1, n-1):
        d[i] = (y[i+1]-y[i-1])/(x[j+1]-x[j-1])
    return d


def fastsmooth(y, w, smooth_type=1, ends_type=0):
    smoothed_y = smooth(y, w, ends_type)
    for i in range(smooth_type-1):
        smoothed_y = smooth(smoothed_y, w, ends)
    return smoothed_y


def smooth(y, smooth_width, ends_type):
    # TODO: Ask why this is. NOTE: numpy behaviour is different than matlab and will round X.5 to nearest even value instead of value farther away from 0
    w = np.round(smooth_width)
    sum_points = np.sum(y[:w])
    s = np.zeros(y.shape)
    half_w = np.round(w/2.0)
    l = len(y)
    for i in range(l-w):
        s[i+half_w-1] = sum_points
        sum_points = sum_points - y[i]
    s[i+half_w] = np.sum(y[l-w:l])
    return s/w

# TODO: Add get_bolowew and power_new
