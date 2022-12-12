from collections import namedtuple

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
# Multiplicative constants (kappa) to get the power radiating in the i^th viewing chord
KAPPA = np.array([1.976e8, 2.060e8, 2.146e8, 2.319e8, 2.277e8, 2.773e8,
                  2.845e8, 2.877e8, 2.780e8, 2.692e8, 2.519e8, 2.509e8,
                  1.004e8, 0.935e8, 0.893e8, 0.790e8, 0.711e8, 0.683e8,
                  0.861e8, 0.834e8, 0.812e8, 0.860e8, 0.941e8, 0.991e8,
                  1.605e8, 1.588e8, 1.592e8, 1.561e8, 1.496e8, 0.700e8,
                  0.698e8, 0.699e8, 0.719e8, 0.774e8, 0.865e8, 0.692e8,
                  0.683e8, 0.698e8, 0.739e8, 1.742e8, 1.771e8, 1.898e8,
                  1.935e8, 1.994e8, 2.035e8, 2.033e8, 2.030e8, 0.848e8,
                  4.632e8, 4.370e8, 4.230e8, 3.792e8, 3.565e8, 3.166e8, 2.836e8], dtype=np.float64)
# Multiplicative constants (etendu, or K from Leonard 1995 formula 2) to get the
# brightness for the i^th viewing chord
# [See bolo_geometric_values_d3d.xlsx or Tony Leonard's prad.pro routine]
ETENDU = np.array([30206, 29034, 28066, 27273, 26635, 40340, 39855, 39488,
                   39235, 39091, 39055, 39126, 7972, 8170, 8498, 7549,
                   7129, 6854, 11162, 11070, 11081, 11196, 11419, 11761,
                   29321, 28825, 28449, 28187, 28033, 7058, 7140, 7334,
                   7657, 8136, 8819, 7112, 6654, 6330, 6123, 29621,
                   29485, 29431, 29458, 29565, 29756, 30032, 30397, 6406]*1e4, dtype=np.float64)  # convert to [m^(-2)]


def power(a):
    Channel = namedtuple(
        'Channel', ['label', 'chanpwr', 'brightness', 'R', 'Z', 'angle'])
    Power = namedtuple('Power', ['pwrmix', 'divl', 'divu', 'chan'])
    c = Channel('', np.zeros((1, 4096)), np.zeros((1, 4096)), 0.0, 0.0, 0.0)
    b = Power(np.zeros((1, 4096)), np.zeros((1, 4096)),
              np.zeros((1, 4096)), np.tile(c, (1, 48)))
    for i in range(48):
        b.chan[i].chanpwr = KAPPA[i]*a.chan[i].pwr
        b.chan[i].brightness = ETENDU[i] * a.chan[i].pwr
        b.chan[i].R = a.chan[i].R
        b.chan[i].Z = a.chan[i].Z
        b.chan[i].angle = a.chan[i].angle
    b.pwrmix = 0.0
    b.divl = 0.0
    b.divu = 0.0
    # Calculate power radiated from lower divertor region
    for i in range(24, 31):
        b.div1 = b.div1 + b.chan[i].chanpwr
    # Calculate power radiated from upper divertor region
    for i in range(21, 24):
        b.divu = b.divu + b.chan[i].chanpwr
    # Calculate total radiated power (based on Tony Leonard's IDL code)
    for i in range(21):
        b.pwrmix = b.pwrmix + b.chan[i].chanpwr
    for i in range(5, 12):
        b.pwrmix = b.pwrmix - KAPPA[i]*b.divl/7.0/KAPPA[i+43]
    b.pwrmix = b.pwrmix + b.divl + b.divu
    return b
