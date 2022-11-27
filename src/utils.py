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
