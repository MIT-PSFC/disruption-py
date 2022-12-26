from dataclasses import dataclass

from scipy.interpolate import interp1d, interp2d
from scipy.optimize import curve_fit
from scipy.signal import lfilter, medfilt
import numpy as np


def interp1(x, y, new_x, kind='linear', bounds_error=True, fill_value=0):
    set_interp = interp1d(
        x, y, kind=kind, bounds_error=bounds_error, fill_value=fill_value)
    return set_interp(new_x)


def interp2(X, Y, V, Xq, Yq, kind='linear'):
    set_interp = interp2d(X, Y, V, kind=kind)
    return set_interp(Xq, Yq)


# TODO: Implement this
def efit_rz_interp():
    pass


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
    coeffs, var_matrix = curve_fit(gauss, x, y)
    return coeffs


def gauss(x, *params):
    """ Guassian function"""
    z, mu, sigma = params
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
        d[i] = (y[i+1]-y[i-1])/(x[i+1]-x[i-1])
    return d


def fastsmooth(y, w, smooth_type=1, ends_type=0):
    smoothed_y = smooth(y, w, ends_type)
    for i in range(smooth_type-1):
        smoothed_y = smooth(smoothed_y, w, ends_type)
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


def power(a):
    # Multiplicative constants (kappa) to get the power radiating in the i^th viewing chord
    kappa = np.array([1.976e8, 2.060e8, 2.146e8, 2.319e8, 2.277e8, 2.773e8,
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
    etendu = np.array([30206, 29034, 28066, 27273, 26635, 40340, 39855, 39488,
                       39235, 39091, 39055, 39126, 7972, 8170, 8498, 7549,
                       7129, 6854, 11162, 11070, 11081, 11196, 11419, 11761,
                       29321, 28825, 28449, 28187, 28033, 7058, 7140, 7334,
                       7657, 8136, 8819, 7112, 6654, 6330, 6123, 29621,
                       29485, 29431, 29458, 29565, 29756, 30032, 30397, 6406]*1e4, dtype=np.float64)  # convert to [m^(-2)]

    @dataclass
    class Channel:
        label: str
        chanpwr: np.ndarray
        brightness: np.ndarray
        R: float
        Z: float
        angle: float

    @dataclass
    class Power:
        pwrmix: np.ndarray
        divl: np.ndarray
        divu: np.ndarray
        chan: np.ndarray
    c = Channel('', np.zeros((1, 4096)), np.zeros((1, 4096)), 0.0, 0.0, 0.0)
    b = Power(np.zeros((1, 4096)), np.zeros((1, 4096)),
              np.zeros((1, 4096)), np.tile(c, (1, 48)))
    for i in range(48):
        b.chan[i].chanpwr = kappa[i]*a.chan[i].pwr
        b.chan[i].brightness = etendu[i] * a.chan[i].pwr
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
        b.pwrmix = b.pwrmix - kappa[i]*b.divl/7.0/kappa[i+43]
    b.pwrmix = b.pwrmix + b.divl + b.divu
    return b


def get_bolo(shot_id, bol_channels, bol_prm, bol_top, bol_time, drtau=50):
    drtau /= 1.e3
    gam = np.zeros((1, 49))
    tau = np.zeros((1, 49))
    kappa = [1.2307e8, 1.4719e8, 1.5750e8, 1.5316e8, 1.6025e8, 1.6418e8,
             1.7073e8, 2.3007e8, 2.2421e8, 2.1928e8, 2.1501e8, 2.1228e8, 2.0982e8,
             2.0762e8, 2.0568e8, 1.8170e8, 1.4038e8, 1.1152e8, 0.98347e8, 0.91370e8,
             0.86452e8, 0.87789e8, 0.93361e8, 0.97920e8, 5.202e8, 5.072e8, 4.865e8,
             4.527e8, 4.180e8, 4.154e8, 4.137e8, 4.248e8, 4.459e8, 4.727e8, 5.195e8,
             6.633e7, 6.544e7, 6.699e7, 7.100e7, 1.749e8, 1.860e8, 1.939e8, 2.016e8,
             2.095e8, 2.169e8, 2.256e8, 2.311e8, 10.39e7, 4.4559e8, 4.3989e8, 4.1525e8,
             3.9185e8, 3.6744e8, 3.4002e8, 3.1598e8]
    scrfact = np.ones((1, 48))
    if shot_id > 97000:
        kappa[24:35] = [1.5398e8, 1.5013e8, 1.4461e8, 1.3399e8, 1.2371e8,
                        1.2297e8, 1.2246e8, 1.2574e8, 1.3198e8, 1.3993e8, 1.5376e8]
        # improved ECH screen
        scrfact = [0.7990, 0.8006, 0.8021, 0.8035, 0.8048, 0.8060,
                   0.8070, 0.8080, 0.8089, 0.8099, 0.8105, 0.8096,
                   0.8080, 0.8058, 0.8035, 0.8011, 0.8038, 0.8062,
                   0.8083, 0.8100, 0.8096, 0.8079, 0.8061, 0.8042,
                   0.8051, 0.8062, 0.8074, 0.8085, 0.8096, 0.8101,
                   0.8079, 0.8056, 0.8031, 0.8004, 0.7973, 0.7991,
                   0.8019, 0.8044, 0.8066, 0.8086, 0.8096, 0.8105,
                   0.8099, 0.8089, 0.8080, 0.8070, 0.8060, 0.8037]
    # Channel positions and angles of orientation
    # Reference: /fusion/pillar-archive/u/leonard/idl/bolo_geom.pro on D3D
    aperx = [234.953, 234.953, 234.953, 234.953, 234.953, 234.953, 234.953,
             234.953, 234.953, 234.953, 234.953, 234.953, 234.953, 234.953,
             234.953, 231.136, 231.136, 231.136, 231.136, 231.136, 231.136,
             231.136, 231.136, 231.136,
             231.894, 231.894, 231.894, 231.894, 231.894, 231.894, 231.894,
             231.894, 231.894, 231.894, 231.894, 234.932, 234.932, 234.932,
             234.932, 234.932, 234.932, 234.932, 234.932, 234.932, 234.932,
             234.932, 234.932, 234.932]*1e-2  # converted to meters
    apery = [72.990, 72.990, 72.990, 72.990, 72.990, 72.990, 72.990,
             72.990, 72.990, 72.990, 72.990, 72.990, 72.990, 72.990,
             72.990, 82.261, 82.261, 82.261, 82.261, 82.261, 82.261,
             82.261, 82.261, 82.261,
             -77.254, -77.254, -77.254, -77.254, -77.254, -77.254, -77.254,
             -77.254, -77.254, -77.254, -77.254, -66.881, -66.881, -66.881,
             -66.881, -66.881, -66.881, -66.881, -66.881, -66.881, -66.881,
             -66.881, -66.881, -66.881]*1e-2  # converted to meters
    angle = [269.449, 265.689, 261.923, 258.163, 254.406, 250.901, 247.898,
             244.898, 241.888, 238.875, 235.879, 232.862, 227.997, 221.241,
             214.504, 208.235, 201.115, 194.026, 187.671, 182.184, 176.689,
             171.193, 165.700, 160.222,
             213.690, 210.190, 206.690, 203.180, 199.690, 194.440, 187.440,
             180.440, 173.420, 166.420, 159.430, 155.960, 149.230, 142.420,
             135.780, 129.600, 126.600, 123.600, 120.600, 117.600, 114.600,
             111.600, 108.600, 101.910]

    @dataclass
    class Channel:
        label: str
        R: float
        Z: float
        angle: float
        ier: int
        pwr: np.ndarray
        raw: np.ndarray
        gam: float
        tau: float
        scrfact: float
    one_channel = Channel('', 0.0, 0.0, 0.0, 0, np.zeros(
        (1, 4096)), np.zeros((1, 4096)), 0.0, 0.0, 0.0)
    channels = np.tile(one_channel, (1, 48))

    @dataclass
    class Bolo:
        shot_id: int
        kappa: np.ndarray
        time: np.ndarray
        raw_time: np.ndarray
        ntimes: int
        tot_pwr: np.ndarray
        channels: np.ndarray
    bolo_shot = Bolo(shot_id, kappa, np.zeros((1, 4096)), np.zeros(
        (1, 16384)), 0, np.zeros((1, 4096)), channels)
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
    if len(bol_time) <= 16384 or bol_time[-1] <= bol_time[0] or min(np.diff(bol_time)) < -1e-5:
        for i in range(48):
            bolo_shot.channels[i].ier = 1
        return bolo_shot
    time = np.linspace(np.min(bol_time[0]), np.min(bol_time[-1]), 16384)
    dt = time[1] - time[0]
    window_size = np.around(drtau/dt)
    smoothing_kernel = (1/window_size) * np.ones((1, window_size))
    bolo_shot.ntimes = int(len(time)/4)
    bolo_shot.time = np.linspace(np.min(time), np.max(time), bolo_shot.ntimes)
    t_del = bolo_shot.time[1] - bolo_shot.time[0]
    bolo_shot.raw_time = time
    m = 2*np.fix(np.fix(1000*drtau)/np.fix(1000*t_del)/2) + 1
    k = np.arange(0, m) - np.fix((m-1)/2)
    nzer = np.where(k != 0)
    k[nzer] = 1.0 / k[nzer]
    k = k/t_del/(np.fix(m/2)*2)
    for i in range(48):
        bolo_shot.chan[i].label = bol_channels[i]
        data = interp1(bol_time, bol_top[i], bolo_shot.raw_time)
        bolo_shot.channels[i].ier = 0
        bolo_shot.channels[i].raw = data
        bolo_shot.channels[i].gam = gam[i+1]
        bolo_shot.channels[i].tau = tau[i+1]
        bolo_shot.channels[i].scrfact = scrfact[i]
        bolo_shot.channels[i].R = aperx[i]
        bolo_shot.channels[i].Z = apery[i]
        bolo_shot.channels[i].angle = angle[i]
        # Subtract baseline offset
        temp = data - np.mean(data[:20])
        # Filter signal using causal moving average filter (i.e. boxcar)
        temp_filtered = lfilter(smoothing_kernel, 1, temp)
        dr_dt = np.gradient(temp_filtered, dt)
        # Calculate power on each detector, P_d(t) [as given in Leonard et al, Rev. Sci. Instr. (1995)]
        bolo_shot.channels[i].pwr = medfilt(
            (gam[i+1]*temp_filtered + tau[i+1]*dr_dt)/scrfact[i], window_size)
    return bolo_shot
