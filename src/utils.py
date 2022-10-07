from scipy.interpolate import interp1d
import numpy as np

def interp1(x,y, new_x, kind='linear',bounds_error = True, fill_value = 0):
    set_interp = interp1d(x,y,kind=kind,bounds_error = bounds_error,fill_value = fill_value)
    return set_interp(new_x)

def smooth(arr, window_size): 
    """ 
    Implements Matlab's smooth function https://www.mathworks.com/help/curvefit/smooth.html.
    """
    mid = np.convolve(arr, np.ones(window_size, dtype=int),'valid')/window_size 
    b_weights = np.arange(1,window_size -1, 2)
    start = np.cumsum(arr[:window_size-1][::2]/b_weights)
    end = (np.cumsum(arr[:-window_size:-1])[::2]/b_weights)[::-1]
    return np.concatenate((start, mid, end))