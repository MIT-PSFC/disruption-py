import xarray as xr
import numpy as np
from scipy.io import loadmat 


def apply_freq_correction(ds: xr.Dataset) -> xr.Dataset:
    complex_fft = ds["mirnov_fft_real"] + 1j * ds["mirnov_fft_imag"]

    sensor_dim = "sensor"
    freq_dim = next((d for d in complex_fft.dims if "freq" in d.lower()), None)
    if freq_dim is None:
        raise KeyError(
            f"Could not find frequency dimension in {complex_fft.dims}. "
            "Expected a dimension name containing 'freq'."
        )

    sensor_name_key = None
    for candidate in ("sensor_names", "sensor_name"):
        if candidate in ds.coords or candidate in ds.variables:
            sensor_name_key = candidate
            break
    if sensor_name_key is None:
        raise KeyError("Could not find sensor name field. Expected 'sensor_names' or 'sensor_name'.")

    other_dims = [d for d in complex_fft.dims if d not in (sensor_dim, freq_dim)]
    corrected = complex_fft.transpose(sensor_dim, freq_dim, *other_dims).copy(deep=True)

    freq_vals = ds.coords[freq_dim].values
    sensor_indices = ds.coords[sensor_dim].values
    sensor_name_arr = ds[sensor_name_key]

    for sensor_idx in sensor_indices:
        sensor_name = sensor_name_arr.sel({sensor_dim: sensor_idx}).item()
        if isinstance(sensor_name, bytes):
            sensor_name = sensor_name.decode()

        _, H = __cal_Correction(str(sensor_name), freq_vals)
        if np.any(np.abs(H) == 0):
            raise ZeroDivisionError(
                f"Calibration transfer function contains zeros for sensor {sensor_name}."
            )

        corrected.loc[{sensor_dim: sensor_idx}] = (
            corrected.sel({sensor_dim: sensor_idx})
            / xr.DataArray(H, dims=[freq_dim], coords={freq_dim: freq_vals})
        )

    corrected = corrected.transpose(*complex_fft.dims)
    ds["mirnov_fft_real"] = corrected.real
    ds["mirnov_fft_imag"] = corrected.imag

    return ds


  
def __cal_Correction(sensor_name,freq):
    CAL_PATH = '/mnt/home/sears/Matlab/Calibration/cal_v2/'
    try:mat = loadmat(CAL_PATH+'451_responses/'+sensor_name +'_cal.mat')
    except: mat = loadmat(CAL_PATH+'451_responses/'+'BP1T_GHK' +'_cal.mat')
    f = mat['f'][0]; H_spline = mat['H_spline'][0]
    try:
        mat = loadmat(CAL_PATH+'fine_tuning/'+sensor_name +'_cal+.mat')
        fine_cal= mat['fine_cal'][0,0]
        mat = loadmat(CAL_PATH+'fine_tuning/'+sensor_name +'_caln.mat')
        fine_caln = mat['fine_caln'][0,0]
    except: fine_cal = fine_caln = 1
    #return f,H_spline,fine_cal,fine_caln 
    
    
    return  f, np.interp(freq,f,H_spline) * fine_cal * fine_caln