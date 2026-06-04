import os
from glob import glob
import numpy as np
import pandas as pd

def read_csv_files(dir, file_pattern):
    csv_files = glob(os.path.join(dir, file_pattern))
    dfs = [pd.read_csv(file) for file in csv_files]
    return pd.concat(dfs, ignore_index=True)

def write_timing_results(df):
    print(f"Timing Avgs:\n---------------")
    print(f"Reading MDSplus: {np.nanmean(df['read_mds']):.6} s")
    print(f"Calculating autocorrelation: {np.nanmean(df['autocorr']):.6} s")
    print(f"Butterworth Filter: {np.nanmean(df['bworth']):.6} s")
    print(f"Finding TQ: {np.nanmean(df['find_tq']):.6} s")
    print(f"Total: {np.nanmean(df['total']):.6} s")
    print("---------------")


if __name__=='__main__':
    timing_dir = '/home/henrycw/projects/label-thermal-quench/disruption-py-label-thermal-quench/disruption-py/tq_timing'
    df_05 = read_csv_files(timing_dir, 'timing_105*.csv')
    print(f"\nTiming results for 2005 shots")
    write_timing_results(df_05)
    df_12_16 = read_csv_files(timing_dir, 'timing_11*.csv')
    print(f"\nTiming results for 2012-2016 shots")
    write_timing_results(df_12_16)
    
