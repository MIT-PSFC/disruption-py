"""
Script for plotting test results from tests/test_thermal_quench_times.py
Author: Henry Wietfeldt
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df_old = pd.read_csv('test_thermal_quench_results_cq_onset90.csv')
df_new = pd.read_csv('test_thermal_quench_results_cq_onset95.csv')
#df_90 = pd.read_csv('test_thermal_quench_results_max90.csv')

dfs = [df_old]
labels = ['0.9*max']
bins = np.linspace(-0.004, 0.004, 60)

for i, df in enumerate(dfs):
    # Print summary statistics
    error = df['onset_error_s'].to_numpy()
    print(f"Stats for {labels[i]}")
    print(f"Mean |Error| = {1e3*np.mean(np.abs(error)):.3f} ms")
    print(f"Median |Error| = {1e3*np.median(np.abs(error)):.3f} ms")
    print(f"Std Dev |Error| = {1e3*np.std(np.abs(error)):.3f} ms")
    print(f"Min Error = {1e3*np.min(error):.3f} ms")
    print(f"Max Error = {1e3*np.max(error):.3f} ms")
    print(f"Num Outliers (outisde TQ [start, end] by >1 ms) = {np.sum(~df['within_tq_range'])} out of {len(df)} shots\n")
    print(df[~df['within_tq_range']].sort_values(by='onset_error_s'))

    plt.hist(df['onset_error_s'], edgecolor='k', bins=bins, label=labels[i], alpha=0.5)

plt.title("Test vs 100 Manually Labeled Shots", fontsize=18)
plt.xlabel('Error in TQ Onset Time (Auto - Manual) [s]', fontsize=16)
plt.xlim(-0.003, 0.00325)
plt.ylabel('Count', fontsize=16)
plt.legend()
plt.show()