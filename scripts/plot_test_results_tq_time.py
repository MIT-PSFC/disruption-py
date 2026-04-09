"""
Script for plotting test results from tests/test_thermal_quench_times.py
Author: Henry Wietfeldt
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('test_thermal_quench_results.csv')

# Print summary statistics
error = df['onset_error_s'].to_numpy()
print(f"Mean |Error| = {1e3*np.mean(np.abs(error)):.3f} ms")
print(f"Median |Error| = {1e3*np.median(np.abs(error)):.3f} ms")
print(f"Std Dev |Error| = {1e3*np.std(np.abs(error)):.3f} ms")
print(f"Min Error = {1e3*np.min(error):.3f} ms")
print(f"Max Error = {1e3*np.max(error):.3f} ms")
print(f"Num Outliers (|error| > 1 ms) = {np.sum(~df['within_tq_range'])} out of {len(df)} shots")

plt.hist(df['onset_error_s'], bins=50)
plt.xlabel('Error in TQ Onset Time (Auto - Manual) [s]', fontsize=16)
plt.ylabel('Count', fontsize=16)
plt.show()