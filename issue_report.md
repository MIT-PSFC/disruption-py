In order to test how different precission casting affect disruption-py calculation we run the following experiment.

Download data from a set of 10429 shots in disruption warning table, casting the dataset to 32 bit precission, 64 bits and no performing any casting.

Then differences between datasets are computed and localized. For example.
```python
# For each shot. Compute absolute differences
df32_shot = df32[df32['shot'] == shot_number]
df64_shot = df64[df64['shot'] == shot_number]

# And for each signal
y64 = df64_shot[signal].values
y32 = df32_shot[signal].values 

# Compute  difference
diff_signal = np.abs(np.abs(y64) - np.abs(y32))

# Locate max difference
max_diff_index = diff_signal.argmax()
```
This way we localize the worst offender for each field out of all shots.

Our finding can be summarize as follows.

- There is no practical difference between the data casted to 32 or 64 bits.
- There is no practical difference between the data casted to 32 and no casting at all.
- Calculated fields like tau are not significantly affected by different casting. 

    $$
    \tau = \frac{w_{mhd}}{p_{rad} - \frac{dw_{mhd}}{dt}}
    $$

- We found a bug with te_width_ece, te_width fields that can be adressed separately.

The worst error is found in tau_rad and calculated tau.

![Worst errors per variable](./images/worst_offenders_abs.png)


As demonstrated this errors don't imply siginificant differences for the calculation. 

The experiment is repeated with the 10429 shots in disruption warning table.

![tau rad worst error](./images/tau_rad_1150626012.png)

A further inspection reveals that these are the outlier cases (1 out of 10429) ! See explore_tau_rad_signal.ipynb



We found out that different nan distributions are generated for the following cases.


| Shot            | Variable          | df32_nan | df64_nan | dataNone_nan |
|-----------------|-------------------|----------|----------|--------------|
| 1050620004      | pressure_peaking  | 77       | 76       | 77           |
| 1120207015      | te_width          | 1        | 0        | 1            |
| 1120221031      | radiated_fraction | 0        | 1        | 0            |
| 1120720021      | prad_peaking      | 3        | 4        | 3            |
| 1120731024      | prad_peaking      | 3        | 4        | 3            |
| 1120821029      | prad_peaking      | 3        | 4        | 3            |
| 1120917009      | prad_peaking      | 4        | 5        | 4            |
| 1140402012      | prad_peaking      | 5        | 6        | 5            |
| 1140403023      | te_width          | 15       | 13       | 15           |
| 1140821002      | prad_peaking      | 2        | 3        | 2            |
| 1140826002      | prad_peaking      | 2        | 3        | 2            |
| 1140827029      | prad_peaking      | 7        | 8        | 7            |
| 1150710003      | prad_peaking      | 6        | 7        | 6            |
| 1160505017      | radiated_fraction | 4        | 5        | 4            |
| 1160512016      | te_width          | 13       | 11       | 13           |
| 1160615016      | te_width          | 0        | 1        | 0            |
| 1160914014      | prad_peaking      | 6        | 7        | 6            |

Here a plot of the data. Arrows on axis indicate NaNs.


![nan discrepancies](./images/nan_differences.png)


To reproduce the experiment, run the following scripts in the branch prepared for this issue:
```bash
python run_experiments_more_data.py
```
Then examine the notebooks to reproduce the analysis like: 
```bash
code ./disruption-py/explore_differences_worst_offender_without_te_buggy_signal_more_data.ipynb
```
.


Fun fact

There is a shot with infinity radiated power!
![infinity](./images/infinity.png)

P_input is to close to 0 at t~1.49

$$
\text{Radiated Fraction} = \frac{P_{\text{rad}}}{P_{\text{input}}}
$$

![infinity](./images/p_rad_and_p_input.png)