# C-Mod `disruption_warning` dataset

The `disruption_warning` SQL table for C-Mod contains important disruption parameters for a large number of shots.
See [SQL Database](../../api/data_io/sql_database.md) for how to retrieve data from it.

The dataset contains unique plasma discharges from MIT's Alcator C-Mod tokamak, from the 2012 to 2016 experimental campaigns, plus additional discharges from 2005.
??? info "Available columns on CMod"

	```
	'dbkey', 'shot', 'time', 'time_until_disrupt', 'ip_error', 'dip_dt',
	'beta_p', 'beta_n', 'li', 'n_equal_1_normalized', 'z_error', 'v_z',
	'z_times_v_z', 'kappa', 'pressure_peaking', 'H98', 'q0', 'qstar', 'q95',
	'v_0', 'v_mid', 'v_edge', 'dn_dt', 'p_rad_slow', 'p_oh_slow', 'p_icrf',
	'p_lh', 'radiated_fraction', 'power_supply_railed', 'v_loop_efit',
	'r_dd', 'lower_gap', 'upper_gap', 'dbetap_dt', 'dli_dt', 'ip', 'zcur',
	'n_e', 'dipprog_dt', 'v_loop', 'p_rad', 'p_oh', 'ssep', 'dWmhd_dt',
	'dprad_dt', 'v_0_uncalibrated', 'Te_width', 'Greenwald_fraction',
	'intentional_disruption', 'Te_width_ECE', 'Wmhd', 'n_over_ncrit',
	'n_equal_1_mode', 'Mirnov', 'Mirnov_norm_btor', 'Mirnov_norm_bpol',
	'Te_peaking', 'ne_peaking', 'Te_peaking_ECE', 'SXR_peaking',
	'kappa_area', 'I_efc', 'SXR', 'H_alpha', 'Prad_peaking_CVA'
	```
