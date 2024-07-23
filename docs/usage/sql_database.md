DisruptionPy uses access to logbook sql databses for convenience when retrieving data from MDSPlus. Users may also use DisruptionPy to directly retrieve data from the logbook database and specifically the `disruption_warning` table.

## The `disruption_warning` table { .doc .doc-heading }
The `disruption_warning` sql tables for CMod and D3D are datasets containing important disruption paramaters for a large number of shots from the respective tokamaks.

### CMod Dataset
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
	'kappa_area', 'I_efc', 'SXR', 'H_alpha', 'Prad_peaking_CVA',
	'commit_hash'
	```

	For more details on computed values please see [parameter reference][disruption-parameter-descriptions].

## Retrieving data from the SQL database { .doc .doc-heading }
Here is an example that uses disruption_py to get shot data from the `disruption_warning` table
for eight shots from the disruption warning shotlist:
```python
--8<--
examples/database_example.py
--8<--
```

## Database Class Reference { .doc .doc-heading }

<!-- ::: disruption_py.databases.database
    handler: python
	options:
	  heading_level: 2
	  show_root_heading: false
	  show_root_toc_entry: false
	  filters: ["!^_[^_]"] -->
        
