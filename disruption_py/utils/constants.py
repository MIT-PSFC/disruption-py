TIME_CONST = 1e-6
MAX_THREADS_PER_SHOT = 10
MAX_PROCESSES = 10

BASE_PROTECTED_COLUMNS = ['time', 'shot']
CMOD_PROTECTED_COLUMNS = ['dbkey', 'shot', 'time', 'time_until_disrupt', 'ip_error', 'dip_dt', 'beta_p', 'beta_n', 'li', 'n_equal_1_normalized', 'z_error', 'v_z', 'z_times_v_z', 'kappa', 'pressure_peaking', 'H98', 'q0', 'qstar', 'q95', 'dn_dt', 'p_rad_slow', 'p_oh_slow', 'p_icrf', 'p_lh', 'radiated_fraction', 'power_supply_railed', 'v_loop_efit', 'lower_gap', 'upper_gap', 'dbetap_dt', 'dli_dt',
						  'ip', 'zcur', 'n_e', 'dipprog_dt', 'v_loop', 'p_rad', 'p_oh', 'ssep', 'dWmhd_dt', 'dprad_dt', 'Te_width', 'Greenwald_fraction', 'intentional_disruption', 'Te_width_ECE', 'Wmhd', 'n_over_ncrit', 'n_equal_1_mode', 'Mirnov', 'Mirnov_norm_btor', 'Mirnov_norm_bpol', 'Te_peaking', 'ne_peaking', 'Te_peaking_ECE', 'SXR_peaking', 'kappa_area', 'I_efc', 'SXR', 'H_alpha', 'Prad_peaking_CVA', 'commit_hash']
D3D_PROTECTED_COLUMNS = []

DEFAULT_THRESHOLD = 0.35  # Time until disrupt threshold for binary classification
DEFAULT_RATIO = .2  # Ratio of test data to total data and validation data to train data
BLACK_WINDOW_THRESHOLD = 5.e-3 # A 'black window' threshold [s]; obscures input data from a window in time on disruptive shots during trianing/testing

MAX_SHOT_TIME = 7.0  # [s] <-- used to detect if shot times are using ms

# Used for testing
VAL_TOLERANCE = 0.01   # Tolerance for comparing values between MDSplus and SQL
MATCH_FRACTION = 0.95   # Fraction of signals that must match between MDSplus and SQL
VERBOSE_OUTPUT = False

DEFAULT_COLS = ['time', 'time_until_disrupt','shot']
PAPER_COLS = [
    'aminor',
    'n_e',
    'ip',
    'delta',
    'li',
    'Wmhd',
    'kappa',
    'squareness',
    'shot',
    'time_until_disrupt',
    'time',
]
DERIVED_PAPER_COLS = [
    'ip-exp-10-none',
    'ip-exp-50-none',
]


# Shot list used for testing
# Mix of disruptive and non-disruptive shots present in SQL and MDSplus
CMOD_TEST_SHOTS = [
    1150805012,   # Flattop Disruption
    1150805013,     # No Disruption
    1150805014,     # No Disruption
    1150805015,     # Rampdown Disruption
    1150805016,     # Rampdown Disruption
    1150805017,     # Rampdown Disruption
    1150805019,     # Rampdown Disruption
    1150805020,     # Rampdown Disruption
    1150805021,     # Rampdown Disruption
    1150805022      # Flattop Disruption
]
CMOD_TEST_COLUMNS = [
    'I_efc', 'sxr', 'time_until_disrupt', 'beta_n', 'beta_p', 'kappa', 'li',
    'upper_gap', 'lower_gap', 'q0', 'qstar', 'q95', 'v_loop_efit', 'Wmhd',
    'ssep', 'n_over_ncrit', 'tritop', 'tribot', 'a_minor', 'rmagx', 'chisq',
    'dbetap_dt', 'dli_dt', 'dWmhd_dt', 'V_surf', 'kappa_area', 'Te_width',
    'ne_peaking', 'Te_peaking', 'pressure_peaking', 'n_e', 'dn_dt',
    'Greenwald_fraction', 'n_equal_1_mode', 'n_equal_1_normalized',
    'n_equal_1_phase', 'BT', 'prad_peaking', 'v_0', 'ip', 'dip_dt',
    'dip_smoothed', 'ip_prog', 'dipprog_dt', 'ip_error', 'z_error',
    'z_prog', 'zcur', 'v_z', 'z_times_v_z', 'p_oh', 'v_loop', 'p_rad',
    'dprad_dt', 'p_lh', 'p_icrf', 'p_input', 'radiated_fraction', 'time',
    'shot', 'commit_hash'
]
CMOD_EXPECTED_FAILURE_COLUMNS = [
    'lower_gap', 'upper_gap', 'ssep', 'dipprog_dt', 'n_over_ncrit', # constant factor scaling error
    'ip_error' # constant error
]


D3D_TEST_SHOTS = [
    161228, # disruptive
    161237, # disruptive
    166177, # non disruptive 
    166253
]
D3D_TEST_COLUMNS = [
    'shot', 'time', 'time_until_disrupt', 'ip_error', 'dip_dt',
    'beta_p', 'beta_n', 'li', 'n_equal_1_mode_IRLM', 'z_error', 'v_z',
    'kappa', 'H98', 'q0', 'qstar', 'q95', 'dn_dt', 'radiated_fraction',
    'power_supply_railed', 'lower_gap', 'upper_gap', 'dbetap_dt', 'dli_dt',
    'ip', 'zcur', 'n_e', 'dipprog_dt', 'v_loop', 'p_rad', 'dWmhd_dt',
    'dprad_dt', 'p_nbi', 'p_ech', 'p_ohm', 'intentional_disruption',
    'Greenwald_fraction', 'Te_HWHM', 'other_hardware_failure', 'Te_HWHM_RT',
    'v_loop_RT', 'n_e_RT', 'Greenwald_fraction_RT', 'ip_error_RT', 'ip_RT',
    'dipprog_dt_RT', 'Wmhd_RT', 'Wmhd', 'n_equal_1_mode',
    'n_equal_1_normalized', 'Te_width_normalized', 'Te_width_normalized_RT',
    'q95_RT', 'li_RT', 'beta_p_RT', 'oeamp1em', 'oeamp1om', 'oefrq1em',
    'oefrq1om', 'oeamp1e', 'oeamp1o', 'oefrq1e', 'oefrq1o', 'delta',
    'squareness', 'zcur_normalized', 'aminor', 'n1rms_normalized',
    'kappa_area', 'Te_peaking_CVA_RT', 'ne_peaking_CVA_RT',
    'Prad_peaking_CVA_RT', 'Prad_peaking_XDIV_RT', 'H_alpha',
]
D3D_EXPECTED_FAILURE_COLUMNS = []
