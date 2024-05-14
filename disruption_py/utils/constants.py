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
CMOD_TEST_SHOTS = {
    "flattop1_fast": 1150805012,
    "no_disrup1_full": 1150805013,
    "no_disrup2_full": 1150805014,
    "rampdown1_full": 1150805015,
    "rampdown2_full": 1150805016,
    "rampdown3_full": 1150805017,
    "rampdown4_full": 1150805019,
    "rampdown5_full": 1150805020,
    "rampdown6_full": 1150805021,
    "flattop2_full": 1150805022,
}
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
CMOD_EXPECTED_FAILURE_COLUMNS = {
    'Te_width': [-1], 'z_error': [-1], 'z_prog': [-1], 'zcur': [-1], 
    'v_z': [-1], 'z_times_v_z': [-1], 
    'dipprog_dt': [-1], # constant factor scaling error
    'ip_error': [-1], # constant error
    'sxr': [-1], 'tritop': [-1], 'tribot': [-1], 'a_minor': [-1], 
    'rmagx': [-1], 'chisq': [-1], 'V_surf': [-1], 'ne_peaking': [-1], 
    'Te_peaking': [-1], 'pressure_peaking': [-1], 
    'Greenwald_fraction': [1150805012, 1150805016, 1150805017, 1150805019, 1150805020, 1150805021, 1150805022], 
    'n_equal_1_phase': [-1], 'BT': [-1], 'prad_peaking': [-1], 
    'dip_smoothed': [-1], 'ip_prog': [-1], 'p_input': [-1]
}


D3D_TEST_SHOTS = {
    "disrup1_fast": 161228,
    "disrup2_full": 161237,
    "no_disrup1_full": 166177,
    "no_disrup2_full": 166253,
}
D3D_TEST_COLUMNS = [
    "H98", "ip", "q95", "squareness", "zcur_normalized", "q0", "ip_error", "beta_p",
    "time_until_disrupt", "z_error", "li_RT", "beta_p_RT", "n1rms_normalized", "qstar",
    "ip_error_RT", "H_alpha", "n_e_RT", "Wmhd_RT", "ip_RT", "dli_dt", "dbetap_dt",
    "dn_dt", "shot", "n_equal_1_mode", "dip_dt", "upper_gap", "n_equal_1_normalized",
    "q95_RT", "zcur", "lower_gap", "Greenwald_fraction_RT", "kappa", "kappa_area",
    "power_supply_railed", "n_e", "delta", "Greenwald_fraction", "dWmhd_dt", "Wmhd",
    "aminor", "time", "li", "beta_n", "dipprog_dt", "dipprog_dt_RT"
]
D3D_EXPECTED_FAILURE_COLUMNS = {
    'kappa': [161237], 'H_alpha': [-1], 'dipprog_dt_RT': [161228, 161237], 'li': [161237], 
    'dWmhd_dt': [161237], 'beta_p': [161237], 'dn_dt': [161237], 'Greenwald_fraction': [161237, 166177], 
    'li_RT': [-1], 'n_equal_1_normalized': [-1], 'zcur': [-1], 'dli_dt': [161237], 
    'Greenwald_fraction_RT': [-1], 'H98': [-1], 'q95_RT': [-1], 
    'zcur_normalized': [-1], 'qstar': [161237], 'Wmhd': [161237], 
    'lower_gap': [161237], 'beta_p_RT': [-1], 'n1rms_normalized': [-1], 
    'dbetap_dt': [161237], 'n_equal_1_mode': [-1], 'q95': [161237], 'upper_gap': [161237], 
    'q0': [161237], 'n_e': [161237], 'beta_n': [161237], 'kappa_area': [-1]
}

