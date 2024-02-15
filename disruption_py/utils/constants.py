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
