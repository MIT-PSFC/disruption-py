#!/usr/bin/env python3

TIME_CONST = 1e-6
MAX_THREADS_PER_SHOT = 10
MAX_PROCESSES = 10

BASE_PROTECTED_COLUMNS = ["time", "shot"]
CMOD_PROTECTED_COLUMNS = [
    "dbkey",
    "shot",
    "time",
    "time_until_disrupt",
    "ip_error",
    "dip_dt",
    "beta_p",
    "beta_n",
    "li",
    "n_equal_1_normalized",
    "z_error",
    "v_z",
    "z_times_v_z",
    "kappa",
    "pressure_peaking",
    "H98",
    "q0",
    "qstar",
    "q95",
    "dn_dt",
    "p_rad_slow",
    "p_oh_slow",
    "p_icrf",
    "p_lh",
    "radiated_fraction",
    "power_supply_railed",
    "v_loop_efit",
    "lower_gap",
    "upper_gap",
    "dbetap_dt",
    "dli_dt",
    "ip",
    "zcur",
    "n_e",
    "dipprog_dt",
    "v_loop",
    "p_rad",
    "p_oh",
    "ssep",
    "dWmhd_dt",
    "dprad_dt",
    "Te_width",
    "Greenwald_fraction",
    "intentional_disruption",
    "Te_width_ECE",
    "Wmhd",
    "n_over_ncrit",
    "n_equal_1_mode",
    "Mirnov",
    "Mirnov_norm_btor",
    "Mirnov_norm_bpol",
    "Te_peaking",
    "ne_peaking",
    "Te_peaking_ECE",
    "SXR_peaking",
    "kappa_area",
    "I_efc",
    "SXR",
    "H_alpha",
    "Prad_peaking_CVA",
    "commit_hash",
]
D3D_PROTECTED_COLUMNS = []

DEFAULT_THRESHOLD = 0.35  # Time until disrupt threshold for binary classification
DEFAULT_RATIO = (
    0.2  # Ratio of test data to total data and validation data to train data
)
BLACK_WINDOW_THRESHOLD = 5.0e-3  # A 'black window' threshold [s]; obscures input data from a window in time on disruptive shots during trianing/testing

MAX_SHOT_TIME = 7.0  # [s] <-- used to detect if shot times are using ms

# Used for testing
VAL_TOLERANCE = 0.01  # Tolerance for comparing values between MDSplus and SQL
MATCH_FRACTION = 0.95  # Fraction of signals that must match between MDSplus and SQL
VERBOSE_OUTPUT = False

DEFAULT_COLS = ["time", "time_until_disrupt", "shot"]
PAPER_COLS = [
    "aminor",
    "n_e",
    "ip",
    "delta",
    "li",
    "Wmhd",
    "kappa",
    "squareness",
    "shot",
    "time_until_disrupt",
    "time",
]
DERIVED_PAPER_COLS = [
    "ip-exp-10-none",
    "ip-exp-50-none",
]


TEST_SHOTS = {
    "cmod": {
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
    },
    "d3d": {
        "disrup1_fast": 161228,
        "disrup2_full": 161237,
        "no_disrup1_full": 166177,
        "no_disrup2_full": 166253,
    },
}

TEST_COLUMNS = {
    "cmod": [
        "I_efc",
        "sxr",
        "time_until_disrupt",
        "beta_n",
        "beta_p",
        "kappa",
        "li",
        "upper_gap",
        "lower_gap",
        "q0",
        "qstar",
        "q95",
        "v_loop_efit",
        "Wmhd",
        "ssep",
        "n_over_ncrit",
        "tritop",
        "tribot",
        "a_minor",
        "rmagx",
        "chisq",
        "dbetap_dt",
        "dli_dt",
        "dWmhd_dt",
        "V_surf",
        "kappa_area",
        "Te_width",
        "ne_peaking",
        "Te_peaking",
        "pressure_peaking",
        "n_e",
        "dn_dt",
        "Greenwald_fraction",
        "n_equal_1_mode",
        "n_equal_1_normalized",
        "n_equal_1_phase",
        "BT",
        "prad_peaking",
        "v_0",
        "ip",
        "dip_dt",
        "dip_smoothed",
        "ip_prog",
        "dipprog_dt",
        "ip_error",
        "z_error",
        "z_prog",
        "zcur",
        "v_z",
        "z_times_v_z",
        "p_oh",
        "v_loop",
        "p_rad",
        "dprad_dt",
        "p_lh",
        "p_icrf",
        "p_input",
        "radiated_fraction",
        "time",
        "shot",
        "commit_hash",
    ],
    "d3d": [
        "H98",
        "ip",
        "q95",
        "squareness",
        "zcur_normalized",
        "q0",
        "ip_error",
        "beta_p",
        "time_until_disrupt",
        "z_error",
        "li_RT",
        "beta_p_RT",
        "n1rms_normalized",
        "qstar",
        "ip_error_RT",
        "H_alpha",
        "n_e_RT",
        "Wmhd_RT",
        "ip_RT",
        "dli_dt",
        "dbetap_dt",
        "dn_dt",
        "shot",
        "n_equal_1_mode",
        "dip_dt",
        "upper_gap",
        "n_equal_1_normalized",
        "q95_RT",
        "zcur",
        "lower_gap",
        "Greenwald_fraction_RT",
        "kappa",
        "kappa_area",
        "power_supply_railed",
        "n_e",
        "delta",
        "Greenwald_fraction",
        "dWmhd_dt",
        "Wmhd",
        "aminor",
        "time",
        "li",
        "beta_n",
        "dipprog_dt",
        "dipprog_dt_RT",
    ],
}

EXPECTED_FAILURE_COLUMNS = {
    "cmod": [
        "Te_width",
        "z_error",
        "z_prog",
        "zcur",
        "v_z",
        "z_times_v_z",
        "dipprog_dt",
        "ip_error",
        "sxr",
        "tritop",
        "tribot",
        "a_minor",
        "rmagx",
        "chisq",
        "V_surf",
        "ne_peaking",
        "Te_peaking",
        "pressure_peaking",
        "Greenwald_fraction",
        "n_equal_1_phase",
        "BT",
        "prad_peaking",
        "dip_smoothed",
        "ip_prog",
        "p_input",
        "commit_hash",
    ],
    "d3d": [
        "kappa",
        "H_alpha",
        "dipprog_dt_RT",
        "li",
        "dWmhd_dt",
        "beta_p",
        "dn_dt",
        "Greenwald_fraction",
        "li_RT",
        "n_equal_1_normalized",
        "zcur",
        "dli_dt",
        "Greenwald_fraction_RT",
        "H98",
        "q95_RT",
        "zcur_normalized",
        "qstar",
        "Wmhd",
        "lower_gap",
        "beta_p_RT",
        "n1rms_normalized",
        "dbetap_dt",
        "n_equal_1_mode",
        "q95",
        "upper_gap",
        "q0",
        "n_e",
        "beta_n",
        "kappa_area",
        "commit_hash",
    ],
}

WRITE_DATABASE_TABLE_NAME = "disruption_warning_test"
DATABASE_CONSTANTS = {
    "cmod": {
        "profile_path": "~/D3DRDB.sybase_login",
        "driver": "ODBC Driver 18 for SQL Server",
        "host": "alcdb2",
        "port": 1433,
        "db_name": "logbook",
        "protected_columns": CMOD_PROTECTED_COLUMNS,
    },
    "d3d": {
        "profile_path": "~/logbook.sybase_login",
        "driver": "FreeTDS",
        "host": "d3drdb",
        "port": 8001,
        "db_name": "d3drdb",
        "protected_columns": D3D_PROTECTED_COLUMNS,
        "additional_databases": {
            "code_rundb": {
                "profile_path": "~/logbook.sybase_login",
                "driver": "FreeTDS",
                "host": "d3drdb",
                "port": 8001,
                "db_name": "code_rundb",
                "protected_columns": [],
            }
        },
    },
}

MDSPLUS_CONNECTION_STRING_CONSTANTS = {
    "cmod": "alcdata-archives",
    "d3d": "atlas",
}
