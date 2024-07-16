#!/usr/bin/env python3

DEFAULT_THRESHOLD = 0.35  # Time until disrupt threshold for binary classification
DEFAULT_RATIO = (
    0.2  # Ratio of test data to total data and validation data to train data
)
BLACK_WINDOW_THRESHOLD = 5.0e-3  # A 'black window' threshold [s]; obscures input data from a window in time on disruptive shots during trianing/testing
DEFAULT_COLS = ["time", "time_until_disrupt", "shot"]
PAPER_COLS = [
    "aminor",
    "n_e",
    "ip",
    "delta",
    "li",
    "wmhd",
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
