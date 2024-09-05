#!/usr/bin/env python3
import json
from itertools import combinations

import numpy as np
import pytest
import xarray as xr
from disruption_py.machine.cmod.tearing import CmodTearingMethods, get_expected_cross_phase


def test_expected_cross_phase():
    """Test to ensure the expected cross phase is calculated correctly.
    Uses synthetic data from POPSIM to compare the expected cross phase with what is actually present
    """

    probe_locations = json.load(open("disruption-py/tests/data/probe_info.json"))

    # This is where each mode is in frequency space so we can compare the expected cross phase at that point
    # Numbers obtained from looking at a spectrogram of the data
    frequency_space = {
        (1, 1): 6000,
        (2, 1): 5750,
        (3, 2): 17250,
    }

    mode_lists = [
        [(1, 1)],
        [(2, 1), (3, 2)],
    ]
    for modes in mode_lists:
        file_string = ""
        for mode in modes:
            file_string += f"_{mode[0]}_{mode[1]}"
        probe_signals = xr.open_dataset(f"disruption-py/tests/data/tearing_sim{file_string}_sfft.nc", chunks={})

        # Get each unique pair of probes
        included_probes = {var.split(".")[0] for var in probe_signals.data_vars}
        probe_pairs = list(combinations(included_probes, 2))

        # Check that we can get the right cross phase for each mode
        for mode in modes:
            frequency_index = np.abs(probe_signals.frequency - frequency_space[mode]).argmin()
            selected_signals = probe_signals.isel(frequency=frequency_index, time=-1)

            for probe_pair in probe_pairs:
                probe_1 = probe_pair[0]
                probe_2 = probe_pair[1]
                # Get the measured cross phase
                probe_1_signal = selected_signals[f"{probe_1}.real"] + 1j * selected_signals[f"{probe_1}.imag"]
                probe_2_signal = selected_signals[f"{probe_2}.real"] + 1j * selected_signals[f"{probe_2}.imag"]

                complex_cross_spectrum = probe_1_signal * probe_2_signal.conj()
                actual_cross_phase = np.angle(complex_cross_spectrum) % (2 * np.pi)
                actual_cross_phase_deg = np.degrees(actual_cross_phase)
                calculated_cross_phase = get_expected_cross_phase(
                    mode,
                    probe_locations[probe_1]["theta"],
                    probe_locations[probe_1]["phi"],
                    probe_locations[probe_2]["theta"],
                    probe_locations[probe_2]["phi"],
                )
                calculated_cross_phase_deg = np.degrees(calculated_cross_phase)

                # The expected cross phase should be within 10 degrees of the actual cross phase
                assert np.abs(actual_cross_phase_deg - calculated_cross_phase_deg) < 10, f"Actual: {actual_cross_phase_deg}, Calculated: {calculated_cross_phase_deg}"