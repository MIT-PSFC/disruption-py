#!/usr/bin/env python3

from disruption_py.workflow import get_shots_data
from disruption_py.settings import (
    LogSettings,
    TimeSetting,
    TimeSettingParams,
    RetrievalSettings,
)


class PRadTime(TimeSetting):
    def _get_times(self, params: TimeSettingParams):
        (time_array,) = params.mds_conn.get_dims(
            r"\twopi_diode", tree_name="spectroscopy"
        )
        time_array = time_array[time_array > 0]
        return time_array


shot_settings = RetrievalSettings(
    time_setting=PRadTime(),
    run_tags=[],  # default is all, so if you do not want all data must set to an empty list
    run_columns=["p_rad"],  # run parameter methods
    only_requested_columns=True,  # only return the column data for run_columns
    signal_domain="flattop",  # retrieve data from the flattop time domain
)

shot_data = get_shots_data(
    tokamak="cmod",
    shotlist_setting="cmod_non_disruptions_ids_not_blacklist_mini",
    shot_settings=shot_settings,
    output_setting="examples/p_rad.csv",
    num_processes=4,
    log_settings=LogSettings(
        log_to_console=False,  # don't log to the console
        log_file_path="examples/last_log.log",  # instead log to this file
    ),
)
