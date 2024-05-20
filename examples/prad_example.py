from disruption_py.handlers import CModHandler
from disruption_py.settings import (
    LogSettings,
    SetTimesRequest,
    SetTimesRequestParams,
    ShotSettings,
)

cmod_handler = CModHandler()


class PRadTime(SetTimesRequest):
    def _get_times(self, params: SetTimesRequestParams):
        (time_array,) = params.mds_conn.get_dims(
            r"\twopi_diode", tree_name="spectroscopy"
        )
        time_array = time_array[time_array > 0]
        return time_array


shot_settings = ShotSettings(
    set_times_request=PRadTime(),
    run_tags=[],  # default is all, so if you do not want all data must set to an empty list
    run_columns=["p_rad"],  # run parameter methods
    only_requested_columns=True,  # only return the column data for run_columns
    signal_domain="flattop",  # retrieve data from the flattop time domain
    log_settings=LogSettings(
        log_to_console=False,  # don't log to the console
        log_file_path="examples/last_log.log",  # instead log to this file
    ),
)

shot_data = cmod_handler.get_shots_data(
    shot_ids_request="cmod_non_disruptions_ids_not_blacklist_mini",
    shot_settings=shot_settings,
    output_type_request="examples/p_rad.csv",
    num_processes=4,
)
