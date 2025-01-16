#!/usr/bin/env python3

"""Example usage of `get_shots_data` demonstrating using a custom time setting."""


from disruption_py.settings import RetrievalSettings, TimeSetting, TimeSettingParams
from disruption_py.workflow import get_shots_data


class PRadTime(TimeSetting):
    """Class for retrieving prad times"""

    def _get_times(self, params: TimeSettingParams):
        """Return prad times"""
        (time_array,) = params.mds_conn.get_dims(
            r"\twopi_diode", tree_name="spectroscopy"
        )
        time_array = time_array[time_array > 0]
        return time_array


retrieval_settings = RetrievalSettings(time_setting=PRadTime())

shot_data = get_shots_data(
    tokamak="cmod",
    shotlist_setting=[1150805012],
    retrieval_settings=retrieval_settings,
)
