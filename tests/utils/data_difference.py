#!/usr/bin/env python3

"""
Module for handling and analyzing data differences between fresh MDSplus data
and cached SQL data.
"""

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import xarray as xr

from disruption_py.config import config
from disruption_py.core.utils.misc import safe_cast


@dataclass
class DataDifference:
    """
    Data difference between fresh MDSplus data and cached SQL data.
    """

    shot_id: int
    data_column: str

    missing_cache_data: bool
    missing_fresh_data: bool

    anomalies: np.ndarray = field(init=False)  # 1 if anomaly, 0 o.w.
    relative_difference: np.ndarray = field(init=False)
    fresh_data_array: xr.DataArray
    cache_data_array: xr.DataArray
    expect_failure: bool

    fresh_time: xr.DataArray
    cache_time: xr.DataArray

    def __post_init__(self):
        self.anomalies, self.relative_difference = self.compute_numeric_anomalies()

    @property
    def num_anomalies(self) -> int:
        """Sum the number of anomalies"""
        return np.sum(self.anomalies)

    @property
    def timebase_length(self) -> int:
        """Get the timebase length"""
        return len(self.anomalies)

    @property
    def missing_data(self) -> bool:
        """Return True if either fresh or cache is missing data."""
        return self.missing_cache_data or self.missing_fresh_data

    @property
    def failed(self) -> str:
        """Return True if missing data or if there are too many anomalies."""
        if self.missing_data:
            return True
        return (
            self.num_anomalies / self.timebase_length
            > 1 - config().tests.match_fraction
        )

    @property
    def failure_ratio_string(self) -> str:
        """Create a string for the fraction of datapoints that are anomalies"""
        return f"{self.num_anomalies / self.timebase_length:.4f}"

    @property
    def column_mismatch_string(self) -> str:
        """Create a string showing the difference between fresh and cache data."""
        # Missing data handled here because difference_df expects data to exist
        s = f"Shot {self.shot_id} column {self.data_column}"
        if self.missing_cache_data or self.missing_fresh_data:
            fresh_str = (
                "Missing fresh data" if self.missing_fresh_data else "Have fresh data"
            )
            cache_str = (
                "Missing cache data" if self.missing_cache_data else "Have cache data"
            )
            return f"{s}: {fresh_str} and {cache_str}"
        return s + f" with arrays:\n{self.difference_df.to_string()}"

    @property
    def difference_df(self) -> pd.DataFrame:
        """
        Create a dataframe with columns for time, fresh data, cache data, the
        ratio between the two data, the absolute difference between them, the relative
        difference, and whether the point is an anomaly.
        """
        indices = (
            np.arange(self.timebase_length)
            if config().tests.verbose_output
            else self.anomalies.flatten()
        )
        anomaly = self.anomalies[indices]
        if config().tests.verbose_output:
            fresh_data = self.fresh_data_array
            cache_data = self.cache_data_array
            fresh_time = self.fresh_time
        else:
            fresh_time = self.fresh_time.where(self.anomalies, drop=True)
            fresh_data = self.fresh_data_array.where(self.anomalies, drop=True)
            cache_data = self.cache_data_array.where(self.anomalies, drop=True)
        return pd.DataFrame(
            {
                "Time": fresh_time,
                "Fresh Data": fresh_data,
                "Cache Data": cache_data,
                "Fresh/Cache": fresh_data / cache_data,
                "Absolute difference": abs(fresh_data - cache_data),
                "Relative difference": self.relative_difference[indices],
                "Anomaly": anomaly,
            }
        )

    def compute_numeric_anomalies(self):
        """
        Get the indices of the data where numeric differences exist between the
        cached and fresh data.
        """

        # handle missing data case
        if self.missing_fresh_data or self.missing_cache_data:
            if self.missing_fresh_data and self.missing_cache_data:
                missing_timebase_length = 0
            elif self.missing_fresh_data:
                missing_timebase_length = len(self.cache_data_array)
            else:
                missing_timebase_length = len(self.fresh_data_array)
            return np.ones(missing_timebase_length, dtype=bool), np.zeros(
                missing_timebase_length
            )

        cache_is_nan = self.cache_data_array.isnull()
        fresh_is_nan = self.fresh_data_array.isnull()

        # handle case where both arrays are all null
        if cache_is_nan.all() and fresh_is_nan.all():
            return np.zeros(len(self.fresh_data_array), dtype=bool), np.zeros(
                len(self.fresh_data_array)
            )

        relative_difference = safe_cast(
            np.where(
                self.cache_data_array != 0,
                np.abs(
                    (self.fresh_data_array - self.cache_data_array)
                    / self.cache_data_array
                ),
                np.where(self.fresh_data_array != 0, np.inf, np.nan),
            ),
            "float64",
        )  # necessary in case all produced values are nan

        numeric_anomalies_mask = np.where(
            np.isnan(relative_difference),
            False,
            relative_difference > config().tests.val_tolerance,
        )

        nan_anomalies_mask: np.ndarray = (cache_is_nan != fresh_is_nan).values
        anomalies: np.ndarray = numeric_anomalies_mask | nan_anomalies_mask
        return anomalies, relative_difference


def assert_frame_equal_unordered(df1: pd.DataFrame, df2: pd.DataFrame):
    """Compare whether two dataframes have the same values."""
    df1_sorted = df1.sort_values(by=config().inout.sql.protected_columns).reset_index(
        drop=True
    )
    df2_sorted = df2.sort_values(by=config().inout.sql.protected_columns).reset_index(
        drop=True
    )
    pd.testing.assert_frame_equal(df1_sorted, df2_sorted, check_like=True)
