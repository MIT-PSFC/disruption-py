#!/usr/bin/env python3

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from disruption_py.config import config
from disruption_py.core.utils.misc import safe_cast


@dataclass
class DataDifference:
    """
    Data difference between MDSplus and SQL.
    """

    shot_id: int
    data_column: str

    missing_sql_data: bool
    missing_mdsplus_data: bool

    anomalies: np.ndarray = field(init=False)  # 1 if anomaly, 0 o.w.
    relative_difference: np.ndarray = field(init=False)
    mdsplus_column_data: pd.Series
    sql_column_data: pd.Series
    expect_failure: bool

    mds_time: pd.Series
    sql_time: pd.Series

    def __post_init__(self):
        self.anomalies, self.relative_difference = self.compute_numeric_anomalies()

    @property
    def num_anomalies(self) -> int:
        return np.sum(self.anomalies)

    @property
    def timebase_length(self) -> int:
        return len(self.anomalies)

    @property
    def missing_data(self) -> bool:
        return self.missing_sql_data or self.missing_mdsplus_data

    @property
    def failed(self) -> str:
        if self.missing_data:
            return True
        return (
            self.num_anomalies / self.timebase_length
            > 1 - config().testing.MATCH_FRACTION
        )

    @property
    def failure_ratio_string(self) -> str:
        return f"{self.num_anomalies / self.timebase_length:.4f}"

    @property
    def column_mismatch_string(self) -> str:
        return f"Shot {self.shot_id} column {self.data_column} with arrays:\n{self.difference_df.to_string()}"

    @property
    def difference_df(self) -> pd.DataFrame:
        indices = (
            np.arange(self.timebase_length)
            if config().testing.VERBOSE_OUTPUT
            else self.anomalies.flatten()
        )
        anomaly = self.anomalies[indices]
        mds_data = self.mdsplus_column_data.iloc[indices]
        sql_data = self.sql_column_data.iloc[indices]
        return pd.DataFrame(
            {
                "MDS Time": self.mds_time[indices],
                "SQL Time": self.sql_time[indices],
                "MDSplus Data": mds_data,
                "Reference Data (SQL)": sql_data,
                "MDS/SQL": mds_data / sql_data,
                "Absolute difference": abs(mds_data - sql_data),
                "Relative difference": self.relative_difference[indices],
                "Anomaly": anomaly,
            }
        )

    def compute_numeric_anomalies(self):
        """
        Get the indices of the data where numeric differences exist between the sql and mdsplus data.
        """

        # handle missing data case
        if self.missing_mdsplus_data or self.missing_sql_data:
            if self.missing_mdsplus_data and self.missing_sql_data:
                missing_timebase_length = 0
            elif self.missing_mdsplus_data:
                missing_timebase_length = len(self.sql_column_data)
            else:
                missing_timebase_length = len(self.mdsplus_column_data)
            return np.ones(missing_timebase_length, dtype=bool), np.zeros(
                missing_timebase_length
            )

        sql_is_nan = pd.isnull(self.sql_column_data)
        mdsplus_is_nan = pd.isnull(self.mdsplus_column_data)

        # handle case where both arrays are all null
        if sql_is_nan.all() and mdsplus_is_nan.all():
            return np.zeros(len(self.mdsplus_column_data), dtype=bool), np.zeros(
                len(self.mdsplus_column_data)
            )

        relative_difference = safe_cast(
            np.where(
                self.sql_column_data != 0,
                np.abs(
                    (self.mdsplus_column_data - self.sql_column_data)
                    / self.sql_column_data
                ),
                np.where(self.mdsplus_column_data != 0, np.inf, np.nan),
            ),
            "float64",
        )  # necessary in case all produced values are nan

        numeric_anomalies_mask = np.where(
            np.isnan(relative_difference),
            False,
            relative_difference > config().testing.VAL_TOLERANCE,
        )
        nan_anomalies_mask = sql_is_nan != mdsplus_is_nan
        anomalies: pd.Series = numeric_anomalies_mask | nan_anomalies_mask

        return anomalies.to_numpy(), relative_difference
