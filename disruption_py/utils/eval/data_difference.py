from dataclasses import dataclass, field
import inspect
from disruption_py.utils.constants import MATCH_FRACTION, VAL_TOLERANCE, VERBOSE_OUTPUT
from typing import Callable, Dict, List
import numpy as np
import pandas as pd

@dataclass
class DataDifference:
    """
    Data difference between MDSplus and SQL.
    """
    shot_id : int
    data_column : str
    anomalies : np.ndarray = field(init=False) # 1 if anomaly, 0 o.w.
    relative_difference : np.ndarray = field(init=False)
    mdsplus_column_data : pd.Series
    sql_column_data : pd.Series
    expect_failure : bool
    
    def __post_init__(self):
        self.anomalies, self.relative_difference = self.compute_numeric_anomalies()
        
    @property
    def num_anomalies(self) -> int:
        return np.sum(self.anomalies)
    
    @property
    def timebase_length(self) -> int:
        return len(self.anomalies)
    
    @property
    def failed(self) -> str:
        return self.num_anomalies / self.timebase_length > 1 - MATCH_FRACTION
    
    @property
    def failure_ratio_string(self) -> str:
        return f"{self.num_anomalies / self.timebase_length:.4f}"
    
    @property
    def column_mismatch_string(self) -> str:
        return f"Shot {self.shot_id} column {self.data_column} with arrays:\n{self.difference_df.to_string()}"
    
    @property
    def difference_df(self) -> pd.DataFrame:
        indexes = np.arange(self.timebase_length) if VERBOSE_OUTPUT else self.anomalies.flatten()
        anomaly = self.anomalies[indexes]
        return pd.DataFrame({
            'MDSplus Data': self.mdsplus_column_data.iloc[indexes], 
            'Reference Data (SQL)': self.sql_column_data.iloc[indexes], 
            'Relative difference': self.relative_difference[indexes], 
            'Anomaly': anomaly
        })
    
    @staticmethod
    def test_shots(
        shot_ids : List[int],
        mdsplus_data : Dict[int, pd.DataFrame],
        sql_data : Dict[int, pd.DataFrame],
        data_columns : List[str],
        fail_quick : bool = False,
        expected_failure_columns : List[str] = None,
    ) -> List["DataDifference"]:
        """
        Test if the difference between the two data is within tolerance.
        """
        if expected_failure_columns is None:
            expected_failure_columns = []
        
        data_differences : List[DataDifference] = []
        for data_column in data_columns:
            for shot_id in shot_ids:
                mdsplus_shot_data, sql_shot_data = mdsplus_data[shot_id], sql_data[shot_id]
                data_difference = DataDifference.test_shot(
                    shot_id=shot_id, 
                    mdsplus_shot_data=mdsplus_shot_data, 
                    sql_shot_data=sql_shot_data, 
                    data_column=data_column,
                    fail_quick = fail_quick,
                    expect_failure = data_column in expected_failure_columns
                )
                data_differences.append(data_difference)
        return data_differences
           
    @staticmethod
    def test_shot(
        shot_id : int,
        mdsplus_shot_data : pd.DataFrame,
        sql_shot_data : pd.DataFrame,
        data_column : str,
        fail_quick : bool = False,
        expect_failure : bool = False,
    ) -> "DataDifference":
        """
        Test if the difference between the two data is within tolerance.
        """
        if data_column not in mdsplus_shot_data:
            raise ValueError(f"Column {data_column} missing from MDSPlus for shot {shot_id}")

        if data_column not in sql_shot_data:
            raise ValueError(f"Column {data_column} missing from SQL for shot {shot_id}")
            
        data_difference = DataDifference(
            shot_id=shot_id,
            data_column=data_column,
            mdsplus_column_data=mdsplus_shot_data[data_column],
            sql_column_data=sql_shot_data[data_column],
            expect_failure=expect_failure,
        )
        
        if fail_quick:
            if expect_failure:
                assert data_difference.failed, "Expected failure but succeeded:\n{}".format(data_difference.column_mismatch_string)
            else:
                assert not data_difference.failed, "Expected success but failed:\n{}".format(data_difference.column_mismatch_string)
        
        return data_difference

    @staticmethod
    def get_failure_statistics_string(data_differences : list["DataDifference"], data_column=None):
        data_difference_by_column = {}
        for data_difference in data_differences:
            data_difference_by_column.setdefault(data_difference.data_column, []).append(data_difference)
        
        failure_strings = {}
        failed_columns, succeeded_columns = set(), set()
        for ratio_data_column, data_differences in data_difference_by_column.items():
            failures = [data_difference.shot_id for data_difference in data_differences if data_difference.failed]
            failed = len(failures) > 0
            
            anomaly_count = sum([data_difference.num_anomalies for data_difference in data_differences])
            timebase_count = sum([data_difference.timebase_length for data_difference in data_differences])
            
            matches_expected_failures = all([data_difference.expect_failure == data_difference.failed for data_difference in data_differences])
            
            
            conditions : Dict[str, Callable[[DataDifference], bool]] = {
                "Shots expected to fail that failed": lambda data_difference: data_difference.expect_failure and data_difference.failed,
                "Shots expected to succeed that failed": lambda data_difference: not data_difference.expect_failure and data_difference.failed,
                "Shots expected to fail that succeeded": lambda data_difference: data_difference.expect_failure and not data_difference.failed,
                "Shots expected to succeed that succeeded": lambda data_difference: not data_difference.expect_failure and not data_difference.failed,
            }
            condition_results = {}
            for condition_name, condition in conditions.items():
                shot_ids = [data_difference.shot_id for data_difference in data_differences if condition(data_difference)]
                if len(shot_ids) > 0:
                    condition_results[condition_name] = shot_ids
            condition_string = "\n".join([f"{condition_name} ({len(condition_result)} shots): {condition_result}" for condition_name, condition_result in condition_results.items()])
            
            
            failure_string = f"""\
            Column {ratio_data_column} {"FAILED" if failed else "SUCCEEDED"}
            Matches expected failures: {matches_expected_failures}
            {condition_string}
            Total Entry Failure Rate: {anomaly_count / timebase_count * 100:.2f}%
            """
            failure_strings[ratio_data_column] = inspect.cleandoc(failure_string)
            
            if failed:
                failed_columns.add(ratio_data_column)
            else:
                succeeded_columns.add(ratio_data_column)
        
        if data_column is not None:
            return failure_strings.get(data_column, "")
        else:
            summary_string = f"""\
            SUMMARY
            Columns with a failure:
            {", ".join(failed_columns)}
            
            Columns without a failure :
            {", ".join(succeeded_columns)}
            """
            return '\n\n'.join(failure_strings.values()) + '\n\n' + inspect.cleandoc(summary_string)
                
    def compute_numeric_anomalies(self):
        """
        Get the indices of the data where numeric differences exist between the sql and mdsplus data.
        """
        
        # handle case where both arrays are all null
        if pd.isnull(self.sql_column_data).all() and pd.isnull(self.mdsplus_column_data).all():
            return np.zeros(len(self.mdsplus_column_data), dtype=bool), np.zeros(len(self.mdsplus_column_data))
        
        relative_difference = np.where(
            self.sql_column_data != 0, 
            np.abs((self.mdsplus_column_data - self.sql_column_data) / self.sql_column_data), 
            np.where(self.mdsplus_column_data != 0, np.inf, np.nan)
        )
        
        numeric_anomalies_mask = (relative_difference > VAL_TOLERANCE)
        
        sql_is_nan_ = pd.isnull(self.sql_column_data)
        mdsplus_is_nan = pd.isnull(self.mdsplus_column_data)
        nan_anomalies_mask = (sql_is_nan_ != mdsplus_is_nan)
        
        anomalies : pd.Series = numeric_anomalies_mask | nan_anomalies_mask
        
        return anomalies.to_numpy(), relative_difference 