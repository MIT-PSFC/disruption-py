"""Unit tests for workflows involving get_dataset_df() for obtaining CMOD data.

Expects to be run on the MFE workstations.
Expects MDSplus to be installed and configured.
Expects SQL credentials to be configured.
"""
import unittest

from disruption_py.ml.preprocessing import get_dataset_df

# Shot list used for testing
# Mix of disruptive and non-disruptive shots present in SQL and MDSplus
TEST_SHOTS = [1150805012,   # Flattop Disruption
            1150805013,     # No Disruption
            1150805014,     # No Disruption
            1150805015,     # Rampdown Disruption
            1150805016,     # Rampdown Disruption
            1150805017,     # Rampdown Disruption
            1150805019,     # Rampdown Disruption
            1150805020,     # Rampdown Disruption
            1150805021,     # Rampdown Disruption
            1150805022]     # Flattop Disruption

TIME_EPSILON = 0.05 # Tolerance for taking the difference between two times [s]
IP_EPSILON = 1e5    # Tolerance for taking the difference between two ip values [A]

MATCH_FRACTION = 0.95   # Fraction of signals that must match between MDSplus and SQL

class TestSQL(unittest.TestCase):
    """Tests for the using the SQL table as a data source."""

    def test_get_all_columns(self):
        """Ensure there is no error when getting all columns from the SQL table for several shots.
        """
        sql_df = get_dataset_df(tokamak='cmod', shot_ids=TEST_SHOTS, data_source=2)

        self.assertIsNotNone(sql_df)



class TestMDSplus(unittest.TestCase):
    """Tests for the using MDSplus as a data source. 
    SQL dataset is considered 'blessed' and is used as the reference."""

    # Tests for matching signals between MDSplus and SQL

    def test_ip(self):
        """Ensure that the ip matches the SQL dataset
        """
        for shot in TEST_SHOTS[:]:
            sql_df = get_dataset_df(tokamak='cmod', shot_ids=[shot], data_source=2, cols=['time', 'ip'])
            mds_df = get_dataset_df(tokamak='cmod', shot_ids=[shot], data_source=3, cols=['time', 'ip'])
            
            num_match = 0
            total_count = 0
            # For each time in the SQL dataset, find the closest time in the MDSplus dataset
            # and compare the ip values
            for sql_time, sql_ip in zip(sql_df['time'], sql_df['ip']):
                mds_time = mds_df['time'].iloc[(mds_df['time'] - sql_time).abs().argsort()[:1]]
                mds_ip = mds_df['ip'].loc[mds_time.index]
                if abs(sql_ip - mds_ip.values[0]) < IP_EPSILON:
                    num_match += 1
                total_count += 1

            
            self.assertGreater(num_match/total_count, MATCH_FRACTION)
    

        

    def test_dipprog_dt(self):
        """Ensure that the dipprog_dt matches the SQL dataset
        """
        for shot in TEST_SHOTS[:]:
            sql_df = get_dataset_df(tokamak='cmod', shot_ids=[shot], data_source=2, cols=['time', 'dipprog_dt'])
            mds_df = get_dataset_df(tokamak='cmod', shot_ids=[shot], data_source=3, cols=['time', 'dipprog_dt'])
            
            num_match = 0
            total_count = 0
            # For each time in the SQL dataset, find the closest time in the MDSplus dataset
            # and compare the ip values
            for sql_time, sql_dipprog_dt in zip(sql_df['time'], sql_df['dipprog_dt']):
                mds_time = mds_df['time'].iloc[(mds_df['time'] - sql_time).abs().argsort()[:1]]
                mds_dipprog_dt = mds_df['dipprog_dt'].loc[mds_time.index]
                if abs(sql_dipprog_dt - mds_dipprog_dt.values[0]) < IP_EPSILON:
                    num_match += 1
                total_count += 1

            
            self.assertGreater(num_match/total_count, MATCH_FRACTION)
    
    # Other tests for MDSplus

    def test_flattop_times(self):
        """Ensure that the flattop time matches the change in programmed ip in the SQL dataset.
        """
        for shot in TEST_SHOTS[:]:
            sql_df = get_dataset_df(tokamak='cmod', shot_ids=[shot], data_source=2, cols=['time', 'dipprog_dt'])
            # Find the first time where dipprog_dt is zero
            sql_flattop_time = sql_df['time'].loc[sql_df['dipprog_dt'] == 0].iloc[0]
            # Find the last time where dipprog_dt is zero
            sql_flattop_end_time = sql_df['time'].loc[sql_df['dipprog_dt'] == 0].iloc[-1]

            flattop_df = get_dataset_df(tokamak='cmod', shot_ids=[shot], timebase_signal='flattop', data_source=3, cols=['time']) 
            
            # Find the first time in the flattop signal
            mds_flattop_time = flattop_df['time'].iloc[0]
            # Find the last time in the flattop signal
            mds_flattop_end_time = flattop_df['time'].iloc[-1]

            self.assertAlmostEqual(mds_flattop_time, sql_flattop_time, delta=TIME_EPSILON)
            self.assertAlmostEqual(mds_flattop_end_time, sql_flattop_end_time, delta=TIME_EPSILON)
            