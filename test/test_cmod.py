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
        self.assertTrue(False)

    def test_dipprog_dt(self):
        """Ensure that the dipprog_dt matches the SQL dataset
        """

        self.assertTrue(False)
    
    # Other tests for MDSplus

    def test_flattop_times(self):
        """Ensure that the flattop time matches the change in programmed ip in the SQL dataset.
        """
        self.assertTrue(False)