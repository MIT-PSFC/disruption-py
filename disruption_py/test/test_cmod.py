"""Unit tests for workflows involving get_dataset_df() for obtaining CMOD data.

Expects to be run on the MFE workstations.
"""
import unittest

from disruption_py.ml.preprocessing import get_dataset_df

class TestSQL(unittest.TestCase):
    """Tests for the using the SQL table as a data source."""

    def test_get_all_columns(self):
        """Ensure there is no error when getting all columns from the SQL table for several shots.
        """
        self.assertTrue(False)



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