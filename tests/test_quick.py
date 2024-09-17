#!/usr/bin/env python3

"""
execute a few simple workflows as tests.
"""

from examples.efit import main as test_efit
from examples.mdsplus import main as test_mdsplus
from examples.sql import main as test_sql

__all__ = ["test_efit", "test_mdsplus", "test_sql"]
