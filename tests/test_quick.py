#!/usr/bin/env python3

"""
execute a few simple workflows as tests.
"""


def test_sql():
    """
    test SQL connection.
    """

    import examples.sql

    print("Success:", examples.sql.__file__)


def test_mdsplus():
    """
    test MDSplus connection.
    """

    import examples.mdsplus

    print("Success:", examples.mdsplus.__file__)


def test_efit():
    """
    test EFIT workflow.
    """

    import examples.efit

    print("Success:", examples.efit.__file__)

FRUITS = {"apple": 1, "pear": 5, "peach": 10}
for fruit in FRUITS.keys():
    print(fruit)