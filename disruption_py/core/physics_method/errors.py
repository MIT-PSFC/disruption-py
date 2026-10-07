#!/usr/bin/env python3

"""
Custom exceptions for the physics methods.
"""


class CustomError(Exception):
    """
    Base custom error.
    """


class DataError(CustomError):
    """
    Generic data error.
    """


class CalculationError(CustomError):
    """
    Generic calculation error.
    """


class FetchDataError(DataError):
    """
    Fetch data error.
    """


class NanDataError(DataError):
    """
    NaN data error.
    """


class MismatchCalculationError(CalculationError):
    """
    Mismatch calculation error.
    """
