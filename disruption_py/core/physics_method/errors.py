#!/usr/bin/env python3

"""
Custom exceptions for the physics methods.
"""


class CalculationError(Exception):
    """
    Custom error specific to physics methods that should be raised when we know
    the result of a calculation is invalid.
    """

    def __init__(self, message):
        super().__init__(message)
