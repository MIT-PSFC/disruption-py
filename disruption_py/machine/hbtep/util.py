#!/usr/bin/env python3

"""
Module for helper, not physics, methods.
"""

import numpy as np


class HbtepUtilMethods:
    """
    A class of helper methods that might fetch and compute data from MDSplus
    but are not physics methods.
    """

    @staticmethod
    def unwrap_phase(phase: np.ndarray):
        """
        Takes in phase array (in radians).  I think it needs to be centered about 0.
        Unwraps phase data so that it is continuous.
        This is important for phase data when you want to take it's derivative to
        get frequency.

        Parameters
        ----------
        phase : numpy.ndarray
            data being unwrapped

        Returns
        ------
        phase_unwrapped : numpy.ndarray
            unwrapped data array

        """
        phase_unwrapped = np.zeros(len(phase))
        offset = 0
        phase_unwrapped[0] = phase[0]
        for i in range(1, len(phase)):
            if phase[i - 1] > np.pi / 4 and phase[i] < -np.pi / 4:
                offset += 2 * np.pi
            elif phase[i - 1] < -np.pi / 4 and phase[i] > np.pi / 4:
                offset -= 2 * np.pi
            phase_unwrapped[i] = phase[i] + offset
        return phase_unwrapped
