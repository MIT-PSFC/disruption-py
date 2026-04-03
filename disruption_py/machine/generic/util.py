#!/usr/bin/env python3

"""
Module for helper, not physics, methods.
"""

import numpy as np
import scipy


class GenericUtilMethods:
    """
    A class of helper methods that might fetch and compute data from MDSplus
    but are not physics methods.
    """

    @staticmethod
    def get_end_of_current(ip, ip_time, threshold=1e5):
        """
        Python translation of the original MATLAB implementation based on end_of_current_d3d.m

        Returns a dictionary of `duration` and `polarity`.
        """
        # Determine if there was any finite plasma current on this shot.  If
        # not, then the shot was a "no plasma" shot, and the end-of-shot is set
        # to 0 s.
        (finite_indices,) = np.where((ip_time >= 0) & (abs(ip) > threshold))
        if len(finite_indices) == 0:
            return {"duration": 0, "ip_max": 0, "polarity": 1}
        dt = np.diff(ip_time)
        duration = sum(dt[finite_indices[:-1]])
        if duration < 0.1:
            # Assume < 100 ms is not a bona fide plasma
            return {"duration": 0, "ip_max": 0, "polarity": 1}

        # Determine the polarity of the plasma current.
        polarity = np.sign(
            scipy.integrate.trapezoid(ip[finite_indices], ip_time[finite_indices])
        )
        ip_upright = ip * polarity

        # Find all the times that Ip is greater than the threshold.  The largest
        # time value is the end of current.  But also check to see if the plasma
        # current signal has been digitized for long enough to capture the end of
        # the discharge.  If not, negate the end of shot value to indicate that
        # that actual value cannot be determined, but it is longer than
        # abs(value).
        (indices,) = np.where((ip_upright >= threshold) & (ip_time > 0))
        # Get the last index
        max_idx = indices[-1] if len(indices) > 0 else None
        duration = ip_time[max_idx]
        if max_idx == len(ip_time) - 1:
            duration = -duration  # TODO: what is this for?

        # Find Ip_max (with correct polarity)
        ip_max = max(ip_upright) * polarity

        return {"duration": duration, "polarity": polarity}
