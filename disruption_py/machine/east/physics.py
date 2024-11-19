#!/usr/bin/env python3

"""
Module for retrieving and calculating data for DIII-D physics methods.
"""

import traceback

import numpy as np
import scipy
from MDSplus import mdsExceptions

from disruption_py.core.physics_method.caching import cache_method
from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.errors import CalculationError
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import (
    interp1,
    matlab_get_bolo,
    matlab_gsastd,
    matlab_power,
)
from disruption_py.machine.tokamak import Tokamak



class EASTPhysicsMethods:
    """
    A class to retrieve and calculate physics-related data for DIII-D.
    """
    
    @staticmethod
    @physics_method(columns=["time_until_disrupt"], tokamak=Tokamak.EAST)
    def get_time_until_disrupt(params: PhysicsMethodParams):
        """
        Calculate the time until the disruption for a given shot. If the shot does
        not disrupt, return NaN.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information.

        Returns
        -------
        dict
            A dictionary containing the time until disruption. If the shot does
            not disrupt, return NaN.
        """
        if params.disrupted:
            return {"time_until_disrupt": params.disruption_time - params.times}
        return {"time_until_disrupt": [np.nan]}

    @staticmethod
    @physics_method(
        columns=["ip", "ip_error", "dip_dt", "dipprog_dt", "power_supply_railed"],
        tokamak=Tokamak.D3D,
    )
    def get_ip_parameters(params: PhysicsMethodParams):
        """
        Retrieve plasma current parameters including measured and programmed values.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the following keys:
            - 'ip' : array
                Measured plasma current [A].
            - 'ip_prog' : array
                Programmed (requested) plasma current [A].
            - 'ip_error' : array
                ip - ip_prog [A].
            - 'dip_dt' : array
                Time derivative of the measured plasma current [A/s].
            - 'dipprog_dt' : array
                Time derivative of the programmed plasma current [A/s]
        
        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/EAST/get_Ip_parameters.m
        
        Original author: Robert Granetz, Dec 2015
        Last major update: 11/19/24 by William Wei
        """
        ip = [np.nan]
        ip_prog = [np.nan]
        dip_dt = [np.nan]
        dipprog_dt = [np.nan]
        # Fill with nans instead of using a single nan because indices are used
        ip_error = np.full(len(params.times), np.nan)
        
        # Read in the measured plasma current, Ip.  There are several
        # different measurements of Ip: IPE, IPG, IPM (all in the EAST tree), and
        # PCRL01 (in the PCS_EAST tree).  At various times in the history of EAST,
        # there are been problems with all of these measurements, such as broken
        # sensors, inverted signals, and shifted timebases.  I think the most
        # reliable one is PCRL01, which is the one used by the Plasma Control
        # System (PCS) for feedback control.  So that is the one I will use for the
        # disruption warning database.
        ip, ip_time = params.mds_conn.get_data_with_dims(
            r"\pcrl01", tree_name="pcs_east"
        )  # [A], [s]
        
        # For shots before year 2014, the PCRL01 timebase needs to be shifted 
        # by 17.0 ms 
        if params.shot_id < 44432:
            ip_time -= 0.0170
            
        # High-frequency noise spikes on some shots can cause a problem with the
        # time derivative and other computations.  Use a median filter to reduce
        # the problem.
        ip = scipy.signal.medfilt(ip, 5)    # Remove noise spikes with median filter

        # Subtract baseline offset
        (base_indices,) = np.where(ip_time <= -5.8)    # time before any PF supplies turn on
        if len(base_indices) > 0:
            baseline = sum(ip[base_indices]) / len(base_indices)
            ip -= baseline

        # Calculate dip_dt
        dip_dt = np.gradient(ip, ip_time)
        
        
        
        # Interpolate all retrieved signals to the requested timebase
        # TODO

        output = {
            'ip': ip,
            'ip_prog': ip_prog,
            'ip_error': ip_error,
            'dip_dt': dip_dt,
            'dipprog_dt': dipprog_dt,
        }
        return output