"""
Module for helper, not physics, methods. 
"""

import numpy as np

from disruption_py.core.physics_method.params import PhysicsMethodParams


class D3DUtilMethods:
    """
    A class of helper methods that might fetch and compute data from MDSPlus
    but are not physics methods.
    """

    @staticmethod
    def get_polarity(params: PhysicsMethodParams):
        """
        Get the plasma current polarity. Accepts PhysicsMethodParams to access
        the MDS connection, but it is not a physics method.

        Returns the first value of polarity array if the polarity is not constant.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information.

        Returns
        -------
        polarity value, -1 or 1.
        """
        polarity = np.unique(
            params.mds_conn.get_data(
                f"ptdata('iptdirect', {params.shot_id})", tree_name="d3d"
            )
        )
        if len(polarity) > 1:
            params.logger.info(
                "Polarity of Ip target is not constant. Using value at first timestep.",
            )
            params.logger.debug("Polarity array: {polarity}", polarity=polarity)
            polarity = polarity[0]
        return polarity
