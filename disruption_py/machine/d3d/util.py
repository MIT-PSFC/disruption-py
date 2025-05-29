"""
Module for helper, not physics, methods.
"""

import numpy as np

from disruption_py.core.physics_method.params import PhysicsMethodParams

from MDSplus import mdsExceptions


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

    @staticmethod
    def check_chisq(params: PhysicsMethodParams, signal: np.ndarray) -> np.ndarray:
        """
        Check chisq and remove unreliable data points

        Consider making this a generic method or move to utils (prefer not to because MDS access)
        """
        chisq_threshold = 50  # to be added to machine specific config

        try:
            chisq = params.mds_conn.get_data(
                r"\efit_a_eqdsk:chisq", tree_name="_efit_tree"
            )
            invalid_indices = np.where(chisq > chisq_threshold)
        except mdsExceptions.MdsException as e:
            params.logger.warning(
                "Failed to obtain chisq to remove unreliable time points.",
            )
            params.logger.opt(exception=True).debug(e)
            return signal

        if len(signal) != len(chisq):
            params.logger.warning(
                "Length of signal does not match length of chisq. Return original signal.",
            )

        signal[invalid_indices] = np.nan
        return signal
