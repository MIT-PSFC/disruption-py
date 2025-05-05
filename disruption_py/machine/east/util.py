"""
Module for helper, not physics, methods.
"""

import numpy as np
import scipy

from disruption_py.inout.mds import MDSConnection


class EastUtilMethods:
    """
    A class of helper methods that might fetch and compute data from MDSplus
    but are not physics methods.
    """

    @staticmethod
    def subtract_ip_baseline_offset(ip, ip_time):
        """
        Subtract the baseline offset from the plasma current signal.

        Parameters
        ----------
        ip : np.ndarray
            Plasma current [A].
        ip_time : np.ndarray
            Time base of plasma current [s].

        Returns
        -------
        np.ndarray
            Offset-subtracted plasma current [A].
        """
        # Get indices of time before any PF supplies turn on
        (base_indices,) = np.where(ip_time <= -5.8)
        # Subtract offset
        if len(base_indices) > 0:
            baseline = sum(ip[base_indices]) / len(base_indices)
            ip -= baseline
        return ip

    @staticmethod
    def retrieve_ip(mds_conn: MDSConnection, shot_id: int):
        """
        Read in the measured plasma current, Ip. There are several different
        measurements of Ip: IPE, IPG, IPM (all in the EAST tree), and PCRL01
        (in the PCS_EAST tree). At various times in the history of EAST, there
        have been problems with all of these measurements, such as broken sensors,
        inverted signals, and shifted timebases. I think the most reliable one is
        PCRL01, which is the one used by the Plasma Control System (PCS) for feedback
        control. So that is the one I will use for the disruption warning database.

        Parameters
        ----------
        mds_conn : MDSConnection
            Connection to MDSplus server.
        shot_id : int
            Shot number.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            Plasma current [A], time base of plasma current [s].
        """

        ip, ip_time = mds_conn.get_data_with_dims(
            r"\pcrl01", tree_name="pcs_east"
        )  # [A], [s]

        # For shots before year 2014, the PCRL01 timebase needs to be shifted
        # by 17.0 ms
        if shot_id < 44432:
            ip_time -= 0.0170

        # High-frequency noise spikes on some shots can cause a problem with the
        # time derivative and other computations.  Use a median filter to reduce
        # the problem.
        ip = scipy.signal.medfilt(ip, 5)  # Remove noise spikes with median filter

        # Subtract baseline offset
        ip = EastUtilMethods.subtract_ip_baseline_offset(ip, ip_time)
        return ip, ip_time

    @staticmethod
    def get_axuv_calib_factors() -> dict:
        """
        Return the calibration factors for the AXUV arrays

        Used by get_radiated_power and get_prad_peaking.

        - fac_1: (unknown)
        - fac_2: factors of Amp.Gain
        - fac_3: cross calibration factors between arrays
        - fac_4: unit convert
        - fac_5: corrected factor by cross calibration with foil bolometer
        - maj_r: major radius (of the machine cross-section?)
        - del_r: (unknown)
        """
        fac_1 = (
            np.array(
                [
                    1.3681,
                    1.3429,
                    1.3215,
                    1.3039,
                    1.2898,
                    1.2793,
                    1.2723,
                    1.2689,
                    1.2689,
                    1.2723,
                    1.2793,
                    1.2898,
                    1.3039,
                    1.3215,
                    1.3429,
                    1.3681,
                ]
            )
            * 1e4
        )
        fac_2 = np.array(
            [
                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            ]
        )
        fac_3 = np.array([1, 1, 1, 1])
        del_r = (
            np.array(
                [
                    3.6,
                    3.6,
                    3.5,
                    3.4,
                    3.4,
                    3.3,
                    3.3,
                    3.2,
                    3.2,
                    3.1,
                    3.1,
                    3.0,
                    3.0,
                    2.9,
                    2.9,
                    2.8,
                    2.9,
                    2.8,
                    2.8,
                    2.8,
                    2.8,
                    2.8,
                    2.8,
                    2.7,
                    2.7,
                    2.7,
                    2.7,
                    2.7,
                    2.6,
                    2.6,
                    2.6,
                    2.6,
                    2.6,
                    2.6,
                    2.6,
                    2.6,
                    2.7,
                    2.7,
                    2.7,
                    2.7,
                    2.7,
                    2.8,
                    2.8,
                    2.8,
                    2.8,
                    2.8,
                    2.8,
                    2.9,
                    2.8,
                    2.9,
                    2.9,
                    3.0,
                    3.0,
                    3.1,
                    3.1,
                    3.2,
                    3.2,
                    3.3,
                    3.3,
                    3.4,
                    3.4,
                    3.5,
                    3.6,
                    3.6,
                ]
            )
            * 0.01
        )
        # Use original names in the MATLAB script
        return {
            "Fac1": fac_1,
            "Fac2": fac_2,
            "Fac3": fac_3,
            "Fac4": 1e-3,
            "Fac5": 2.5,
            "Maj_R": 1.85,
            "Del_r": del_r,
        }

    @staticmethod
    def load_rmp_saddle_coeff_matrix():
        """
        Load the rmp saddle coefficient matrix.
        Source: rmp_saddle_coeff_matrix.mat
        """
        coeff_matrix = [
            [
                -0.000325933950706964,
                -5.34427695316486e-05,
                -9.60845947944944e-06,
                -5.04741218778060e-06,
                -5.59732219661124e-06,
                -5.91438750088332e-06,
                -7.67820045826030e-06,
                -5.20559984855693e-05,
                -0.000302471162280883,
                -5.35568562149143e-05,
                -9.56803697152602e-06,
                -5.88885912077860e-06,
                -7.99784778858085e-06,
                -6.80525414488422e-06,
                -8.71081837989025e-06,
                -5.76306168766297e-05,
            ],
            [
                -5.22709907162920e-05,
                -0.000320133097404118,
                -5.50071761972670e-05,
                -7.37274061429450e-06,
                -6.90578359290429e-06,
                -5.02997540459485e-06,
                -5.10697946272394e-06,
                -8.99898260398151e-06,
                -5.66360083538718e-05,
                -0.000311386125033681,
                -5.23092310075144e-05,
                -8.27585100002495e-06,
                -8.00023807246884e-06,
                -6.44383893574099e-06,
                -6.54598718347685e-06,
                -1.00724151886658e-05,
            ],
            [
                -8.08294716760164e-06,
                -5.55046044711059e-05,
                -0.000333809266236893,
                -5.60440355539882e-05,
                -9.50511197149052e-06,
                -5.91247287479164e-06,
                -5.20599513677853e-06,
                -6.16780683251725e-06,
                -1.15430117601721e-05,
                -5.51203237240259e-05,
                -0.000316394172222160,
                -5.48875912260367e-05,
                -1.06537242928572e-05,
                -6.47369008200614e-06,
                -6.27085788458336e-06,
                -6.55919579400570e-06,
            ],
            [
                -5.12806076473295e-06,
                -8.49939605799587e-06,
                -5.92476041414216e-05,
                -0.000328925331027763,
                -5.69684251594807e-05,
                -8.70337525570129e-06,
                -4.92573011225709e-06,
                -7.82605756687085e-06,
                -8.55723442676295e-06,
                -8.33796262503122e-06,
                -5.77746097768847e-05,
                -0.000325606660148107,
                -5.93152260417858e-05,
                -8.66043101910145e-06,
                -6.37549893634100e-06,
                -6.85515607000913e-06,
            ],
            [
                -4.65502958467029e-06,
                -4.85211383504998e-06,
                -1.02244946378116e-05,
                -4.93116559918267e-05,
                -0.000332458254995670,
                -5.51888967708278e-05,
                -7.85611622708372e-06,
                -6.92063279447296e-06,
                -8.36497581142532e-06,
                -5.98112870895933e-06,
                -1.04474088707701e-05,
                -5.44799256210205e-05,
                -0.000326342886258557,
                -5.27273498884329e-05,
                -9.18704874796283e-06,
                -6.78335745990761e-06,
            ],
            [
                -4.24256096566107e-06,
                -4.36018675342360e-06,
                -6.89150573826226e-06,
                -6.63787474197891e-06,
                -5.61987082604549e-05,
                -0.000301238014687208,
                -5.52453316551176e-05,
                -9.13919157975075e-06,
                -8.37854203786939e-06,
                -5.93322139470238e-06,
                -6.58104941723929e-06,
                -8.33865462199638e-06,
                -5.68131341626709e-05,
                -0.000317650462941753,
                -5.63542945554439e-05,
                -9.17221707740544e-06,
            ],
            [
                -7.79176701405456e-06,
                -4.90936998737812e-06,
                -6.98520489721429e-06,
                -4.40637734472976e-06,
                -9.65641938226777e-06,
                -5.33787503140721e-05,
                -0.000318500478130682,
                -5.72159746691176e-05,
                -1.15055865602172e-05,
                -5.31875434926984e-06,
                -6.76325866931821e-06,
                -5.68472690837507e-06,
                -9.76815563125627e-06,
                -5.75939865428040e-05,
                -0.000325379112579523,
                -5.60371003488067e-05,
            ],
            [
                -5.72759424001505e-05,
                -8.13153861612048e-06,
                -7.30722987804662e-06,
                -5.64219477471097e-06,
                -6.40689982720806e-06,
                -1.06255695592301e-05,
                -6.29862623902986e-05,
                -0.000330999771429404,
                -5.87929295332690e-05,
                -7.60812324388172e-06,
                -7.52054522920347e-06,
                -5.45222900075003e-06,
                -8.27237279295883e-06,
                -9.78221302986545e-06,
                -6.81254152056869e-05,
                -0.000336854492397579,
            ],
        ]
        return np.array(coeff_matrix)
