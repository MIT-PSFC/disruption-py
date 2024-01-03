import traceback
import logging
from disruption_py.mdsplus_integration.tree_manager import EnvModifications
from disruption_py.settings.enum_options import SignalDomain
from disruption_py.settings.set_times_request import SetTimesRequest
from disruption_py.shots.cmod.cmod_data_requests import BasicCmodRequests 

from disruption_py.shots.shot import Shot, ShotSetupParams
from disruption_py.shots.shot_data_request import ShotDataRequestParams
from disruption_py.utils.method_caching import parameter_cached_method, cached_method
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.utils import without_duplicates

import pandas as pd
import numpy as np

import MDSplus
from MDSplus import *


class CModShot(Shot):
    """
    Class for a single CMod shot.
    """

    def __init__(
        self, 
        shot_id : str,
        num_threads_per_shot : int,
        override_exising_data : bool,
        set_times_request : SetTimesRequest,
        signal_domain : SignalDomain,
        existing_data : pd.DataFrame,
        disruption_time : float,
        efit_tree_name : str,
        attempt_local_efit_env : EnvModifications,
        **kwargs
    ):
        efit_names_to_test = without_duplicates([
            efit_tree_name,
            "analysis",
            *[f"efit0{i}" for i in range(1, 10)],
            *[f"efit{i}" for i in range(10, 19)],
        ])
        tree_nicknames = { "efit" : (efit_names_to_test, attempt_local_efit_env) }
        
        shot_setup_params = ShotSetupParams(
            shot_id=shot_id, 
            tokamak=Tokamak.CMOD, 
            num_threads_per_shot=num_threads_per_shot, 
            override_exising_data=override_exising_data, 
            set_times_request=set_times_request, 
            signal_domain=signal_domain,
            existing_data=existing_data,
            disruption_time=disruption_time,
            tree_nicknames=tree_nicknames,
        )
        super().__init__(setup_params=shot_setup_params)

    def get_shot_data_request_params(self):
        return ShotDataRequestParams(self, self.logger, Tokamak.CMOD)

    def set_flattop_timebase(self):
        shot_data_requests_params = self.get_shot_data_request_params()
        ip_parameters = BasicCmodRequests._get_ip_parameters(params=shot_data_requests_params)
        ipprog, dipprog_dt = ip_parameters['ip_prog'], ip_parameters['dipprog_dt']
        indices_flattop_1 = np.where(np.abs(dipprog_dt) <= 1e3)[0]
        indices_flattop_2 = np.where(np.abs(ipprog) > 1.e5)[0]
        indices_flattop = np.intersect1d(indices_flattop_1, indices_flattop_2)
        if len(indices_flattop) == 0:
            self.logger.warning(
                f"[Shot {self._shot_id}]:Could not find flattop timebase. Defaulting to full shot(efit) timebase.")
            return
        self._times = self._times[indices_flattop]

    def set_rampup_and_flattop_timebase(self):
        shot_data_requests_params = self.get_shot_data_request_params()
        ip_parameters = BasicCmodRequests._get_ip_parameters(params=shot_data_requests_params)
        ipprog, dipprog_dt = ip_parameters['ip_prog'], ip_parameters['dipprog_dt']
        indices_flattop_1 = np.where(np.abs(dipprog_dt) <= 6e4)[0]
        indices_flattop_2 = np.where(np.abs(ipprog) > 1.e5)[0]
        indices_flattop = np.intersect1d(indices_flattop_1, indices_flattop_2)
        if len(indices_flattop) == 0:
            self.logger.warning(
                f"[Shot {self._shot_id}]:Could not find flattop timebase. Defaulting to full shot(efit) timebase.")
            return
        end_index = np.max(indices_flattop)
        self._times = self._times[:end_index]

