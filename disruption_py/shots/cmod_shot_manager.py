import numpy as np
from dataclasses import replace
from disruption_py.settings.shot_data_request import ShotDataRequestParams
from disruption_py.shots.parameter_functions.cmod.cmod_data_requests import BasicCmodRequests
from disruption_py.shots.shot_props import ShotProps
from disruption_py.utils.mappings.tokamak import Tokamak


def CModShotManager(ShotManager):
    
    @classmethod
    def modify_times_flattop_timebase(cls, shot_props : ShotProps, **kwargs):
        shot_data_requests_params = ShotDataRequestParams(shot_props=shot_props, logger=cls.logger, tokamak=Tokamak.CMOD)
        ip_parameters = BasicCmodRequests._get_ip_parameters(params=shot_data_requests_params)
        ipprog, dipprog_dt = ip_parameters['ip_prog'], ip_parameters['dipprog_dt']
        indices_flattop_1 = np.where(np.abs(dipprog_dt) <= 1e3)[0]
        indices_flattop_2 = np.where(np.abs(ipprog) > 1.e5)[0]
        indices_flattop = np.intersect1d(indices_flattop_1, indices_flattop_2)
        if len(indices_flattop) == 0:
            cls.logger.warning(
                f"[Shot {shot_props.shot_id}]:Could not find flattop timebase. Defaulting to full shot(efit) timebase.")
            return
        shot_props.times = shot_props.times[indices_flattop]
        shot_props._cached_results.clear() #TODO: Make this only modify the cached results for new times
        return shot_props
        
    @classmethod
    def modify_times_rampup_and_flattop_timebase(cls, shot_props : ShotProps, **kwargs):
        shot_data_requests_params = ShotDataRequestParams(shot_props=shot_props, logger=cls.logger, tokamak=Tokamak.CMOD)
        ip_parameters = BasicCmodRequests._get_ip_parameters(params=shot_data_requests_params)
        ipprog, dipprog_dt = ip_parameters['ip_prog'], ip_parameters['dipprog_dt']
        indices_flattop_1 = np.where(np.abs(dipprog_dt) <= 6e4)[0]
        indices_flattop_2 = np.where(np.abs(ipprog) > 1.e5)[0]
        indices_flattop = np.intersect1d(indices_flattop_1, indices_flattop_2)
        if len(indices_flattop) == 0:
            cls.logger.warning(
                f"[Shot {shot_props.shot_id}]:Could not find flattop timebase. Defaulting to full shot(efit) timebase.")
            return
        end_index = np.max(indices_flattop)
        shot_props.times = shot_props.times[:end_index]
        shot_props._cached_results.clear() #TODO: Make this only modify the cached results for new times
        return shot_props