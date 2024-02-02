import numpy as np
from dataclasses import replace

import pandas as pd
from disruption_py.mdsplus_integration.tree_manager import TreeManager
from disruption_py.settings.shot_data_request import ShotDataRequestParams
from disruption_py.settings.shot_settings import ShotSettings
from disruption_py.shots.parameter_methods.cmod.basic_parameter_methods import BasicCmodRequests
from disruption_py.shots.shot_manager import ShotManager
from disruption_py.shots.shot_props import ShotProps
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.utils import without_duplicates


class CModShotManager(ShotManager):
    
    @classmethod
    def cmod_setup_shot_props(
        cls,
        shot_id : int,
        existing_data : pd.DataFrame,
        disruption_time : float,
        shot_settings : ShotSettings,
        **kwargs
    ) -> ShotProps:
        """
        Sets up the shot properties for cmod.
        """
        
        efit_names_to_test = without_duplicates([
            shot_settings.efit_tree_name,
            "analysis",
            *[f"efit0{i}" for i in range(1, 10)],
            *[f"efit{i}" for i in range(10, 19)],
        ])
        efit_envs_to_test = [shot_settings.attempt_local_efit_env, ()] if shot_settings.attempt_local_efit_env is not None else [()]
        tree_nicknames = { "efit_tree" : (efit_names_to_test, efit_envs_to_test) }
            
        tree_manager = TreeManager(shot_id)
        try:
            shot_props = cls.setup(
                shot_id=shot_id,
                tree_manager=tree_manager,
                initial_existing_data=existing_data,
                disruption_time=disruption_time,
                tree_nicknames=tree_nicknames,
                shot_settings=shot_settings,
                tokamak=Tokamak.CMOD,
                **kwargs
            )
            return shot_props
        except Exception as e:
            cls.logger.info(f"[Shot {shot_id}]: Caught failed to setup shot {shot_id}, cleaning up tree manager.")
            tree_manager.cleanup()
            raise e
    
    @classmethod
    def _modify_times_flattop_timebase(cls, shot_props : ShotProps, **kwargs):
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
    def _modify_times_rampup_and_flattop_timebase(cls, shot_props : ShotProps, **kwargs):
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