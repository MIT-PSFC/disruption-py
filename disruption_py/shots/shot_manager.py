from abc import ABC, abstractmethod
from typing import Dict, Tuple, List
import pandas as pd
import numpy as np
import logging
import traceback
from disruption_py.mdsplus_integration.tree_manager import TreeManager, EnvModifications
from disruption_py.settings.enum_options import InterpolationMethod, SignalDomain
from disruption_py.settings.set_times_request import SetTimesRequest, SetTimesRequestParams
from disruption_py.settings.shot_data_request import ShotDataRequestParams
from disruption_py.settings.shot_settings import ShotSettings
from disruption_py.shots.helpers.populate_shot import populate_shot
from disruption_py.shots.shot_props import ShotProps
from disruption_py.utils.command_utils import get_commit_hash
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.math_utils import interp1
from disruption_py.utils.utils import without_duplicates

MAX_SHOT_TIME = 7.0  # [s]

class ShotManager(ABC):
    logger = logging.getLogger('disruption_py')
    
    @classmethod
    @abstractmethod
    def _modify_times_flattop_timebase(cls, shot_props : ShotProps, **kwargs) -> ShotProps:
        pass
    
    @classmethod
    @abstractmethod
    def _modify_times_rampup_and_flattop_timebase(cls, shot_props : ShotProps, **kwargs) -> ShotProps:
        pass
    
    @classmethod
    def setup(
        cls,
        shot_id : str,
        tokamak: Tokamak,
        existing_data : pd.DataFrame,
        disruption_time : float,
        tree_nicknames : Dict[str, Tuple[List[str], List[EnvModifications]]],
        shot_settings : ShotSettings,
        **kwargs
    ) -> ShotProps:
    
        tree_manager = TreeManager(shot_id)
        metadata = {
            'labels': {},
            'commit_hash': get_commit_hash(),
            'timestep': {},
            'duration': {},
            'description': "",
            'disrupted': 100  # TODO: Fix
        }
        
        cls._init_nicknames(tree_manager, tree_nicknames)
        
        interpolation_method  = interp1 # TODO: fix
        
        times = cls._init_times(
            shot_id=shot_id, 
            existing_data=existing_data, 
            tree_manager=tree_manager, 
            tokamak=tokamak, 
            disruption_time=disruption_time,
            shot_settings=shot_settings
        )
        
        populated_existing_data = cls._init_data(
            times=times,
            existing_data=existing_data,
        )
        
        shot_props = ShotProps(
            shot_id=shot_id,
            tokemak=tokamak,
            num_threads_per_shot=shot_settings.num_threads_per_shot,
            disruption_time = disruption_time,
            tree_manager = tree_manager,
            initial_existing_data = existing_data,
            populated_existing_data = populated_existing_data,
            interpolation_method = interpolation_method,
            metadata = metadata,
        )
        
        # modify already existing shot props, such as modifying timebase
        shot_props = cls._modify_shot_props_for_settings(
            shot_props, 
            shot_settings, 
            **kwargs
        )
        
        return shot_props
    
    @classmethod
    def run_data_retrieval(cls, shot_props : ShotProps, shot_settings : ShotSettings):
        shot_data_request_params = ShotDataRequestParams(shot_props, cls.logger, shot_props.tokamak)
        return populate_shot(shot_settings=shot_settings, params=shot_data_request_params)
    
    @classmethod
    def cleanup(cls, shot_props : ShotProps,):
        shot_props.cleanup()
    
    @classmethod
    def _modify_shot_props_for_settings(
        cls, 
        shot_props : ShotProps,
        shot_settings : ShotSettings,
        **kwargs
    ) -> ShotProps:
        if shot_settings.signal_domain is SignalDomain.FLATTOP:
            shot_props = cls._modify_times_flattop_timebase(shot_props)
        elif shot_settings.signal_domain is SignalDomain.RAMP_UP_AND_FLATTOP:
            shot_props = cls._modify_times_rampup_and_flattop_timebase(shot_props)
        
        if shot_props is None:
            cls.logger.error(f"Shot_props set to None in modify_shot_props()")
        
        return shot_props
    
    @classmethod
    def _init_nicknames(
        cls, 
        tree_manager : TreeManager, 
        tree_nicknames : Dict[str, Tuple[List[str], List[EnvModifications]]],
    ):
        for nickname, (tree_names, env_modifications) in tree_nicknames.items():
            efit_names_to_test = without_duplicates(tree_names)
            tree_manager.nickname(nickname, efit_names_to_test, env_modifications)
    
    @classmethod
    def _init_times(
        cls,
        shot_id : str,
        existing_data : pd.DataFrame, 
        tree_manager : TreeManager,
        tokamak : Tokamak,
        disruption_time : float,
        shot_settings : ShotSettings,
    ):
        """
        Initialize the timebase of the shot.
        """
        if existing_data is not None and shot_settings.override_exising_data is False:
            # set timebase to be the timebase of existing data
            try:
                times = existing_data['time'].to_numpy()
                # Check if the timebase is in ms instead of s
                if times[-1] > MAX_SHOT_TIME:
                    times /= 1000  # [ms] -> [s]
                return times
            except KeyError as e:
                cls.logger.warning(
                    f"[Shot {shot_id}]: Shot constructor was passed data but no timebase.")
                cls.logger.debug(
                    f"[Shot {shot_id}]:{traceback.format_exc()}")
        else:
            request_params = SetTimesRequestParams(tree_manager=tree_manager, tokamak=tokamak, logger=cls.logger, disruption_time=disruption_time)
            return shot_settings.set_times_request.get_times(request_params)
    
    @classmethod
    def _init_data(times : np.ndarray, existing_data : pd.DataFrame):
        '''
        Intialize the shot with data, if existing data matches the shot timebase.
        '''
        if existing_data is not None:
            time_df = pd.DataFrame(times, columns=['time'])
            flagged_existing_data = existing_data.assign(merge_success_flag=1)
            timed_existing_data = pd.merge_asof(time_df, flagged_existing_data, on='time', direction='nearest', tolerance=TIME_CONST)
            if timed_existing_data['merge_success_flag'].isna().any():
                populated_existing_data = None
            else:
                populated_existing_data = timed_existing_data.drop(columns=['merge_success_flag'])
        
        return populated_existing_data