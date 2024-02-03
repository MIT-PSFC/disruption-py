from abc import ABC, abstractmethod
from typing import Dict, Tuple, List
import pandas as pd
import numpy as np
import logging
from disruption_py.databases.database import ShotDatabase
from disruption_py.mdsplus_integration.tree_manager import TreeManager, EnvModifications
from disruption_py.settings.enum_options import InterpolationMethod, SignalDomain
from disruption_py.settings.existing_data_request import ExistingDataRequestParams
from disruption_py.settings.set_times_request import SetTimesRequest, SetTimesRequestParams
from disruption_py.settings.shot_data_request import ShotDataRequestParams
from disruption_py.settings.shot_settings import ShotSettings
from disruption_py.shots.helpers.populate_shot import populate_shot
from disruption_py.shots.shot_props import ShotProps
from disruption_py.utils.command_utils import get_commit_hash
from disruption_py.utils.constants import TIME_CONST
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.math_utils import interp1
from disruption_py.utils.utils import without_duplicates

class ShotManager(ABC):
    logger = logging.getLogger('disruption_py')
    
    def __init__(self, database : ShotDatabase, mds_connection_manager):
        self.database = database
        self.mds_connection_manager = mds_connection_manager
        
    @classmethod
    @abstractmethod
    def _modify_times_flattop_timebase(cls, shot_props : ShotProps, **kwargs) -> ShotProps:
        pass
    
    @classmethod
    @abstractmethod
    def _modify_times_rampup_and_flattop_timebase(cls, shot_props : ShotProps, **kwargs) -> ShotProps:
        pass
    
    def shot_setup(
        self,
        shot_id : int,
        tree_manager : TreeManager,
        disruption_time : float,
        tree_nicknames : Dict[str, Tuple[List[str], List[EnvModifications]]],
        shot_settings : ShotSettings,
        tokamak: Tokamak,
        **kwargs
    ) -> ShotProps:
        
        
        self._init_nicknames(tree_manager, tree_nicknames)
        
        existing_data = self._retrieve_existing_data(
            shot_id=shot_id,
            tokamak=tokamak,
            shot_settings=shot_settings,
        )

        interpolation_method  = interp1 # TODO: fix
        
        times = self._init_times(
            shot_id=shot_id, 
            existing_data=existing_data, 
            tree_manager=tree_manager, 
            tokamak=tokamak, 
            disruption_time=disruption_time,
            shot_settings=shot_settings
        )
        
        pre_filled_shot_data = self._pre_fill_shot_data(
            times=times,
            existing_data=existing_data,
        )
        
        metadata = {
            'labels': {},
            'commit_hash': get_commit_hash(),
            'timestep': {},
            'duration': {},
            'description': "",
            'disrupted': 100  # TODO: Fix
        }
        
        shot_props = ShotProps(
            shot_id=shot_id,
            tokamak=tokamak,
            num_threads_per_shot=shot_settings.num_threads_per_shot,
            disruption_time = disruption_time,
            tree_manager = tree_manager,
            times = times,
            existing_data = existing_data,
            pre_filled_shot_data = pre_filled_shot_data,
            interpolation_method = interpolation_method,
            metadata = metadata,
        )
        
        # modify already existing shot props, such as modifying timebase
        shot_props = self._modify_shot_props_for_settings(
            shot_props, 
            shot_settings, 
            **kwargs
        )
        
        return shot_props
    
    @classmethod
    def shot_data_retrieval(cls, shot_props : ShotProps, shot_settings : ShotSettings):
        shot_data_request_params = ShotDataRequestParams(shot_props=shot_props, logger=cls.logger, tokamak=shot_props.tokamak)
        return populate_shot(shot_settings=shot_settings, params=shot_data_request_params)
    
    @classmethod
    def shot_cleanup(cls, shot_props : ShotProps,):
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
            
    def _retrieve_existing_data(
        self,
        shot_id : int,
        tokamak : Tokamak,
        shot_settings : ShotSettings,
    ) -> pd.DataFrame:
        if shot_settings.existing_data_request is not None:
            existing_data_request_params = ExistingDataRequestParams(
                shot_id=shot_id,
                database=self.database,
                tokamak=tokamak, 
                logger=self.logger,
            )
            existing_data = shot_settings.existing_data_request.get_existing_data(existing_data_request_params)
            existing_data['shot'] = existing_data['shot'].astype(int)
            existing_data = existing_data[existing_data['shot'] == shot_id]
        else:
            existing_data = None
        return existing_data
    
    def _init_times(
        self,
        shot_id : int,
        existing_data : pd.DataFrame, 
        tree_manager : TreeManager,
        tokamak : Tokamak,
        disruption_time : float,
        shot_settings : ShotSettings,
    ) -> np.ndarray:
        """
        Initialize the timebase of the shot.
        """
        request_params = SetTimesRequestParams(
            shot_id=shot_id, 
            tree_manager=tree_manager, 
            existing_data=existing_data,
            database=self.database, 
            disruption_time=disruption_time, 
            tokamak=tokamak, 
            logger=self.logger,
        )
        return shot_settings.set_times_request.get_times(request_params)
    
    @classmethod
    def _pre_fill_shot_data(cls, times : np.ndarray, existing_data : pd.DataFrame) -> pd.DataFrame:
        '''
        Intialize the shot with data, if existing data matches the shot timebase.
        '''
        if existing_data is not None:
            time_df = pd.DataFrame(times, columns=['time'])
            flagged_existing_data = existing_data.assign(merge_success_flag=1)
            timed_existing_data = pd.merge_asof(time_df, flagged_existing_data, on='time', direction='nearest', tolerance=TIME_CONST)
            if not timed_existing_data['merge_success_flag'].isna().any():
                return timed_existing_data.drop(columns=['merge_success_flag'])
            else:
                return None
        else:
            return None