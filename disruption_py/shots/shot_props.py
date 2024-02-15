import logging
from typing import Dict, Any
import pandas as pd
import numpy as np
from dataclasses import dataclass, field
from disruption_py.mdsplus_integration.mds_connection import MDSConnection
from disruption_py.utils.mappings.tokamak import Tokamak


@dataclass
class ShotProps:
    logger = logging.getLogger('disruption_py')
    
    shot_id : int
    tokamak : Tokamak
    disruption_time : float
    mds_conn : MDSConnection
    times : np.ndarray
    existing_data : pd.DataFrame # existing data passed to shot class
    pre_filled_shot_data : pd.DataFrame # existing data after changed to times domain
    interpolation_method : Any # Fix
    metadata : dict
    
    _cached_results : Dict[str, Any] = field(default_factory=dict)
 
    @property
    def disrupted(self):
        return self.disruption_time is not None
    
    def cleanup(self):
        self.mds_conn.close_all_trees()
        self.times = None
        self._cached_results.clear()
        