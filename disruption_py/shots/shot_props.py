from typing import Dict, Any
import pandas as pd
from dataclasses import dataclass
from disruption_py.mdsplus_integration.tree_manager import TreeManager
from disruption_py.utils.mappings.tokamak import Tokamak


@dataclass
class ShotProps:
    shot_id : str
    tokamak : Tokamak
    num_threads_per_shot : int
    disruption_time : float
    tree_manager : TreeManager
    initial_existing_data : pd.DataFrame # existing data passed to shot class
    processed_existing_data : pd.DataFrame # existing data after changed to times domain
    interpolation_method : Any # Fix
    metadata : dict
    
    _cached_results : Dict[str, Any] = {}
 
    @property
    def disrupted(self):
        return self.disruption_time is not None
    
    def cleanup(self):
        self.tree_manager.cleanup()
        self.times = None
        self._cached_results.clear()