import pandas as pd
from abc import ABC, abstractmethod
from dataclasses import dataclass
import os
from typing import Dict, List
from logging import Logger
from disruption_py.utils.mappings.tokemak import Tokemak
from disruption_py.shots import Shot

@dataclass
class ShotDataRequestParams:
    shot : Shot
    existing_data : pd.DataFrame
    tokemak : Tokemak
    logger : Logger
    
class ShotDataRequest(ABC):
    
    def get_shot_data(self, params : ShotDataRequestParams):
        if hasattr(self, 'tokemak_overrides'):
            if params.tokemak in self.tokemak_overrides:
                return self.tokemak_overrides[params.tokemak](params)
        return self._get_shot_data(params)
    
    @abstractmethod
    def _get_shot_data(self, params : ShotDataRequestParams) -> List[int]:
        pass