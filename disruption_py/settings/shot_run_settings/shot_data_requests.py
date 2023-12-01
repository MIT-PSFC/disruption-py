import pandas as pd
from abc import ABC, abstractmethod
from dataclasses import dataclass
import os
from typing import List, Callable
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
    
    def get_request_methods_for_tokemak(self, tokemak: Tokemak) -> List[Callable]:
        request_methods = []
        for method in dir(self):
            if hasattr(method, "tokemaks") and (tokemak in method.tokemaks or tokemak is method.tokemaks):
                request_methods.append(method)
        return request_methods
