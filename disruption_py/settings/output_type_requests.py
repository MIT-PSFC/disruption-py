import pandas as pd
from abc import ABC, abstractmethod
from dataclasses import dataclass
import os
from typing import Dict, List, Type, Union
from logging import Logger
from disruption_py.utils.mappings.tokamak import Tokemak


@dataclass
class ResultOutputTypeRequestParams:
    result : pd.DataFrame
    tokamak : Tokemak
    logger : Logger
    
@dataclass
class FinishOutputTypeRequestParams:
    tokamak : Tokemak
    logger : Logger

OutputTypeRequestType = Union['OutputTypeRequest', str, Dict[str, 'OutputTypeRequestType'], List['OutputTypeRequestType']]

class OutputTypeRequest(ABC):
    
    def output_shot(self, params : ResultOutputTypeRequestParams):
        if hasattr(self, 'tokamak_overrides'):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._output_shot(params)
    
    @abstractmethod
    def _output_shot(self, params : ResultOutputTypeRequestParams):
        pass
    
    def stream_output_cleanup(self, params: FinishOutputTypeRequestParams):
        pass
    
    @abstractmethod
    def get_results(self, params: FinishOutputTypeRequestParams):
        pass

class OutputTypeRequestList(OutputTypeRequest):
    '''
    List of output type requests
    '''
    def __init__(self, output_type_request_list : List[OutputTypeRequestType]):
        self.output_type_request_list = [resolve_output_type_request(individual_type_request) for individual_type_request in output_type_request_list]

    def _output_shot(self, params : ResultOutputTypeRequestParams):
        all_results = []
        for individual_type_request in self.output_type_request_list:
            sub_result = individual_type_request.output_shot(params)
            all_results.append(sub_result)
                
        return all_results
    
    def stream_output_cleanup(self, params: FinishOutputTypeRequestParams):
         for individual_type_request in self.output_type_request_list:
            individual_type_request.stream_output_cleanup()
            
    def get_results(self, params: FinishOutputTypeRequestParams):
        return [individual_type_request.get_results() for individual_type_request in self.output_type_request_list]
    

class OutputTypeRequestDict(OutputTypeRequest):
    '''
    Dict of output type requests
    '''
    def __init__(self, output_type_request_dict : Dict[Tokemak, OutputTypeRequestType]):
        resolved_output_type_request_dict = {
            tokamak: resolve_output_type_request(individual_type_request) 
            for tokamak, individual_type_request in output_type_request_dict.items()
        }
        self.output_type_request_dict = resolved_output_type_request_dict

    def _output_shot(self, params : ResultOutputTypeRequestParams):
        chosen_request = self.output_type_request_dict.get(params.tokamak, None)
        if chosen_request is not None:
            return chosen_request.output_shot(params)
        else:
            params.logger.warning(f'No output type request for tokamak {params.tokamak}')
            return None
    
    def stream_output_cleanup(self, params: FinishOutputTypeRequestParams):
        chosen_request = self.output_type_request_dict.get(params.tokamak, None)
        if chosen_request is not None:
            return chosen_request.stream_output_cleanup(params)
        else:
            params.logger.warning(f'No output type request for tokamak {params.tokamak}')
            return None
            
    def get_results(self, params: FinishOutputTypeRequestParams):
        chosen_request = self.output_type_request_dict.get(params.tokamak, None)
        if chosen_request is not None:
            return chosen_request.stream_output_cleanup(params.tokamak, params.logger)
        else:
            params.logger.warning(f'No output type request for tokamak {params.tokamak}')
            return None    

class ListOutputRequest(OutputTypeRequest):
    def __init__(self):
        self.results = []
        
    def _output_shot(self, params : ResultOutputTypeRequestParams):
        self.results.append(params.result)
    
    def get_results(self, params: FinishOutputTypeRequestParams):
        return self.results
    
    
class HDF5OutputRequest(OutputTypeRequest):
    def __init__(self, filepath):
        self.filepath = filepath
        self.store = pd.HDFStore(filepath, mode='w')
        self.output_shot_count = 0

    def _output_shot(self, params : ResultOutputTypeRequestParams):
        shot_id = params.result['shot'].iloc[0] if (not params.result.empty and ('shot' in params.result.columns)) else self.output_shot_count
        self.store.append(f'df_{shot_id}', params.result, format='table', data_columns=True)
        self.output_shot_count += 1
    
    def stream_output_cleanup(self, params: FinishOutputTypeRequestParams):
        self.store.close()
    
    def get_results(self, params: FinishOutputTypeRequestParams):
        return self.output_shot_count
    
    
class CSVOutputRequest(OutputTypeRequest):
    def __init__(self, filepath, flexible_columns=False, clear_file=True):
        self.filepath = filepath
        self.flexible_columns = flexible_columns
        self.output_shot_count = 0
        if clear_file is True and os.path.exists(filepath):
            os.remove(filepath)

    def _output_shot(self, params : ResultOutputTypeRequestParams):
        file_exists = os.path.isfile(self.filepath)
        if self.flexible_columns:
            if file_exists:
                existing_df = pd.read_csv(self.filepath)
                combined_df = pd.concat([existing_df, params.result], ignore_index=True, sort=False)
            else:
                combined_df = params.result

            combined_df.to_csv(self.filepath, index=False)
        else: 
            params.result.to_csv(self.filepath, mode='a', index=False, header=(not file_exists))
        self.output_shot_count += 1

    def get_results(self, params: FinishOutputTypeRequestParams):
        return self.output_shot_count

_output_type_request_mappings: Dict[str, Type[OutputTypeRequest]] = {
    # do not include classes that require initialization arguments
    "list" : ListOutputRequest,
} 

_file_suffix_to_output_type_request : Dict[str, Type[OutputTypeRequest]] = {
    ".h5" : HDF5OutputRequest,
    ".hdf5" : HDF5OutputRequest,
    ".csv" : CSVOutputRequest,
} 

def resolve_output_type_request(output_type_request : OutputTypeRequestType) -> OutputTypeRequest:
    if isinstance(output_type_request, OutputTypeRequest):
        return output_type_request
    
    if isinstance(output_type_request, str):
        output_type_request_type = _output_type_request_mappings.get(output_type_request, None)
        if output_type_request_type is not None:
            return output_type_request_type()
        
    if isinstance(output_type_request, str):
        # assume that it is a file path
       for suffix, output_type_request_type in _file_suffix_to_output_type_request.items():
           if output_type_request.endswith(suffix):
               return output_type_request_type(output_type_request)
            
    if isinstance(output_type_request, dict):
        return OutputTypeRequestDict(output_type_request)
        
    if isinstance(output_type_request, list):        
        return OutputTypeRequestList(output_type_request)
    
    raise ValueError(f"Invalid output processror {output_type_request}")
    

