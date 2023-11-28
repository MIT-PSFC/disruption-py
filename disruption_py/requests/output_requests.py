import pandas as pd
from abc import ABC, abstractmethod
from dataclasses import dataclass
import os
from typing import Dict, Type
from logging import Logger
from disruption_py.utils.mappings.tokemak import Tokemak

@dataclass
class OutputTypeRequestParams:
    result : pd.DataFrame
    tokemak : Tokemak
    logger : Logger
    
class OutputTypeRequest(ABC):
    
    def output_shot(self, params : OutputTypeRequestParams):
        if hasattr(self, 'tokemak_overrides'):
            if params.tokemak in self.tokemak_overrides:
                return self.tokemak_overrides[params.tokemak](params)
        return self._output_shot(params)
    
    @abstractmethod
    def _output_shot(self, params : OutputTypeRequestParams):
        pass
    
    def stream_output_cleanup(self):
        pass
    
    @abstractmethod
    def get_results(self):
        pass


class ListOutputRequest(OutputTypeRequest):
    def __init__(self):
        self.results = []
        
    def _output_shot(self, params : OutputTypeRequestParams):
        self.results.append(params.result)
    
    def get_results(self):
        return self.results
    
    
class HDF5OutputRequest(OutputTypeRequest):
    def __init__(self, filepath):
        self.filepath = filepath
        self.store = pd.HDFStore(filepath, mode='w')
        self.output_shot_count = 0

    def _output_shot(self, params : OutputTypeRequestParams):
        shot_id = params.result['shot'].iloc[0] if (not params.result.empty and ('shot' in params.result.columns)) else self.output_shot_count
        self.store.append(f'df_{shot_id}', params.result, format='table', data_columns=True)
        self.output_shot_count += 1
    
    def stream_output_cleanup(self):
        self.store.close()
    
    def get_results(self):
        return self.output_shot_count
    
    
class CSVOutputRequest(OutputTypeRequest):
    def __init__(self, filepath, flexible_columns=False):
        self.filepath = filepath
        self.flexible_columns = flexible_columns
        self.output_shot_count = 0

    def _output_shot(self, params : OutputTypeRequestParams):
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

    def get_results(self):
        return self.output_shot_count

_output_type_request_mappings: Dict[str, OutputTypeRequest] = {
    # do not include classes that require initialization arguments
    "list" : ListOutputRequest(),
} 

_file_suffix_to_output_type_request : Dict[str, Type[OutputTypeRequest]] = {
    ".h5" : HDF5OutputRequest,
    ".hdf5" : HDF5OutputRequest,
    ".csv" : CSVOutputRequest,
} 

def output_type_request_runner(output_type_request, params : OutputTypeRequestParams):
    if isinstance(output_type_request, OutputTypeRequest):
        return output_type_request.output_shot(params)
    
    if isinstance(output_type_request, str):
        output_type_request_object = _output_type_request_mappings.get(output_type_request, None)
        if output_type_request_object is not None:
            return output_type_request_object.output_shot(params)
        
    if isinstance(output_type_request, str):
        # assume that it is a file path
       for suffix, output_type_request_type in _file_suffix_to_output_type_request.items():
           if output_type_request.endswith(suffix):
               output_type_request_object = output_type_request_type(output_type_request)
               return output_type_request_type().output_shot(params)       
            
    if isinstance(output_type_request, dict):
        chosen_request = output_type_request.get(params.tokemak, None)
        if chosen_request is not None:
            return output_type_request_runner(chosen_request, params)
        
    if isinstance(output_type_request, list):
        all_results = []
        for individual_type_request in output_type_request:
            sub_result = output_type_request_runner(individual_type_request, params)
            if sub_result is not None:
                all_results.append(sub_result)
        
        return all_results
    
    raise ValueError(f"Invalid output processror {output_type_request}")

