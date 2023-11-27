import pandas as pd
from abc import ABC, abstractmethod
from dataclasses import dataclass
import os
from typing import Dict
from logging import Logger
from disruption_py.utils.mappings.tokemak import Tokemak

@dataclass
class OutputProcessorParams:
    result : pd.DataFrame
    tokemak : Tokemak
    logger : Logger
    
class OutputProcessor(ABC):
    
    def output_shot(self, params : OutputProcessorParams):
        if hasattr(self, 'tokemak_overrides'):
            if params.tokemak in self.tokemak_overrides:
                return self.tokemak_overrides[params.tokemak](params)
        return self._get_shot_numbers(params)
    
    @abstractmethod
    def _output_shot(self, params : OutputProcessorParams):
        pass
    
    def stream_output_cleanup(self):
        pass
    
    @abstractmethod
    def get_results(self):
        pass


class ListOutputProcessor(OutputProcessor):
    def __init__(self):
        self.results = []
        
    def _output_shot(self, params : OutputProcessorParams):
        self.results.append(params.result)
    
    def get_results(self):
        return self.results
    
    
class HDF5OutputProcessor(OutputProcessor):
    def __init__(self, filepath):
        self.filepath = filepath
        self.store = pd.HDFStore(filepath, mode='w')
        self.output_shot_count = 0

    def _output_shot(self, params : OutputProcessorParams):
        shot_id = params.result['shot'].iloc[0] if (not params.result.empty and ('shot' in params.result.columns)) else self.output_shot_count
        self.store.append(f'df_{shot_id}', params.result, format='table', data_columns=True)
        self.output_shot_count += 1
    
    def stream_output_cleanup(self):
        self.store.close()
    
    def get_results(self):
        return self.output_shot_count
    
    
class CSVOutputProcessor(OutputProcessor):
    def __init__(self, filepath, flexible_columns=False):
        self.filepath = filepath
        self.flexible_columns = flexible_columns
        self.output_shot_count = 0

    def _output_shot(self, params : OutputProcessorParams):
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

_output_processor_mappings: Dict[str, OutputProcessor] = {
    # do not include classes that require initialization arguments
    "list" : ListOutputProcessor(),
} 

_file_suffix_mappings_to_output_processor = Dict[str, OutputProcessor] = {
    ".h5" : HDF5OutputProcessor,
    ".hdf5" : HDF5OutputProcessor,
    ".csv" : CSVOutputProcessor,
} 

def output_processor_runner(output_processor, params : OutputProcessorParams):
    if isinstance(output_processor, OutputProcessor):
        return output_processor.output_shot(params)
    
    if isinstance(output_processor, str):
        output_processor_object = _output_processor_mappings.get(output_processor, None)
        if output_processor_object is not None:
            return output_processor_object.output_shot(params)
        
    if isinstance(output_processor, str):
        # assume that it is a file path
       for suffix, output_processor_object in _file_suffix_mappings_to_output_processor.items():
           if output_processor.endswith(suffix):
               return _file_suffix_mappings_to_output_processor[output_processor].output_shot(params)       
            
    if isinstance(output_processor, dict):
        chosen_request = output_processor.get(params.tokemak, None)
        if chosen_request is not None:
            return output_processor_runner(chosen_request, params)
        
    if isinstance(output_processor, list):
        all_results = []
        for individual_processor in output_processor:
            sub_result = output_processor_runner(individual_processor, params)
            if sub_result is not None:
                all_results.append(sub_result)
        
        return all_results
    
    raise ValueError(f"Invalid output processror {output_processor}")

