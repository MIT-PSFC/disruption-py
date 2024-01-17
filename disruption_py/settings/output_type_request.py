import pandas as pd
from abc import ABC, abstractmethod
from dataclasses import dataclass
import os
from typing import Any, Dict, List, Type, Union
from logging import Logger
from disruption_py.databases.database import ShotDatabase
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.mappings.mappings_helpers import map_string_to_enum


@dataclass
class ResultOutputTypeRequestParams:
    """Params passed by disruption_py to _output_shot() method.

    Attributes
    ----------
    result : pd.Dataframe
        The dataframe of results for a single shot.
    database : ShotDatabase
        Database object to use for getting existing data.
        A different database connection is used by each thread/process.
    tokamak : Tokemak
        The tokamak for which the data request is made.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """
    result : pd.DataFrame
    database : ShotDatabase
    tokamak : Tokamak
    logger : Logger
    
@dataclass
class FinishOutputTypeRequestParams:
    """Params passed by disruption_py to stream_output_cleanup() and get_results methods.

    Attributes
    ----------
    tokamak : Tokemak
        The tokamak for which the data request is made.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """
    tokamak : Tokamak
    logger : Logger

OutputTypeRequestType = Union['OutputTypeRequest', str, Dict[str, 'OutputTypeRequestType'], List['OutputTypeRequestType']]

class OutputTypeRequest(ABC):
    """OutputTypeRequest abstract class that should be inherited by all output type request classes.
    """
    
    def output_shot(self, params : ResultOutputTypeRequestParams):
        if hasattr(self, 'tokamak_overrides'):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._output_shot(params)
    
    @abstractmethod
    def _output_shot(self, params : ResultOutputTypeRequestParams):
        """Abstract method implemented by subclasses to handle data output for a single shot.
        This method is called by disruption_py with the shots dataframe in the params object 
        once the data has been retrieved.
        
        Parameters
        ----------
        params : ResultOutputTypeRequestParams
            Params containing the data retrieved for a shot in a dataframe and other utility parameters.
        """
        pass
    
    def stream_output_cleanup(self, params: FinishOutputTypeRequestParams):
        """Empty method optionally overriden by subclasses to handle cleanup after all shots have been 
        output. This may include closing files or other cleanup.
        
        Parameters
        ----------
        params : FinishOutputTypeRequestParams
            Utility parameters such as the tokamak and logger.
        """
        pass
    
    @abstractmethod
    def get_results(self, params: FinishOutputTypeRequestParams) -> Any:
        """Abstract method implemented by subclasses to handle the output of the data from
        calls to `get_shots_data`. This method is called by disruption_py once `output_shot()` has been
        called for all shots ids in the shot ids request.
        
        Parameters
        ----------
        params : FinishOutputTypeRequestParams
            Utility parameters such as the tokamak and logger.
        
        Returns
        -------
        Any
            The desired output of the call to `get_shots_data` potentially containing the data for all shots,
            some aggregation of that data, or nothing.
        """
        pass

class OutputTypeRequestList(OutputTypeRequest):
    """
    Utility class that is automatically used when a list is passed as the `output_type_request` parameter in `ShotSettings.
    
    All listed output types will be output to in the order listed.
    Similarly, results will be returned in the order listed.
    
    Parameters
    ----------
    output_type_request_list : list[OutputTypeRequestType]
        A python list of any other output type request option that can be passed as the `output_type_request` parameter in `ShotSettings`.
        Any other option passable to the `output_type_request` parameter in `ShotSettings` may be used.
    """
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
    """
    Utility class that is automatically used when a dicationary is passed as the `output_type_request` parameter in `ShotSettings.
    
    Parameters
    ----------
    output_type_request_dict : dict[Tokamak, OutputTypeRequestType]
        A dictionary mapping tokamak type strings to the desired `OutputTypeRequestType` for that tokamak.  E.g. `{'cmod': 'list'}`.
    """
    def __init__(self, output_type_request_dict : Dict[Tokamak, OutputTypeRequestType]):
        resolved_output_type_request_dict = {
            map_string_to_enum(tokamak, Tokamak): resolve_output_type_request(individual_type_request) 
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
    """
    Output all retrieved shot data as a list of dataframes, once retrieval complete.
    """
    def __init__(self):
        self.results = []
        
    def _output_shot(self, params : ResultOutputTypeRequestParams):
        self.results.append(params.result)
    
    def get_results(self, params: FinishOutputTypeRequestParams):
        return self.results
    
    
class HDF5OutputRequest(OutputTypeRequest):
    """
    Stream outputted data to an HDF5 file.
    """
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
    """
    Stream outputted data to a single csv file.
    """
    def __init__(self, filepath, flexible_columns=True, clear_file=True):
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

# --8<-- [start:output_type_request_dict]
_output_type_request_mappings: Dict[str, OutputTypeRequest] = {
    "list" : ListOutputRequest(),
}
# --8<-- [end:output_type_request_dict]

# --8<-- [start:file_suffix_to_output_type_request_dict]
_file_suffix_to_output_type_request : Dict[str, Type[OutputTypeRequest]] = {
    ".h5" : HDF5OutputRequest,
    ".hdf5" : HDF5OutputRequest,
    ".csv" : CSVOutputRequest,
} 
# --8<-- [end:file_suffix_to_output_type_request_dict]

def resolve_output_type_request(output_type_request : OutputTypeRequestType) -> OutputTypeRequest:
    if isinstance(output_type_request, OutputTypeRequest):
        return output_type_request
    
    if isinstance(output_type_request, str):
        output_type_request_object = _output_type_request_mappings.get(output_type_request, None)
        if output_type_request_object is not None:
            return output_type_request_object
        
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
    

