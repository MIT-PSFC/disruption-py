from dataclasses import dataclass, field
import pandas as pd
from typing import List, Union, Callable, Tuple
from disruption_py.settings.enum_options import InterpolationMethod, SignalDomain
from disruption_py.settings.log_settings import LogSettings
from disruption_py.settings.existing_data_request import ExistingDataRequest, resolve_existing_data_request
from disruption_py.settings.shot_data_request import ShotDataRequest
from disruption_py.settings.set_times_request import SetTimesRequest, resolve_set_times_request
from disruption_py.settings.output_type_request import OutputTypeRequest, resolve_output_type_request
from disruption_py.utils.mappings.mappings_helpers import map_string_attributes_to_enum
from disruption_py.utils.utils import instantiate_classes, without_duplicates

def default_tags():
    return ["all"]
    
@dataclass
class ShotSettings:
    """Settings to be used for retrieving data for a single shot.

    Attributes
    ----------
    log_settings : LogSettings
        Settings for logging.
    existing_data_request : ExistingDataRequest
        The existing data request to be used when prefilling data for the shot. Can pass any 
        ExistingDataRequestType that resolves to a ExistingDataRequest. See ExistingDataRequest for more 
        details. Set to None if no data should be prefilled. Defaults to None.
    num_threads_per_shot : int
        Use not recommended. Alternatively, please see num_processes in get_shots_data.
    efit_tree_name : str
        The name of the tree to first try for the efit environment. Other tree names will be tried if 
        opening this tree name fails. Default is 'analysis'.
    attempt_local_efit_env : Tuple[Tuple[str, str]]
        A list of tuples of the form (efit_env_name, efit_env_path) that will be used to set the 
        environment variables when trying to open the efit tree. If opening the efit tree with the 
        local environment variables fails, will try to open the efit tree with the  regular environment 
        variables. Default is None.
    run_methods : List[str]
        A list of parameter method names to be run. Named methods will be run when retrieving data 
        from  mdsplus for the shot. Named methods must have with the parameter_cached_method 
        decorator and can either be located in the shot object or in a shot_data_request. Defaults 
        to an empty list.
    run_tags : List[str]
        A list of parameter method tags to be run. Methods used for retrieving data from mdsplus can be
        tagged with the parameter_cached_method decorator and can be located in either the shot 
        class or in an included shot_data_request. All methods with at least one included tag will 
        be run. Defaults to ["all"].
    run_columns : List[str]
        A list of columns to be retrieved. All methods with parameter_cached_method decorator referenced 
        as containing an included column will be run and all columns returned by those methods will be used. 
        Methods can either be located in the shot class or in an included shot_data_request. If you wish to
        only return the requested columns, set only_requested_columns to true in the shot_settings.
    only_requested_columns : bool
        Whether only columns requested in run_columns should be included in the produced dataframe.
        Even if not all requested columns exist in the produced dataframe only the requested columns will
        be produced. Otherwise all columns returned by all methods run will be included. Default false.
    shot_data_requests : List[ShotDataRequest]
        A list of ShotDataRequest objects. Methods from these objects are run when retrieving data 
        from mdsplus if the method is included through either the run_methods or run_tags setting.
        Defaults to an empty list.
    output_type_request : OutputTypeRequest
        The output type request to be used when outputting the retrieved data for each shot. Note that data
        is streamed to the output type request object as it is retrieved. Can pass any OutputTypeRequestType 
        that resolves to an OutputTypeRequest. See OutputTypeRequest for more details. Defaults to "list".
    set_times_request : SetTimesRequest
        The set times request to be used when setting the timebase for the shot. The retrieved data will
        be interpolated to this timebase. Can pass any SetTimesRequestType that resolves to a SetTimesRequest.
        See SetTimesRequest for more details. Defaults to "efit".
    signal_domain : SignalDomain
        The domain of the timebase that should be used when retrieving data for the shot. Either "full", 
        "flattop", or "rampup_and_flattop". Can pass either a SignalDomain or the associated string. Defaults 
        to "full".
    override_exising_data : bool
        If false and existing data exists for the shot, the timebase from the existing data will be used instead
        of the timebase from the set_times_request. If true the timebase from the set_times_request is always used.
        Defaults to True.
    interpolation_method : InterpolationMethod
        The interpolation method to be used when retrieving data for the shot. CURRENTLY UNIMPLEMENTED.
    
    Methods
    -------
    resolve()
        Resolve all attributes that can be resolved to a class. This primarily refers to passed strings
        lists and dictinoaries that can be resolved to a specific request type or a specific enum.
    """
    # General Settings
    log_settings : LogSettings = field(default_factory=LogSettings.default)
    
    # Prefill data settings
    existing_data_request : ExistingDataRequest = None
        
    # Shot creation settings
    num_threads_per_shot: int = 1
    efit_tree_name: str = 'analysis'
    attempt_local_efit_env: Tuple[Tuple[str, str]] = None
    
    # Shot run settings
    run_methods : List[str] = field(default_factory=list)
    run_tags : List[str] = field(default_factory=default_tags)
    run_columns : List[str] = field(default_factory=list)
    only_requested_columns : bool = False
    shot_data_requests : List[ShotDataRequest] = field(default_factory=list)
    
    # Shot Output settings
    output_type_request: OutputTypeRequest = "list"

    # Timebase setting
    set_times_request : SetTimesRequest = "efit"
    signal_domain : SignalDomain = "full"
    override_exising_data : bool = True
    interpolation_method : InterpolationMethod = "linear"
    
    additional_args : dict = field(default_factory=dict)
    
    
    def __post_init__(self):
        self.resolve()
        
    def resolve(self):
        self.existing_data_request = resolve_existing_data_request(self.existing_data_request)
        self.output_type_request = resolve_output_type_request(self.output_type_request)
        self.set_times_request = resolve_set_times_request(self.set_times_request)
        
        map_string_attributes_to_enum(self, {
            "signal_type": SignalDomain,
            "interpolation_method": InterpolationMethod
        })
        
        if self.attempt_local_efit_env is not None:
            for idx, (option, value) in enumerate(self.attempt_local_efit_env):
                if not option.endswith("_path"):
                    self.attempt_local_efit_env[idx] = (option + "_path", value)
        
        # we can also setup logging on resolve
        self.log_settings.setup_logging()