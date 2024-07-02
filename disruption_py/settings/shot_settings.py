#!/usr/bin/env python3

from dataclasses import dataclass, field
from typing import List, Tuple

from disruption_py.settings.enum_options import InterpolationMethod, SignalDomain
from disruption_py.settings.existing_data_request import (
    ExistingDataRequest,
    resolve_existing_data_request,
)
from disruption_py.settings.output_type_request import OutputTypeRequest
from disruption_py.settings.set_times_request import (
    ExistingDataSetTimesRequest,
    SetTimesRequest,
    resolve_set_times_request,
)
from disruption_py.utils.mappings.mappings_helpers import map_string_attributes_to_enum
from disruption_py.machine.tokamak import Tokamak, is_tokamak_indexed


def default_tags():
    return ["all"]


@dataclass
class ShotSettings:
    """Settings to be used for retrieving data for a single shot.

    Attributes
    ----------
    existing_data_request : ExistingDataRequest
        The existing data request to be used when prefilling data for the shot. Can pass any
        ExistingDataRequestType that resolves to a ExistingDataRequest. See ExistingDataRequest for more
        details. Set to None if no data should be prefilled. Defaults to None.
    efit_tree_name : str
        The name of the tree to first try for the efit environment. Other tree names will be tried if
        opening this tree name fails. Default is 'analysis'.
    run_methods : List[str]
        A list of parameter method names to be run. Named methods will be run when retrieving data
        from  mdsplus for the shot. Named methods must have the parameter_method decorator and either
        be passed in the `custom_parameter_methods` argument or included in the built-in list. Defaults
        to an empty list.
    run_tags : List[str]
        A list of parameter method tags to be run. Methods used for retrieving data from mdsplus can be
        tagged with the parameter_method decorator and can either be passed in the `custom_parameter_methods`
        argument or included in the built-in list. All methods with at least one included tag will be run.
        Defaults to ["all"].
    run_columns : List[str]
        A list of columns to be retrieved. All methods with the parameter_method decorator referenced
        as containing an included column will be run and all columns returned by those methods will be used.
        Methods can either be passed in the `custom_parameter_methods` argument or included in the built-in
        list. If you wish to only return the requested columns, set only_requested_columns to true in the
        shot_settings.
    only_requested_columns : bool
        Whether only columns requested in run_columns should be included in the produced dataframe.
        Even if not all requested columns exist in the produced dataframe only the requested columns will
        be produced. Otherwise all columns returned by all methods run will be included. Default false.
    custom_parameter_methods : list
        A list of parametered methods and objects containing registred methods. The Methods are
        collected and run when retrieving data from mdsplus if the method is included through
        either the run_methods, run_tags, run_columns setting. Defaults to an empty list.
    set_times_request : SetTimesRequest
        The set times request to be used when setting the timebase for the shot. The retrieved data will
        be interpolated to this timebase. Can pass any SetTimesRequestType that resolves to a SetTimesRequest.
        See SetTimesRequest for more details. Defaults to "disruption_warning".
    signal_domain : SignalDomain
        The domain of the timebase that should be used when retrieving data for the shot. Either "full",
        "flattop", or "rampup_and_flattop". Can pass either a SignalDomain or the associated string. Defaults
        to "full".
    use_existing_data_timebase : bool
        If true and existing data exists for the shot, the timebase from the existing data will be used instead
        of the timebase from the set_times_request. Wraps the set_times_request with ExistingDataSetTimesRequest.
        Defaults to False.
    interpolation_method : InterpolationMethod
        The interpolation method to be used when retrieving data for the shot. CURRENTLY UNIMPLEMENTED.
    output_type_request : OutputTypeRequest
        DEPRECTATED. output_type_request has moved to a parameter in the get_shots_data method.
        Will error if used, please set to None.
    attempt_local_efit_env : Tuple[Tuple[str, str]]
        DEPRECTATED. Support no longer exists. Please reach out to maintainers with questions.
        A list of tuples of the form (efit_env_name, efit_env_path) that will be used to set the
        environment variables when trying to open the efit tree. If opening the efit tree with the
        local environment variables fails, will try to open the efit tree with the  regular environment
        variables. Default is None.
    """

    # Prefill data settings
    existing_data_request: ExistingDataRequest = None

    # Shot creation settings
    efit_tree_name: str = "analysis"

    # Shot run settings
    run_methods: List[str] = field(default_factory=list)
    run_tags: List[str] = field(default_factory=default_tags)
    run_columns: List[str] = field(default_factory=list)
    only_requested_columns: bool = False
    custom_parameter_methods: list = field(default_factory=list)

    # Timebase setting
    set_times_request: SetTimesRequest = "disruption_warning"
    signal_domain: SignalDomain = "full"
    use_existing_data_timebase: bool = False
    interpolation_method: InterpolationMethod = "linear"

    additional_args: dict = field(default_factory=dict)

    # DEPRECATED
    output_type_request: OutputTypeRequest = None  # moved to get_shots_data
    attempt_local_efit_env: Tuple[Tuple[str, str]] = None  # support removed

    def __post_init__(self):
        self.resolve()

    def _check_deprecated_settings(self):
        """
        Check that no deprectated settings are set in the shot settings.
        """
        DEPRECTATED_SETTINGS = [
            (
                "output_type_request",
                """
                output_type_request no longer set in shot_settings. 
                Please set output_type_request in get_shots_data. 
                To not throw error please set output_type_request to None.
                """
                "attempt_local_efit_env",
                """
                attempt_local_efit_env support no longer exists. Please reach out to maintainers with questions.
                """,
            ),
        ]

        for setting in DEPRECTATED_SETTINGS:
            if getattr(self, setting[0]) is not None:
                raise ValueError(f"{setting[1]}")

    @classmethod
    def from_dict(cls, prop_dict, tokamak: Tokamak):
        """
        Create a ShotSettings object from a dictionary.
        """
        if is_tokamak_indexed(prop_dict):
            if tokamak.value not in prop_dict:
                raise ValueError(
                    f"Tokamak {tokamak.value} not found in shot settings. Available tokamaks are {prop_dict.keys()}"
                )

            prop_dict = prop_dict[tokamak.value]

        return cls(**prop_dict)

    def resolve(self):
        """
        Take parameters that are passed preset values, and convert to value usable by disruption_py

        This primarily refers to passed strings lists and dictinoaries that can be resolved to a specific request type or a specific enum.
        """
        self._check_deprecated_settings()

        self.existing_data_request = resolve_existing_data_request(
            self.existing_data_request
        )
        self.set_times_request = resolve_set_times_request(self.set_times_request)

        map_string_attributes_to_enum(
            self,
            {
                "signal_domain": SignalDomain,
                "interpolation_method": InterpolationMethod,
            },
        )

        if self.use_existing_data_timebase and not isinstance(
            self.set_times_request, ExistingDataSetTimesRequest
        ):
            self.set_times_request = ExistingDataSetTimesRequest(self.set_times_request)

        if self.attempt_local_efit_env is not None:
            for idx, (option, value) in enumerate(self.attempt_local_efit_env):
                if not option.endswith("_path"):
                    self.attempt_local_efit_env[idx] = (option + "_path", value)
