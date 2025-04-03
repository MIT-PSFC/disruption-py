#!/usr/bin/env python3

"""
This module defines the RetrievalSettings class, which is used to configure
settings for retrieving data for a single shot.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import List

from disruption_py.core.utils.enums import map_string_attributes_to_enum
from disruption_py.machine.tokamak import Tokamak, is_tokamak_indexed
from disruption_py.settings.domain_setting import DomainSetting, resolve_domain_setting
from disruption_py.settings.nickname_setting import (
    NicknameSetting,
    resolve_nickname_setting,
)
from disruption_py.settings.time_setting import TimeSetting, resolve_time_setting


class InterpolationMethod(Enum):
    """Enum for specifying interpolation methods."""

    LINEAR = "linear"
    LOG = "log"
    LOG_LINEAR = "log_linear"
    EXP = "exp"
    EXP_LINEAR = "exp_linear"
    SIN = "sin"
    SIN_LINEAR = "sin_linear"
    SIN_EXP = "sin_exp"


@dataclass
class RetrievalSettings:
    """
    Settings for retrieving data for a single shot.

    Attributes
    ----------
    efit_nickname_setting : NicknameSetting, optional
        Nickname setting for retrieving efit tree data (default is "disruption").
    run_methods : list of str, optional
        List of physics methods to run (default is None). If None, and run_columns
        is None, all methods will be run. Named methods will be run when retrieving
        data from  MDSplus for the shot. Named methods must have the physics_method
        decorator and either be passed in the `custom_physics_methods` argument
        or included in the built-in method holders.
    run_columns : list of str, optional
        List of columns to retrieve (default is None). If None, and run_methods is
        None, all methods will be run. If specified, all methods with the physics_method
        decorator referencing the specified column will be run and all columns returned
        by those methods will be used. If you wish to only return the requested columns,
        set only_requested_columns to True in the retrieval_settings.
    only_requested_columns : bool, optional
        Whether to only include requested columns in the result (default is False).
    custom_physics_methods : list, optional
        List of custom physics methods (default is an empty list). The Methods are
        collected and run when retrieving data from MDSplus if the method is included
        through either the run_methods or run_columns setting.
    time_setting : TimeSetting, optional
        Time setting for the shot (default is "disruption_warning"). The retrieved
        data will be interpolated to this timebase. Can pass any `TimeSettingType`
        that resolves to a TimeSetting. See TimeSetting for more details.
    domain_setting : DomainSetting, optional
        Domain setting for the timebase (default is "full"). Either "full", "flattop",
        or "rampup_and_flattop". Can pass any `DomainSettingType` that resolves
        to a `DomainSetting` such as the listed strings.
    interpolation_method : InterpolationMethod, optional
        Interpolation method to be used (default is "linear"). CURRENTLY UNIMPLEMENTED.
    """

    # Shot creation settings
    efit_nickname_setting: NicknameSetting = "disruption"

    # Shot run settings
    run_methods: List[str] | None = None
    run_columns: List[str] | None = None
    only_requested_columns: bool = False
    custom_physics_methods: list = field(default_factory=list)

    # Timebase setting
    time_setting: TimeSetting = "disruption_warning"
    domain_setting: DomainSetting = "full"
    interpolation_method: InterpolationMethod = "linear"

    def __post_init__(self):
        """Resolve settings after initialization."""
        self.resolve()

    @classmethod
    def from_dict(cls, prop_dict, tokamak: Tokamak):
        """
        Create a RetrievalSettings object from a dictionary.

        Parameters
        ----------
        prop_dict : dict
            Dictionary containing properties for RetrievalSettings.
        tokamak : Tokamak
            Tokamak for which to create settings.

        Returns
        -------
        RetrievalSettings
            A configured RetrievalSettings object.
        """
        if is_tokamak_indexed(prop_dict):
            if tokamak.value not in prop_dict:
                raise ValueError(
                    f"Tokamak {tokamak.value} not found in shot settings. Available"
                    f" tokamaks are {prop_dict.keys()}"
                )

            prop_dict = prop_dict[tokamak.value]

        return cls(**prop_dict)

    def resolve(self):
        """
        Resolve preset values into specific objects or enums.

        This method resolves passed strings, lists, and dictionaries into specific
        request types or enums.
        """

        self.time_setting = resolve_time_setting(self.time_setting)
        self.domain_setting = resolve_domain_setting(self.domain_setting)
        self.efit_nickname_setting = resolve_nickname_setting(
            self.efit_nickname_setting
        )

        map_string_attributes_to_enum(
            self,
            {
                "interpolation_method": InterpolationMethod,
            },
        )

        if self.run_columns is not None:
            self.run_columns = [col.lower() for col in self.run_columns]
