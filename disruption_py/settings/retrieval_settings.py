#!/usr/bin/env python3

from dataclasses import dataclass, field
from enum import Enum
from typing import List

from disruption_py.core.utils.enums import map_string_attributes_to_enum
from disruption_py.machine.tokamak import Tokamak, is_tokamak_indexed
from disruption_py.settings.cache_setting import CacheSetting, resolve_cache_setting
from disruption_py.settings.domain_setting import DomainSetting, resolve_domain_setting
from disruption_py.settings.nickname_setting import (
    NicknameSetting,
    resolve_nickname_setting,
)
from disruption_py.settings.time_setting import (
    CacheTimeSetting,
    TimeSetting,
    resolve_time_setting,
)


def default_tags():
    return ["all"]


class InterpolationMethod(Enum):
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
    """RetrievalSettings to be used for retrieving data for a single shot.

    Attributes
    ----------
    cache_setting : CacheSetting
        The cache setting to be used when prefilling data for the shot. Can pass any
        CacheSettingType that resolves to a CacheSetting. See CacheSetting for more
        details. Set to None if no data should be prefilled. Defaults to None.
    efit_nickname_setting : NicknameSetting
        The nickname setting to be used when retrieving data for the shot. Defines the nicknames that
        should be created for the efit_tree.
    run_methods : List[str]
        A list of physics method names to be run. Named methods will be run when retrieving data
        from  mdsplus for the shot. Named methods must have the physics_method decorator and either
        be passed in the `custom_physics_methods` argument or included in the built-in list. Defaults
        to an empty list.
    run_tags : List[str]
        A list of physics method tags to be run. Methods used for retrieving data from mdsplus can be
        tagged with the physics_method decorator and can either be passed in the `custom_physics_methods`
        argument or included in the built-in list. All methods with at least one included tag will be run.
        Defaults to ["all"].
    run_columns : List[str]
        A list of columns to be retrieved. All methods with the physics_method decorator referenced
        as containing an included column will be run and all columns returned by those methods will be used.
        Methods can either be passed in the `custom_physics_methods` argument or included in the built-in
        list. If you wish to only return the requested columns, set only_requested_columns to true in the
        retrieval_settings.
    only_requested_columns : bool
        Whether only columns requested in run_columns should be included in the produced dataframe.
        Even if not all requested columns exist in the produced dataframe only the requested columns will
        be produced. Otherwise all columns returned by all methods run will be included. Default false.
    custom_physics_methods : list
        A list of physics methods and objects containing physics methods. The Methods are
        collected and run when retrieving data from mdsplus if the method is included through
        either the run_methods, run_tags, run_columns setting. Defaults to an empty list.
    time_setting : TimeSetting
        The time setting to be used when setting the timebase for the shot. The retrieved data will
        be interpolated to this timebase. Can pass any `TimeSettingType` that resolves to a TimeSetting.
        See TimeSetting for more details. Defaults to "disruption_warning".
    domain_setting : DomainSetting
        The domain of the timebase that should be used when retrieving data for the shot. Either "full",
        "flattop", or "rampup_and_flattop". Can pass any `DomainSettingType` that resolves to a `DomainSetting`
        such as the listed strings. Defaults to "full".
    use_cache_setting_timebase : bool
        If true and cache data exists for the shot, the timebase from the cached data will be used instead
        of the timebase from the time_setting. Wraps the time_setting with CacheTimeSetting.
        Defaults to False.
    interpolation_method : InterpolationMethod
        The interpolation method to be used when retrieving data for the shot. CURRENTLY UNIMPLEMENTED.
    """

    # Prefill data settings
    cache_setting: CacheSetting = None

    # Shot creation settings
    efit_nickname_setting: NicknameSetting = "disruption"

    # Shot run settings
    run_methods: List[str] = field(default_factory=list)
    run_tags: List[str] = field(default_factory=default_tags)
    run_columns: List[str] = field(default_factory=list)
    only_requested_columns: bool = False
    custom_physics_methods: list = field(default_factory=list)

    # Timebase setting
    time_setting: TimeSetting = "disruption_warning"
    domain_setting: DomainSetting = "full"
    use_cache_setting_timebase: bool = False
    interpolation_method: InterpolationMethod = "linear"

    def __post_init__(self):
        self.resolve()

    @classmethod
    def from_dict(cls, prop_dict, tokamak: Tokamak):
        """
        Create a RetrievalSettings object from a dictionary.
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

        self.cache_setting = resolve_cache_setting(self.cache_setting)
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

        if self.use_cache_setting_timebase and not isinstance(
            self.time_setting, CacheTimeSetting
        ):
            self.time_setting = CacheTimeSetting(self.time_setting)
