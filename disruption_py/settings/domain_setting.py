#!/usr/bin/env python3

"""
This module provides classes and methods to specify different domain settings
for the timebase used, including full, flattop, and ramp-up domains. 
"""

import traceback
from abc import ABC, abstractmethod
from dataclasses import dataclass
from logging import Logger
from typing import Dict, Union

import numpy as np
from MDSplus import mdsExceptions

from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.enums import map_string_to_enum
from disruption_py.core.utils.math import interp1
from disruption_py.machine.cmod.physics import CmodPhysicsMethods
from disruption_py.machine.tokamak import Tokamak

DomainSettingType = Union["DomainSetting", str, Dict[Tokamak, "DomainSettingType"]]


@dataclass
class DomainSettingParams:
    """
    Parameters passed to the `get_domain` method.

    Attributes
    ----------
    physics_method_params : PhysicsMethodParams
        The physics method parameters that include timebase information.
    tokamak : Tokamak
        The tokamak instance for which the domain is set.
    logger : Logger
        Logger instance for logging relevant messages.
    """

    physics_method_params: PhysicsMethodParams
    tokamak: Tokamak
    logger: Logger


class DomainSetting(ABC):
    """
    Abstract base class for domain settings that modify the timebase
    domain in physics method parameters.
    """

    def get_domain(self, params: DomainSettingParams) -> np.ndarray:
        """
        Get the domain of the timebase to use.

        Parameters
        ----------
        params : DomainSettingParams
            The parameters including the physics method and tokamak information.

        Returns
        -------
        np.ndarray
            Array representing the timebase domain.
        """
        if hasattr(self, "tokamak_overrides"):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._get_domain(params)

    @abstractmethod
    def _get_domain(self, params: DomainSettingParams) -> np.ndarray:
        """
        Abstract method to get the modified times attribute of the physics method
        parameters based on the domain.

        Parameters
        ----------
        params : DomainSettingParams
            The parameters for domain modification.

        Returns
        -------
        np.ndarray
            Array of modified times.
        """


class DomainSettingDict(DomainSetting):
    """
    Domain setting class that handles a dictionary of tokamak-domain mappings.

    Parameters
    ----------
    domain_setting_dict : dict[Tokamak, DomainSettingType]
        A dictionary mapping tokamak type strings to the desired domain setting for
        that tokamak.  E.g. `{'cmod': 'flattop'}`. Any other option passable to the
        `domain_setting` parameter in `RetrievalSettings` may be used.
    """

    def __init__(self, domain_setting_dict: Dict[Tokamak, DomainSettingType]):
        """
        Initialize DomainSettingDict with a dictionary of tokamak-domain settings.

        Parameters
        ----------
        domain_setting_dict : dict
            Dictionary of tokamak to domain setting mappings.
        """
        resolved_domain_setting_dict = {
            map_string_to_enum(tokamak, Tokamak): resolve_domain_setting(
                individual_setting
            )
            for tokamak, individual_setting in domain_setting_dict.items()
        }
        self.resolved_domain_setting_dict = resolved_domain_setting_dict

    def _get_domain(self, params: DomainSettingParams) -> np.ndarray:
        """
        Get the domain for the given tokamak based on the resolved settings.

        Parameters
        ----------
        params : DomainSettingParams
            Parameters containing physics method and tokamak information.

        Returns
        -------
        np.ndarray
            The resolved timebase domain for the given tokamak.
        """
        chosen_setting = self.resolved_domain_setting_dict.get(params.tokamak, None)
        if chosen_setting is not None:
            return chosen_setting.get_domain(params)
        params.logger.warning("No domain setting for tokamak %s", params.tokamak)
        return None


class FullDomainSetting(DomainSetting):
    """
    Domain setting that uses the full timebase.
    """

    def _get_domain(self, params: DomainSettingParams) -> np.ndarray:
        """
        Get the full timebase domain.

        Parameters
        ----------
        params : DomainSettingParams
            Parameters containing physics method and tokamak information.

        Returns
        -------
        np.ndarray
            The full timebase.
        """
        return params.physics_method_params.times


class FlattopDomainSetting(DomainSetting):
    """
    Specifies that the flattop domain of the timebase should be used.
    """

    def __init__(self):
        """
        Initialize the FlattopDomainSetting with tokamak-specific overrides.
        """
        self.tokamak_overrides = {
            Tokamak.CMOD: self._get_domain_cmod,
            Tokamak.D3D: self._get_domain_d3d,
        }

    def _get_domain(self, params: DomainSettingParams) -> np.ndarray:
        """
        Get the flattop domain for the given tokamak.

        Parameters
        ----------
        params : DomainSettingParams
            Parameters containing physics method and tokamak information.

        Returns
        -------
        np.ndarray
            The flattop timebase domain.
        """
        raise ValueError(f"flattop domain not defined for tokamak: {params.tokamak}")

    def _get_domain_cmod(self, params: DomainSettingParams) -> np.ndarray:
        """
        Get the flattop domain for CMOD tokamak.

        Parameters
        ----------
        params : DomainSettingParams
            Parameters containing physics method and tokamak information.

        Returns
        -------
        np.ndarray
            The flattop timebase domain for CMOD.
        """
        ip_parameters = CmodPhysicsMethods.get_ip_parameters(
            params=params.physics_method_params
        )
        ipprog, dipprog_dt = ip_parameters["ip_prog"], ip_parameters["dipprog_dt"]
        indices_flattop_1 = np.where(np.abs(dipprog_dt) <= 1e3)[0]
        indices_flattop_2 = np.where(np.abs(ipprog) > 1.0e5)[0]
        indices_flattop = np.intersect1d(indices_flattop_1, indices_flattop_2)
        if len(indices_flattop) == 0:
            params.logger.warning(
                (
                    "[Shot %s]: Could not find flattop timebase. "
                    "Defaulting to full shot(efit) timebase."
                ),
                params.physics_method_params.shot_id,
            )
            return None
        return params.physics_method_params.times[indices_flattop]

    def _get_domain_d3d(self, params: DomainSettingParams) -> np.ndarray:
        """
        Get the flattop domain for D3D tokamak.

        Parameters
        ----------
        params : DomainSettingParams
            Parameters containing physics method and tokamak information.

        Returns
        -------
        np.ndarray
            The flattop timebase domain for D3D.
        """
        try:
            (
                ip_prog,
                t_ip_prog,
            ) = params.physics_method_params.mds_conn.get_data_with_dims(
                f"ptdata('iptipp', {params.physics_method_params.shot_id})",
                tree_name="d3d",
            )
            t_ip_prog = t_ip_prog / 1.0e3  # [ms] -> [s]
            polarity = np.unique(
                params.physics_method_params.mds_conn.get_data(
                    f"ptdata('iptdirect', {params.physics_method_params.shot_id})",
                    tree_name="d3d",
                )
            )
            if len(polarity) > 1:
                params.logger.info(
                    (
                        "[Shot %s]: Polarity of Ip target is not constant. "
                        "Using value at first timestep."
                    ),
                    params.physics_method_params.shot_id,
                )
                params.logger.debug(
                    "[Shot %s]: Polarity array %s",
                    params.physics_method_params.shot_id,
                    polarity,
                )
                polarity = polarity[0]
            ip_prog = ip_prog * polarity
            dipprog_dt = np.gradient(ip_prog, t_ip_prog)
            ip_prog = interp1(
                t_ip_prog, ip_prog, params.physics_method_params.times, "linear"
            )
            dipprog_dt = interp1(
                t_ip_prog, dipprog_dt, params.physics_method_params.times, "linear"
            )
        except mdsExceptions.MdsException:
            params.logger.warning(
                "[Shot %s]: Could not find flattop timebase. Defaulting to full timebase.",
                params.physics_method_params.shot_id,
            )
            params.logger.debug(
                "[Shot %s]: %s",
                params.physics_method_params.shot_id,
                traceback.format_exc(),
            )
            return None
        epsoff, t_epsoff = params.physics_method_params.mds_conn.get_data_with_dims(
            f"ptdata('epsoff', {params.physics_method_params.shot_id})", tree_name="d3d"
        )
        # [ms] -> [s] # Avoid problem with simultaneity of epsoff being triggered
        # exactly on the last time sample
        t_epsoff = t_epsoff / 1.0e3 + 0.001
        epsoff = interp1(t_epsoff, epsoff, params.physics_method_params.times, "linear")
        railed_indices = np.where(np.abs(epsoff) > 0.5)
        power_supply_railed = np.zeros(len(params.physics_method_params.times))
        power_supply_railed[railed_indices] = 1
        indices_flattop = np.where(
            (np.abs(dipprog_dt) <= 2.0e3)
            & (np.abs(ip_prog) > 100e3)
            & (power_supply_railed != 1)
        )
        return params.physics_method_params.times[indices_flattop]


class RampupAndFlattopDomainSetting(DomainSetting):
    """
    Specifies that the ramp up and flattop domain of the timebase should be used.
    """

    def __init__(self):
        """
        Initialize the RampupAndFlattopDomainSetting with tokamak-specific overrides.
        """
        self.tokamak_overrides = {
            Tokamak.CMOD: self._get_domain_cmod,
        }

    def _get_domain(self, params: DomainSettingParams) -> np.ndarray:
        """
        Get the ramp-up and flattop domain for the given tokamak.

        Parameters
        ----------
        params : DomainSettingParams
            Parameters containing physics method and tokamak information.

        Returns
        -------
        np.ndarray
            The ramp-up and flattop timebase domain.
        """
        raise ValueError(
            f"ramp up and flattop domain not defined for tokamak: {params.tokamak}"
        )

    def _get_domain_cmod(self, params: DomainSettingParams) -> np.ndarray:
        """
        Get the ramp-up and flattop domain for CMOD tokamak.

        Parameters
        ----------
        params : DomainSettingParams
            Parameters containing physics method and tokamak information.

        Returns
        -------
        np.ndarray
            The ramp-up and flattop timebase domain for CMOD.
        """
        ip_parameters = CmodPhysicsMethods.get_ip_parameters(
            params=params.physics_method_params
        )
        ipprog, dipprog_dt = ip_parameters["ip_prog"], ip_parameters["dipprog_dt"]
        indices_flattop_1 = np.where(np.abs(dipprog_dt) <= 6e4)[0]
        indices_flattop_2 = np.where(np.abs(ipprog) > 1.0e5)[0]
        indices_flattop = np.intersect1d(indices_flattop_1, indices_flattop_2)
        if len(indices_flattop) == 0:
            params.logger.warning(
                "[Shot %s]: Could not find flattop timebase. Defaulting to full timebase.",
                params.physics_method_params.shot_id,
            )
            return None
        end_index = np.max(indices_flattop)
        return params.physics_method_params.times[:end_index]


# --8<-- [start:domain_setting_dict]
_domain_setting_mappings: Dict[str, DomainSetting] = {
    "full": FullDomainSetting(),
    "flattop": FlattopDomainSetting(),
    "rampup_and_flattop": RampupAndFlattopDomainSetting(),
}
# --8<-- [end:domain_setting_dict]


def resolve_domain_setting(
    domain_setting: DomainSettingType,
) -> DomainSetting:
    """
    Resolve the given domain setting into a DomainSetting object.

    Parameters
    ----------
    domain_setting : DomainSettingType
        The domain setting to resolve, which can be a string, DomainSetting, or dictionary.

    Returns
    -------
    DomainSetting
        The resolved domain setting.
    """
    if domain_setting is None:
        return FullDomainSetting()

    if isinstance(domain_setting, DomainSetting):
        return domain_setting

    if isinstance(domain_setting, str):
        domain_setting_object = _domain_setting_mappings.get(domain_setting, None)
        if domain_setting_object is not None:
            return domain_setting_object

    if isinstance(domain_setting, dict):
        return DomainSettingDict(domain_setting)

    raise ValueError(f"Invalid domain setting: {domain_setting}")
