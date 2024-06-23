#!/usr/bin/env python3

from abc import ABC
from collections.abc import Iterable
from dataclasses import dataclass
from logging import Logger
from typing import Callable, List

from disruption_py.mdsplus_integration.mds_connection import MDSConnection
from disruption_py.shots.helpers.method_metadata import (
    get_method_metadata,
    is_registered_method,
)
from disruption_py.shots.shot_props import ShotProps
from disruption_py.utils.mappings.tokamak import Tokamak


@dataclass
class ShotDataRequestParams:
    """Params passed by disruption_py to decorated methods.

    Attributes
    ----------
    mds_conn : ShotConnection
        The shot connection object containing the connection to MDSPlus. The same as shot_props.mds_conn.
    shot_props : ShotProps
                A reference to the shot props object containing the setup information, such as the shot id,
        timebase, and disruption time, for the shot data retrieval from MDSPlus.
    tokamak : Tokemak
        The tokamak for which the set times request is made.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """

    mds_conn: MDSConnection
    shot_props: ShotProps
    tokamak: Tokamak
    logger: Logger


class ShotDataRequest(ABC):

    def setup(self, shot_data_request_params: ShotDataRequestParams):
        """Do any setup for the shot such as giving"""
        pass

    def get_request_methods_for_tokamak(self, tokamak: Tokamak) -> List[Callable]:
        """Method used to determine which methods should be considered for execution given the tokamak.

        The default implementation returns all methods that have had the provided tokamak included in the tokamak parameter
        of there `cached_method` or `parameter_cached_method` decorator. This method may be overridden by subclasses if an
        alternate scheme of determining which methods can be executed for each tokamak is required.

        Parameters
        ----------
        tokamak : Tokamak
            The tokamak for which we are retrieving elligible methods.

        Returns
        -------
        List[Callable]
            A list of methods that can be considered for execution for the given tokamak.
        """
        request_methods = []
        for method_name in dir(self):
            method = getattr(self, method_name, None)
            if method is None or not is_registered_method(method):
                continue
            method_metadata = get_method_metadata(method, should_throw=True)
            if (
                method_metadata.tokamaks is None
                or tokamak is method_metadata.tokamaks
                or (
                    isinstance(method_metadata.tokamaks, Iterable)
                    and tokamak in method_metadata.tokamaks
                )
            ):
                request_methods.append(method_name)
        return request_methods
