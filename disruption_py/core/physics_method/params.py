#!/usr/bin/env python3

"""
Module for defining parameters used in physics methods for DisruptionPy.
"""
from dataclasses import dataclass, field
from typing import Any, Dict

import numpy as np
import xarray as xr
from loguru import logger

from disruption_py.core.utils.misc import shot_msg_patch, to_tuple
from disruption_py.inout.mds import MDSConnection
from disruption_py.machine.tokamak import Tokamak


@dataclass
class PhysicsMethodParams:
    """
    Holder for useful variables for the physics methods like an MDSplus connection
    and the timebase for the data.
    """

    shot_id: int
    tokamak: Tokamak
    disruption_time: float
    mds_conn: MDSConnection
    times: np.ndarray

    def __post_init__(self):
        self.logger = shot_msg_patch(logger, self.shot_id)

    cached_results: Dict[str, Any] = field(default_factory=dict)

    @property
    def disrupted(self) -> bool:
        """
        Check if the disruption time is set.

        Returns
        -------
        bool
            True if disruption time is not None, False otherwise.
        """
        return self.disruption_time is not None

    def cleanup(self) -> None:
        """
        Clean up resources used by the physics method parameters.
        """
        self.mds_conn.cleanup()
        self.times = None
        self.cached_results.clear()

    def to_dataset(self, data: Dict[str, np.ndarray]) -> xr.Dataset:
        """
        Create a dataset of variables, including coordinates based on the parameters.

        Returns
        -------
        xr.Dataset
            A dataset of assumed dimension `idx`, with `shot` and `time` as coordinates.
        """
        data_vars = to_tuple(data=data, dim="idx")
        coords = to_tuple(
            data={
                "shot": len(self.times) * [self.shot_id],
                "time": self.times,
            },
            dim="idx",
        )
        return xr.Dataset(data_vars=data_vars, coords=coords)
