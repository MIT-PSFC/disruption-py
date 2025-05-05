#!/usr/bin/env python3

"""
Module for defining parameters used in physics methods for DisruptionPy.
"""
from dataclasses import dataclass, field
from typing import Any, Dict, Tuple

import numpy as np
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

    def to_coords(self) -> Dict[str, Tuple[str, np.ndarray]]:
        """
        Create a dictionary of coordinates based on the parameters.

        Returns
        -------
        Dict[str, Tuple[str, np.ndarray]]
            A dictionary with `shot` and `time` as coordinates for dimension `idx`.
        """
        return to_tuple(
            data={
                "shot": len(self.times) * [self.shot_id],
                "time": self.times,
            },
            dim="idx",
        )
