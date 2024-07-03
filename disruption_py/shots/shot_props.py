#!/usr/bin/env python3

import logging
from dataclasses import dataclass, field
from typing import Any, Dict

import numpy as np
import pandas as pd

from disruption_py.io.mds import MDSConnection
from disruption_py.machine.tokamak import Tokamak


@dataclass
class ShotProps:
    logger = logging.getLogger("disruption_py")

    shot_id: int
    tokamak: Tokamak
    disruption_time: float
    mds_conn: MDSConnection
    times: np.ndarray
    input_data: pd.DataFrame  # input data passed to shot class
    pre_filled_shot_data: pd.DataFrame  # input data after changed to times domain
    interpolation_method: Any  # Fix
    metadata: dict

    _cached_results: Dict[str, Any] = field(default_factory=dict)

    @property
    def disrupted(self):
        return self.disruption_time is not None

    def cleanup(self):
        self.mds_conn.cleanup()
        self.times = None
        self._cached_results.clear()
