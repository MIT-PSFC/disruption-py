#!/usr/bin/env python3

"""
DIII-D CER retrieval methods.

Provides routines to retrieve CER (Charge Exchange Recombination) data from
DIII-D MDSplus, focusing on the CERAUTO analysis (mode 'a').
"""

import numpy as np
import xarray as xr
from typing import List, Dict, Tuple, Optional

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
from disruption_py.machine.tokamak import Tokamak


# chord lists
TANGENTIAL_CHORDS = [f"t{i}" for i in range(1, 57)]
VERTICAL_CHORDS = [f"v{i}" for i in range(1, 33)]
CER_CHORDS = TANGENTIAL_CHORDS + VERTICAL_CHORDS

# CER data type codes (as described in the DIII-D CER naming conventions)
CER_DATA_TYPES = [
    "r",
    "phi",
    "z",
    "amp",
    "samp",
    "samp_ps",
    "ti",
    "tic",
    "sti",
    "sti_ps",
    "rot",
    "srot",
    "srot_ps",
    "rotc",
    "nz",
    "fz",
    "zeff",
    "loc",
    "sloc",
    "sloc_ps",
    "vb",
    "svb",
    "svb_ps",
    "time",
    "stime",
    "nave",
]


class D3DCerMethods:
    """Class for retrieving CER (Charge Exchange Recombination) data for DIII-D."""

    @staticmethod
    def _build_pointname(dtype: str, chord: str, mode: str = "a") -> str:
        return f"cer{mode}{dtype}{chord}"

    @staticmethod
    def get_cer_node(params: PhysicsMethodParams, dtype: str, chord: str, mode: str = "a"):
        """Retrieve a single CER node using the standard naming convention.

        Returns (data, dims) where dims is a list (possibly empty) of dimension
        arrays. On failure returns (None, None) and logs a warning.
        """
        point = D3DCerMethods._build_pointname(dtype, chord, mode)
        path = f"ptdata('{point}', {params.shot_id})"
        try:
            res = params.mds_conn.get_data_with_dims(path=path, tree_name="d3d")
            if isinstance(res, tuple):
                data = res[0]
                dims = list(res[1:]) if len(res) > 1 else []
            else:
                data = res
                dims = []
            data = np.asarray(data)
            return data, dims
        except Exception as e:  # pragma: no cover - defensive on MDS failures
            params.logger.warning(f"[Shot {params.shot_id}] CER node not available: {point} | {e}")
            return None, None

    @staticmethod
    @physics_method(tokamak=Tokamak.D3D)
    def get_all_cer_data(params: PhysicsMethodParams):
        """Retrieve CER CERAUTO data for all chords.

        Returns an `xarray.Dataset` with dimensions (`chord`, `time`) and
        chord spatial coordinates (`r`, `phi`, `z`). If no CER data is found
        for the shot returns `None`.
        """

        chords = CER_CHORDS
        n_chords = len(chords)
        times_target = np.asarray(params.times)
        if times_target.size == 0:
            params.logger.warning(f"[Shot {params.shot_id}] No target times available in params.times")
            return None

        # Pre-fetch per-chord time arrays (if available) and convert ms->s when needed
        per_chord_times: Dict[str, Optional[np.ndarray]] = {}
        for chord in chords:
            data, dims = D3DCerMethods.get_cer_node(params, "time", chord)
            if data is not None:
                try:
                    times_ch = np.asarray(data, dtype=float)
                    if np.nanmax(times_ch) > 1e3:
                        times_ch = times_ch / 1000.0
                    per_chord_times[chord] = times_ch
                except Exception:
                    params.logger.warning(f"[Shot {params.shot_id}] Failed to parse CER time for {chord}")
                    per_chord_times[chord] = None
            else:
                per_chord_times[chord] = None

        # spatial coords (per chord)
        r = np.full(n_chords, np.nan, dtype=np.float32)
        phi = np.full(n_chords, np.nan, dtype=np.float32)
        z = np.full(n_chords, np.nan, dtype=np.float32)

        data_vars: Dict[str, Tuple[Tuple[str, str], np.ndarray]] = {}
        found_any = False

        # Collect time-series data types
        for dtype in CER_DATA_TYPES:
            if dtype in ("r", "phi", "z", "time"):
                continue

            arr = np.full((n_chords, times_target.size), np.nan, dtype=np.float32)
            any_for_dtype = False

            for idx, chord in enumerate(chords):
                node_data, dims = D3DCerMethods.get_cer_node(params, dtype, chord)
                if node_data is None:
                    continue
                node_data = np.asarray(node_data, dtype=float)

                # determine time axis for this node
                times_ch = per_chord_times.get(chord)
                if times_ch is None and dims and len(dims) >= 1 and dims[0] is not None:
                    try:
                        times_ch = np.asarray(dims[0], dtype=float)
                        if np.nanmax(times_ch) > 1e3:
                            times_ch = times_ch / 1000.0
                    except Exception:
                        times_ch = None

                try:
                    if times_ch is None:
                        # if shape matches target times assign directly
                        if node_data.size == times_target.size:
                            arr[idx, :] = node_data.astype(np.float32)
                            any_for_dtype = True
                        else:
                            params.logger.warning(
                                f"[Shot {params.shot_id}] Cannot determine time axis for CER {dtype}{chord}; skipping"
                            )
                            continue
                    else:
                        if node_data.ndim == 0 or node_data.size == 1:
                            arr[idx, :] = float(node_data)
                        else:
                            try:
                                arr[idx, :] = interp1(times_ch, node_data, times_target, axis=-1).astype(np.float32)
                            except Exception:
                                try:
                                    arr[idx, :] = np.interp(times_target, times_ch, node_data).astype(np.float32)
                                except Exception:
                                    params.logger.warning(f"[Shot {params.shot_id}] Failed to interpolate CER {dtype}{chord}")
                                    continue
                        any_for_dtype = True
                except Exception as e:  # pragma: no cover - defensive
                    params.logger.warning(f"[Shot {params.shot_id}] Error processing CER {dtype}{chord}: {e}")
                    continue

            if any_for_dtype:
                found_any = True
                data_vars[dtype] = (("chord", "time"), arr)

        # fetch spatial coords r, phi, z
        for idx, chord in enumerate(chords):
            d, _ = D3DCerMethods.get_cer_node(params, "r", chord)
            if d is not None:
                try:
                    r[idx] = float(np.asarray(d).flatten()[0])
                except Exception:
                    r[idx] = np.nan
            d, _ = D3DCerMethods.get_cer_node(params, "phi", chord)
            if d is not None:
                try:
                    phi[idx] = float(np.asarray(d).flatten()[0])
                except Exception:
                    phi[idx] = np.nan
            d, _ = D3DCerMethods.get_cer_node(params, "z", chord)
            if d is not None:
                try:
                    z[idx] = float(np.asarray(d).flatten()[0])
                except Exception:
                    z[idx] = np.nan

        if not found_any:
            params.logger.warning(f"[Shot {params.shot_id}] No CER CERAUTO data found.")
            return None

        coords = {
            "chord": ("chord", np.array(chords, dtype=object)),
            "time": ("time", times_target),
            "r": ("chord", r),
            "phi": ("chord", phi),
            "z": ("chord", z),
        }

        ds = xr.Dataset(data_vars=data_vars, coords=coords)
        ds.attrs["shot"] = int(params.shot_id)

        params.logger.info(f"[Shot {params.shot_id}] Retrieved CER CERAUTO data: {list(data_vars.keys())}")
        return ds
