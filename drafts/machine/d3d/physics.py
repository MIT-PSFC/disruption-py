#!/usr/bin/env python3

import traceback

import numpy as np

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.machine.tokamak import Tokamak


class D3DDraftPhysicsMethods:

    @staticmethod
    @physics_method(
        tags=["experimental"],
        tokamak=Tokamak.D3D,
        columns=[
            "te_core",
            "ne_core",
            "te_core",
            "ne_edge",
            "te_edge_80to85",
            "ne_edge_80to85",
            "te_edge_85to90",
            "ne_edge_85to90",
            "te_edge_90to95",
            "ne_edge_90to95",
            "te_edge_95to100",
            "ne_edge_95to100",
            "te_sep",
            "ne_sep",
        ],
    )
    def get_core_edge_vals(params: PhysicsMethodParams):
        ##################################################
        # Settings
        ts_data_type = "blessed"  # either 'blessed', 'unblessed', or 'ptdata'
        # metric to use for core/edge binning (either 'psin' or 'rhovn')
        ts_radius = "rhovn"
        # ts_radius value defining boundary of 'core' region (between 0 and 1)
        ts_core_margin = 0.3
        # ts_radius value defining inner and outer side of 'edge' region
        # (between ts_core_margin and 1)
        ts_edge_inner = 0.85
        ts_edge_outer = 0.95
        # All data outside this range excluded. For example, psin=0 at magnetic axis
        # and 1 at separatrix.
        ts_radial_range = (0, 1)
        # set to true to interpolate ts_channel data onto equispaced radial grid
        ts_equispaced = True
        ###################################################

        # Initialize arrays
        te_core = np.full(len(params.times), np.nan)
        ne_core = np.full(len(params.times), np.nan)
        # Averaged over edge region
        te_edge = np.full(len(params.times), np.nan)
        ne_edge = np.full(len(params.times), np.nan)
        # Averaged over 85th to 88th surface
        te_edge_80to85 = np.full(len(params.times), np.nan)
        ne_edge_80to85 = np.full(len(params.times), np.nan)
        te_edge_85to90 = np.full(len(params.times), np.nan)
        ne_edge_85to90 = np.full(len(params.times), np.nan)
        te_edge_90to95 = np.full(len(params.times), np.nan)
        ne_edge_90to95 = np.full(len(params.times), np.nan)
        te_edge_95to100 = np.full(len(params.times), np.nan)
        ne_edge_95to100 = np.full(len(params.times), np.nan)
        # Separatrix
        te_sep = np.full(len(params.times), np.nan)
        ne_sep = np.full(len(params.times), np.nan)

        # Try to get data via _get_ne_te()
        try:
            ts = D3DPhysicsMethods._get_ne_te(params)
            efit_dict = D3DPhysicsMethods._get_efit_dict(params)
            ts["psin"], ts["rhovn"] = D3DPhysicsMethods.efit_rz_interp(ts, efit_dict)
        except Exception as e:
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
            ts = 0
        if ts == 0:
            params.logger.info(
                f"[Shot {params.shot_id}]:Both TS data missing for shot #{params.shot_id}"
            )
        if ts != 0:
            # Drop data outside of valid range #ADM: this looks unfinished
            invalid_indices = np.where(
                (ts[ts_radius] < ts_radial_range[0])
                | (ts[ts_radius] > ts_radial_range[1])
            )

        # TODO: 1) Interpolate in core and edge regions, 2) compute average in
        # these regions and store in respective array. Note that we may need to
        # expand the available indices beyond 1

        return {
            "te_core": te_core,
            "ne_core": ne_core,
            "te_core": te_edge,
            "ne_edge": ne_edge,
            "te_edge_80to85": te_edge_80to85,
            "ne_edge_80to85": ne_edge_80to85,
            "te_edge_85to90": te_edge_85to90,
            "ne_edge_85to90": ne_edge_85to90,
            "te_edge_90to95": te_edge_90to95,
            "ne_edge_90to95": ne_edge_90to95,
            "te_edge_95to100": te_edge_95to100,
            "ne_edge_95to100": ne_edge_95to100,
            "te_sep": te_sep,
            "ne_sep": ne_sep,
        }
