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
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
            ts = 0
        if ts == 0:
            params.logger.info("[Shot %s]: Both TS data missing", params.shot_id)
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

    @staticmethod
    @physics_method(
        columns=["n_equal_1_normalized", "n_equal_1_mode"],
        tokamak=Tokamak.D3D,
    )
    def get_n1_bradial_parameters(params: PhysicsMethodParams):
        """
        Retrieve normalized n=1 bradial parameters for a given shot.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information.

        Returns
        -------
        dict
            A dictionary containing the following keys:
            - 'n_equal_1_normalized' : array
                Normalized n=1 bradial parameter values.
            - 'n_equal_1_mode' : array
                Raw n=1 bradial parameter values.
        """
        output = {
            "n_equal_1_normalized": [np.nan],
            "n_equal_1_mode": [np.nan],
        }
        # The following shots are missing bradial calculations in MDSplus and
        # must be loaded from a separate datafile
        # TODO: This implementation needs revising.
        if 176030 <= params.shot_id <= 176912:
            raise NotImplementedError
            # TODO: Move to a folder like "/fusion/projects/disruption_warning/data"
            filename = "/fusion/projects/disruption_warning/matlab_programs/recalc.nc"
            ncid = nc.Dataset(filename, "r")
            brad = ncid.variables["dusbradial_calculated"][:]
            t_n1 = ncid.variables["times"][:] * 1.0e-3  # [ms] -> [s]
            shots = ncid.variables["shots"][:]
            shot_indices = np.where(shots == params.shot_id)
            if len(shot_indices) == 1:
                dusbradial = brad[shot_indices, :] * 1.0e-4  # [T]
            else:
                params.logger.info(
                    f"Shot {params.shot_id} not found in {filename}.  Returning NaN."
                )
                dusbradial = np.full(len(params.times), np.nan)
            ncid.close()
        # Check ONFR than DUD(legacy)
        else:
            try:
                # TODO: TREE NAME?
                dusbradial, t_n1 = params.mds_conn.get_data_with_dims(
                    f"ptdata('onsbradial', {params.shot_id})",
                    tree_name="d3d",
                )
                dusbradial = interp1(t_n1, dusbradial, params.times)
                dusbradial *= 1.0e-4  # [T]
            except mdsExceptions.MdsException:
                params.logger.debug(
                    "[Shot %s]: %s", params.shot_id, traceback.format_exc()
                )
                try:
                    dusbradial, t_n1 = params.mds_conn.get_data_with_dims(
                        f"ptdata('dusbradial', {params.shot_id})",
                        tree_name="d3d",
                    )
                    dusbradial = interp1(t_n1, dusbradial, params.times)
                    dusbradial *= 1.0e-4  # [T]
                except mdsExceptions.MdsException:
                    params.logger.info(
                        "[Shot %s]: Failed to get n1 bradial signal. Returning NaN.",
                        params.shot_id,
                    )
                    params.logger.debug(
                        "[Shot %s]: %s", params.shot_id, traceback.format_exc()
                    )
                    return output
        n_equal_1_mode = interp1(dusbradial, t_n1, params.times)
        # Get toroidal field Btor
        b_tor, t_b_tor = params.mds_conn.get_data_with_dims(
            f"ptdata('bt', {params.shot_id})", tree_name="d3d"
        )
        b_tor = interp1(t_b_tor, b_tor, params.times)  # [T]
        n_equal_1_normalized = n_equal_1_mode / b_tor
        output = {
            "n_equal_1_normalized": n_equal_1_normalized,
            "n_equal_1_mode": n_equal_1_mode,
        }
        return output
