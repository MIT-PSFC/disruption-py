# import MDSplus
import pandas as pd
import os
import sys

import sys
from datetime import datetime
sys.path.insert(0, "/home/spangher/projects/disruption-py-github/disruption-py/")

from disruption_py.handlers import cmod_handler
from disruption_py.settings.shot_settings import ShotSettings
from disruption_py.settings import shot_id_requests
from disruption_py.settings.log_settings import LogSettings
from disruption_py.settings.output_type_requests import CSVOutputRequest
import logging
import multiprocessing
import pyarrow.parquet as pq

def process_shots(shot_ids_subset, flattop_shots):

    # split the shot_ids_subset into chunks of 10
    shot_ids_subset = list(chunk_list(shot_ids_subset, 10))

    processed_chunks = pd.read_csv("chunks/processed_chunks.csv")

    # for each chunk, run the following
    for chunk in shot_ids_subset:
        formatted_datetime = datetime.now().strftime("%Y%m%d_%H%M%S")
        settings = ShotSettings(
            existing_data_request=None,
            # output_type_request=CSVOutputRequest(filepath=f"chunks/chunked_test_output_{formatted_datetime}.csv", flexible_columns=True),
            efit_tree_name="efitspangher",
            attempt_local_efit_env=[("efitspangher", "/home/spangher/projects/data/CMOD_data_gathering/EFIT_runs/trees")],
            signal_domain="rampup_and_flattop",
            log_settings=LogSettings(console_log_level=logging.DEBUG),
            set_times_request="magnetics004",
            run_columns=features,
            # only_requested_columns=True,
            )
        
        handler = cmod_handler.CModHandler()
        result = handler.get_shots_data(
            shot_id_requests.ListShotIdRequest(chunk),
            shot_settings=settings)

        # remove shot_ids_subset from flattop_shots
        shot_list = flattop_shots.copy()

        # append chunks to a csv called "processed_chunks.csv"
        temp_processed_chunks = pd.DataFrame(chunk, columns=["shot"])
        processed_chunks = pd.concat([processed_chunks, temp_processed_chunks], ignore_index=True, sort=False)

        # remove chunk from shot_list

        shot_list = shot_list[~shot_list["shot"].isin(processed_chunks["shot"])]
        shot_list.to_csv("/home/spangher/projects/disruption-py-github/disruption-py/shots_to_run_disruption_py_on.csv", index=False)

        print("------------------------------------")
        print("------------------------------------")
        print("------------------------------------")
        print("------------------------------------")
        print("------------------------------------")
        print("------------------------------------")
        print("finished processing shots in this one chunk")
        print("------------------------------------")
        print("------------------------------------")
        print("------------------------------------")
        print("------------------------------------")
        print("------------------------------------")
        print("------------------------------------")

    # You might want to return the result or save it in a file

def chunk_list(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


if __name__ == "__main__":
    filepath = "/home/spangher/projects/data/CMOD_data_gathering/EFIT_runs/processed_shots.csv"
    flattop_shots = pd.read_csv("/home/spangher/projects/disruption-py-github/disruption-py/shots_to_run_disruption_py_on.csv")

    # remove nans from flattop_shots
    flattop_shots = flattop_shots.dropna()

    shot_list = flattop_shots["shot"]
    
    # remove any entries from shot_list that are smaller than 100000
    non_efit_derived_features = [
        "ts_width"
    ]
    efit_derived_features = [
        "beta_p", "kappa", "li", "lower_gap", "q95", "zcur",
        "Te_width", "time", "shot", 'btor', 'beta_N', 'beta_p',
        'Wmhd',  "zmag", "qstar",
    ] 
    features = non_efit_derived_features + efit_derived_features

    
   # Split shot_ids into chunks for multiprocessing
    shot_ids = shot_list.values.flatten().tolist()
    num_processes = multiprocessing.cpu_count() - 2 # or set a fixed number
    shot_id_chunks = list(chunk_list(shot_ids, len(shot_ids) // num_processes))

    # Create a pool of processes
    with multiprocessing.Pool(num_processes) as pool:
        pool.starmap(process_shots, [(chunk, flattop_shots) for chunk in shot_id_chunks])
