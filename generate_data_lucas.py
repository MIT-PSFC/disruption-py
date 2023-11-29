import gc
import MDSplus
from scipy.interpolate import interp1d
from disruption_py.shots.cmod_shot import CModShot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import logging
import sys
from concurrent.futures import ThreadPoolExecutor
import pdb
import os

VAR_TYPE = "nefit"

def process_batch(batch_shots, features, mag_interp, batch_num, efit):
    batch_data = []

    for shot in batch_shots:
        try:
            df = pull_data(shot, mag_interp=mag_interp)
            filtered_df = df.reindex(columns=features)
            batch_data.append(filtered_df)
        except Exception as e:
            logging.warning(f"Error processing shot {shot}: {e}")
            continue

    try:
        batch_df = pd.concat(batch_data, ignore_index=True)
    except Exception as e:
        logging.warning(f"Error concatenating batch data: {e}")
    
    remove_shots_from_csv(batch_shots)
        
    return batch_df


def get_last_number(basefilename):
    """Get the last used number based on existing filenames."""
    existing_files = [f for f in os.listdir() if f.startswith(basefilename)]
    existing_numbers = [int(f.replace(basefilename, '').replace('.csv', '')) for f in existing_files if f.replace(basefilename, '').replace('.csv', '').isdigit()]
    
    if existing_numbers:
        return max(existing_numbers)
    else:
        return -1  # If no files exist, start from 0

def get_next_number(batch_num, basefilename):
    """Calculate the next number, considering the last used number in existing files."""
    last_number = get_last_number(basefilename)
    return last_number + 1 + batch_num  # Increment by one and add current batch number


# Set up logging
logging.basicConfig(filename="pull_data.log",
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    level=logging.INFO)


# Define a custom handler to catch errors
class TerminateAfterThresholdHandler(logging.Handler):
    def __init__(self, threshold=500, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.count = 0
        self.threshold = threshold

    def emit(self, record):
        if "Failed to open efit tree analysis" in record.getMessage():
            self.count += 1
            if self.count >= self.threshold:
                print("Error threshold reached. Exiting...")
                sys.exit(1)

# Attach the custom handler to the root logger
root_logger = logging.getLogger()
handler = TerminateAfterThresholdHandler()
root_logger.addHandler(handler)


def pull_data(shot, mag_interp):
    cmod_shot = CModShot(shot_id=shot, mag_interp=mag_interp, interp_scheme="previous")
    return cmod_shot.data


def remove_shots_from_csv(pulled_shots):
    """Remove the shots that have been tried from the list of shots to try."""

    # print pulled_shots aggregated into one string
    print("--" * 20)
    print("pulled_shots eliminated from remaining_shot.csv:")
    print(pulled_shots)
    print("--" * 20)

    current_shots = pd.read_csv(data_dir + "/" + f'remaining_{VAR_TYPE}_shots.csv')["shot"]
    remaining_shots = current_shots[~current_shots.isin(pulled_shots)]
    remaining_shots.to_frame(name='shot').to_csv(data_dir + "/" + f'remaining_{VAR_TYPE}_shots.csv', index=False)  
       
    return 


if __name__ == "__main__":

    data_dir = "/home/spangher/projects/disruption-py-github/disruption-py"
    shots = pd.read_csv(data_dir + "/" + f'remaining_{VAR_TYPE}_shots.csv')
    shots = shots["shot"]

    ########
    # EFIT Pulling

    batch_size = 5

    num_batches = len(shots) // batch_size + (1 if len(shots) % batch_size != 0 else 0)
    
    # Define your features for each category
    non_efit_derived_features = [
        "Greenwald_fraction", "n_equal_1_normalized", "v_loop",
        "radiated_fraction", "time", "shot", 'p_lh', 'p_icrf',
        "ip", "ip_error", "ip_prog", "dip_dt", "p_rad", "time_until_disrupt",
    ]
    efit_derived_features = [
        "beta_p", "kappa", "li", "lower_gap", "q95", "zcur",
        "Te_width", "time", "shot", 'btor', 'beta_N', 'beta_p',
        'Wmhd',  "zmag", "qstar",
    ] 

    # "dli_dt", "dWmhd_dt",

    if VAR_TYPE == "nefit":
        feature_group = {"features": non_efit_derived_features,
                         "mag_interp": True, "filename": "shots_9_22_nefit_",
                         "efit": "nefit"}
    else:
        feature_group = {"features": efit_derived_features,
                         "mag_interp": False,
                         "filename": "shots_9_22_efit_",
                         "efit": "efit"}

    with ThreadPoolExecutor(max_workers=12) as executor:
        logging.info(f"Processing feature group {feature_group['filename']}")
        futures = []
        for batch_num in range(num_batches):
            logging.info(f"Processing batch {batch_num}...")
            start_idx = batch_num * batch_size
            end_idx = start_idx + batch_size
            batch_shots = shots[start_idx:end_idx]

            # Schedule the batch processing to run
            future = executor.submit(process_batch, batch_shots, feature_group["features"],
                                        feature_group["mag_interp"], batch_num, feature_group["efit"])
            futures.append(future)
        
        logging.info("Waiting for batches to finish...")

        for f, batch_num in zip(futures, range(num_batches)):
            logging.info(f)
            batch_df = f.result()

            # test if there's more than 10 rowsin the batch_df
            if batch_df.shape[0] > 10:
                batch_df.to_csv(f"{feature_group['filename']}{get_next_number(batch_num=batch_num, basefilename=feature_group['filename'])}.csv", index=False, header=True)

            # Collect all shot numbers from the batch_df
            pulled_shots = batch_df['shot'].unique().tolist()

            # Save the reduced shot list back to 'remaining_{VAR_TYPE}_shots.csv'
            print("--" * 20)
            print("Saving remaining shots to CSV file...")
            remove_shots_from_csv(pulled_shots=pulled_shots)
            print("--" * 20)

    logging.info("All batches processed!")
