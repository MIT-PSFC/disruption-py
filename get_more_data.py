"""
This script retrieves shot data from the database and saves it to a CSV file.
It uses the `disruption_py` library to handle the data retrieval process.

Call examples:
    python get_data.py float32
    python get_data.py float64
    python get_data.py None

"""

import os
from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.workflow import get_shots_data
from disruption_py.inout.sql import DummyDatabase
# import numpy as np
# import pickle

timebase = "efit"  # np.arange(0.9, 1., 0.001)
retrieval_settings = RetrievalSettings(
    time_setting=timebase,
)



def main(experimen_option):
    """
    Main function to retrieve shot data from the database.
    Parameters
    ----------
    experimen_option : str
        The precision option for the data retrieval.
        Options are "float32", "float64", or "None".

    Returns
    -------
    None
        The function retrieves the data and saves it to a CSV file.

    Example
    -------
    >>> main("float32")
    """
    if experimen_option == "float32":
        output_setting = "data32_disruption_warning_table.csv"
    elif experimen_option == "float64":
        output_setting = "data64_disruption_warning_table.csv"
    elif experimen_option == "None":
        output_setting = "dataNone_disruption_warning_table.csv"
    else:
        raise ValueError("Invalid option. Choose 'float32', 'float64', or 'None'.")

    # Check if the output file already exists
    if os.path.exists(output_setting):
        print(f"Output file '{output_setting}' already exists. Removing.")
        # remove the file for overwrite it
        os.remove(output_setting)
    else:
        print(
            f"Output file '{output_setting}' does not exist. Proceeding with data retrieval."
        )

    # Check if the output disruption-py/output_cast_float32.txt file exists
    if os.path.exists(f"./output_cast_{experimen_option}.txt"):
        print(
            f"Output file 'disruption-py/output_{experimen_option}.txt' already exists. Removing."
        )
        # remove the file for overwrite it
        os.remove(f"./output_cast_{experimen_option}.txt")
    else:
        print(
            f"Output file './output_cast_{experimen_option}.txt' does not exist. Proceeding with data retrieval."
        )

    
    shot_data = get_shots_data(
        tokamak="cmod",
        shotlist_setting="disruption_warning",
        num_processes=4,
        output_setting=output_setting
    )

    # Save the shot data to a pickle file (optional)
    # with open("shot_data1.pkl", "wb") as f:
    #     pickle.dump(shot_data, f)


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        option = sys.argv[1]
    else:
        option = "float64"  # Default option

    main(option)
