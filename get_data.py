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

# list_of_shots = [1100204004]
list_of_shots = [
    950112023,
    950112027,
    950208012,
    950208015,
    950208017,
    950208024,
    950210024,
    950210026,
    950328011,
    950329013,
    950329014,
    950329017,
    950329019,
    950331038,
    950523007,
    950523017,
    950602008,
    950602009,
    950602015,
    950602023,
    950602024,
    951129009,
    951129010,
    951130011,
    951130014,
    951201025,
    951205007,
    951212025,
    960112003,
    960112031,
    960112039,
    960117002,
    960117003,
    960117005,
    960118001,
    960126003,
    960126004,
    960126005,
    960126006,
    960126017,
    960126028,
    960131004,
    960131010,
    990225009,
    990225017,
    990225018,
    990225022,
    1000511004,
    1000511005,
    1000511006,
    1000511012,
    1000511013,
    1000511018,
    1000511019,
    1000511026,
    1000613022,
    1031204013,
    1031204016,
    1040106010,
    1040219006,
    1040219024,
    1040220006,
    1040220008,
    1040325027,
    1040325030,
    1040507011,
    1050621015,
    1050623008,
    1050623010,
    1070531007,
    1070531012,
    1070531014,
    1070531018,
    1070605011,
    1070726018,
    1070807008,
    1070807016,
    1070829008,
    1070829028,
    1070829029,
    1070830004,
    1070830005,
    1070830017,
    1080415017,
    1080416006,
    1080416007,
    1080523018,
    1080523020,
    1090911011,
    1090914006,
    1090914015,
    1091014007,
    1091014015,
    1091014027,
    1091016023,
    1091016029,
    1091016032,
    1091016033,
    1091021008,
    1091021014,
    1091021020,
    1091021032,
    1091021035,
    1091112004,
    1091112005,
    1091203012,
    1091203020,
    1091215012,
    1091215023,
    1100121018,
    1100204004,
    1100204015,
    1100204016,
    1100204017,
    1100204018,
    1100204019,
    1100204021,
    1100204022,
    1100210012,
    1100311003,
    1100311004,
    1100311014,
    1100817008,
    1100817009,
    1100817012,
    1100817014,
    1100817016,
    1100817017,
    1100817019,
    1100817020,
    1100910021,
    1100910022,
    1100910023,
    1101209010,
    1101209011,
    1101209012,
    1101209013,
    1101209014,
    1101209019,
    1101209021,
    1101209024,
    1101209029,
    1101209031,
    1110106010,
    1110106023,
    1110215005,
    1110215010,
    1110215012,
    1110215014,
    1110215016,
    1110215017,
    1110215021,
    1110215031,
    1110309007,
    1110309013,
    1110309022,
    1110309024,
    1110309025,
    1110309030,
    1110309032,
    1110310022,
    1110316005,
    1110316020,
    1110317015,
    1110317026,
    1120824016,
    1120829010,
    1120829014,
    1120829016,
    1120829017,
    1120829022,
    1120829025,
    1120829026,
    1120907008,
    1120907009,
    1120907010,
    1120907016,
    1120907025,
    1120907028,
    1120907032,
    1120907035,
    1120920011,
    1120920017,
    1120920026,
    1140328017,
    1140328029,
    1140522016,
    1140618013,
    1140625012,
    1140625015,
    1140625018,
    1150616024,
    1150826008,
    1150826029,
    1150928016,
    1150929010,
    1160616018,
    1160803022,
    1160803023,
    1160826012,
    1160826014,
    1160826027,
    1160826029,
    1160826030,
    1160909009,
    1160909010,
    1160909011,
    1160921011,
    1160921019,
    1160921020,
    1160921021,
    1160923003,
]


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
        output_setting = "data32.csv"
    elif experimen_option == "float64":
        output_setting = "data64.csv"
    elif experimen_option == "None":
        output_setting = "dataNone.csv"
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
        shotlist_setting=list_of_shots,
        retrieval_settings=retrieval_settings,
        num_processes=4,
        output_setting=output_setting,
        database_initializer=DummyDatabase.initializer,
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
