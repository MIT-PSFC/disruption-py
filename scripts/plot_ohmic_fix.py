"""Grab """

import argparse
import os

from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import pandas as pd


def read_dpy_data(shot, local_data_dir='p_oh_cache', efit_name='efit21', bypass=False, save_data=True):
    data_path = os.path.join(local_data_dir, f"dpy_{efit_name}.csv")
    print(data_path)

    # check if requested data already stored locally
    if os.path.exists(data_path) and not bypass:
        print('Local dpy data already exists')
        data = pd.read_csv(data_path)
        return data
    retrieval_settings = RetrievalSettings(
        efit_nickname_setting="efit21",
        # method/column selection
        # default None: all methods/columns
        run_columns=["ip", "p_rad", "wmhd", "p_icrf", "p_oh", "p_oh_fix", "v_inductive_fix", "inductance", "inductance_smooth", "li"],
        only_requested_columns=True,
        custom_physics_methods=[],
    )

    data = get_shots_data(
        # required argument
        shotlist_setting=[shot],
        # default None: detect from environment
        tokamak=None,
        # default None: standard backend connection(s)
        database_initializer=None,
        connection_initializer=None,
        retrieval_settings=retrieval_settings,
        output_setting="dataset",
        num_processes=1,
        log_settings=LogSettings(
            # default None: "output.log" in temporary session folder
            file_path=None,
            file_level="DEBUG",
            # default None: read from the configuration, "INFO" out of the box
            console_level=None,
        ),
    )
    data = data.to_dataframe()
    if save_data:
        print(f"Saving data to {data_path}")
        data.to_csv(data_path)
    return data

def plot_ohmic_power(df):
    fig, axs = plt.subplots(4,1, sharex=True)
    df['p_net'] = df['p_icrf'] + df['p_oh'] - df['p_rad']
    df['p_net_fix'] = df['p_icrf'] + df['p_oh_fix'] - df['p_rad']
    axs[0].plot(df['time'], df['ip']/1e6, c='k')
    axs[1].plot(df['time'], df['wmhd']/1e3, c='k', marker='.')
    axs[2].plot(df['time'], df['p_rad']/1e6, c='k', marker='.')
    axs[2].plot(df['time'], df['p_icrf']/1e6, c='grey', marker='.')
    axs[2].plot(df['time'], df['p_oh']/1e6, c='b', marker='.')
    axs[2].plot(df['time'], df['p_oh_fix']/1e6, c='r', marker='.')
    axs[2].set_ylim(0,10)
    axs[3].plot(df['time'], df['p_net']/1e6, c='b', marker='o', linestyle='-')
    axs[3].set_ylim(-5,5)
    axs[3].plot(df['time'], df['p_net_fix']/1e6, c='r', marker='.', linestyle='--')
    axs[-1].set_xlabel('Time [s]')
    plt.show()

if __name__=='__main__':
    parser = argparse.ArgumentParser(
        description="Calculate Dpy Ohmic data for a shot and plot"
    )
    parser.add_argument(
        "--shot", "-s", required=True,
        help="Shot number to get data and plot. Only one shot allowed"
    )
    args = parser.parse_args()
    df = read_dpy_data(args.shot, bypass=True)
    #plot_ohmic_power(df)