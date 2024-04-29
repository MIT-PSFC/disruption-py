from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.settings.shot_settings import ShotSettings
import numpy as np
import matplotlib.pyplot as plt

#shot_list = [1150805012, 1150805013, 1150805014]
shot_list = [1120119001, 1070406023, \
             1100302003, 1100302005, 1160511013,]
#shot_list = [1100302001]

cmod_handler = CModHandler()
shot_settings = ShotSettings(
    # uses the efit timebase when returning data 
    set_times_request="efit",
     efit_tree_name="efit18",
    
    # run all available methods
    run_tags=[],
    run_columns=["zcur","lower_gap","upper_gap","z_error","ip"],
)
shot_data = cmod_handler.get_shots_data(
    # Retrieve data for the desired shots
    shot_ids_request=shot_list,
    shot_settings=shot_settings,
    
    # automatically stream retrieved data to a csv file by passing in a file path ending in .csv
    output_type_request="dataframe",#"data.csv",
    
    num_processes = 1
)

print(shot_data.head())
print(np.shape(shot_data["zcur"]))

for shot in shot_list:
    shot_indx = np.where(shot_data["shot"] == shot)[0]

    fig1,ax1 = plt.subplots(3,1,sharex=True)
    ax1[0].set_title(str(shot))
    ax1[0].plot(np.array(shot_data["time"][shot_indx]), np.array(shot_data["ip"][shot_indx])*1e-6,label='Ip')
    ax1[1].plot(np.array(shot_data["time"][shot_indx]), np.array(shot_data["zcur"][shot_indx]),label='zcur')
    ax1[1].plot(np.array(shot_data["time"][shot_indx]), np.array(shot_data["z_prog"][shot_indx]),linestyle='--',label='zprog')
    ax1[1].plot(np.array(shot_data["time"][shot_indx]), np.array(shot_data["z_error"][shot_indx]),linestyle='dotted',label='z_error')
    ax1[2].plot(np.array(shot_data["time"][shot_indx]), np.array(shot_data["upper_gap"][shot_indx]),label='upper gap')
    ax1[2].plot(np.array(shot_data["time"][shot_indx]), np.array(shot_data["lower_gap"][shot_indx]),label='lower gap')
    
    ax1[0].legend()
    ax1[1].legend()
    ax1[2].legend()

    ax1[1].set_ylabel("[MA]")
    ax1[1].set_ylabel("[m]")
    ax1[2].set_ylabel("[m]")
    ax1[-1].set_xlabel("time [s]")

plt.show()