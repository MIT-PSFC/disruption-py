from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.settings.set_times_request import SetTimesRequest, SetTimesRequestParams
from disruption_py.settings.shot_settings import ShotSettings

class MagneticsSetTimesRequest(SetTimesRequest):
    """Get times request for using the start and end times of the magnetics tree, with a 
    custom timestep passed in the constructor.
    
    Parameters
    ----------
    timestep : float
        The timestep to use for the magnetics timebase.
    """
    def __init__(self, timestep = 0.004):
        self.timestep = timestep
        
    def _get_times(self, params : SetTimesRequestParams) -> np.ndarray:
        try:
            magnetics_tree = params.tree_manager.open_tree(tree_name='magnetics')
            magnetic_times = magnetics_tree.getNode('\ip').dim_of().data().astype('float64', copy=False)
            # interpolate self._times to be every [timestep] seconds
            return np.arange(magnetic_times[0], magnetic_times[-1], self.timestep).astype('float64', copy=False)
        except:
            params.logger.info("Failed to set up magnetic timebase")
            return None


cmod_handler = CModHandler()
shot_settings = ShotSettings(
    # retrieve existing data from the disruption_warnings table use it to prepopulate queries
    existing_data_request="sql",
    
    # use the efit timebase preset for set_times_request
    # uses the efit timebase when returning data 
    set_times_request="efit",
        
    # only run the get_ip_parameters method
    run_methods=["_get_ip_parameters"],
    run_tags=[],
    run_columns=[],
    # automatically uses the CSVOutputRequest preset because of the .csv file descriptor
    # streams outputted data to the ip_data.csv file
    output_type_request="ip_data.csv", 
)
shot_data = cmod_handler.get_shots_data(
    shot_ids_request="shot_ids_of_interest.txt", # use the shot ids listed in the 
    shot_settings=shot_settings,
    num_processes = 4,
)