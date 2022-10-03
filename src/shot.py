from MDSplus import *
from pandas import pd
import numpy as np 

DEFAULT_SHOT_COLUMNS = ['time','shot','time_until_disrupt','ip']

class Shot:
    """
    Base shot class to represent a single shot of a fusion experiment.

    Attributes
    ---
    metadata->dict: 
        labels: {} /* Dictionary where key is the label name(e.g. "disruption period")
                    and value is a list of starting and ending time tuples(e.g.[(10,100)]*/
        commit_hash: // Commit hash of code version that was used to generate shot data 
        timestep: // Time delta. Every array must conform to this time step
        duration: // Length in milliseconds of shot
        description: // Optional text summary to add to shot
    data->pandas.DataFrame:
    """

    # TODO: Populate metadata dict
    def __init__(self, mdsplus_name,shot_id,data_columns = DEFAULT_SHOT_COLUMNS,data = None):
        self._tree = Tree(mdsplus_name,shot_id)
        self._shot_id = shot_id
        self._metadata = {
            'labels': {}
            'commit_hash': {}
            'timestep': {}
            'duration': {}
            'description': ""
            'disrupted': 100  #TODO: Fix
        }
        try:
            self._times = self.tree.getNode(r"\efit18::efit_aeqdsk:time")
        except Exception as e:
            print(e)
            print("WARNING: No EFIT18 data for this shot.")
        self._populate_shot_data()
    
    def _populate_shot_data(self):
        self. = self.get_times
        self.data['time_until_disrupt'] = self.get_time_until_disrupt(self, times) 

    def _calc_time_until_disrupt(self):
        if self.metadata['disrupted']:
            if self.metadata['disrupted'] == -1:
                raise ValueError("Shot disrupted but no t_disrupt found.")   
            time_until_disrupt = times - self.metadata['disrupted'] 
        return time_until_disrupt

    def _calc_Ip_parameters(self, times):
        ip = self._tree.getNode(r"\ip")
        magtime = len(ip)
        dip_dt = np.gradient(ip, magtime)

    def _calc_Z_parameters(self, times):

    def _calc_power(self, times):

    def _calc_EFIT_parameters(self, times):
    
    def _calc_kappa_area(self, times):

    def _calc_rotation_velocity(self, times):
    
    def _calc_n_equal_1_amplitude(self, times):
         
    def _calc_TS_data_cmod(self, times):

    def _calc_densities(self, times):

    def _calc_efc_current(self, times):

    def _calc_sxr_data(self, times):


        


class D3DShot(Shot):
    pass

class EASTShot(Shot):
    pass

class CMODShot(Shot):
    pass 
    